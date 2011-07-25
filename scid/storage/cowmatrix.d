/** Implementation of the CowMatrix container, a copy-on-write matrix which is the default container for matrix storage
    types.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.storage.cowmatrix;

import scid.matrix;

import scid.ops.common;
import scid.storage.arraydata;
import scid.storage.cowarray;
import scid.common.meta;

import scid.internal.assertmessages;
import scid.bindings.blas.dblas;
import std.algorithm, std.array, std.typecons;

/** A copy-on-write matrix. Used as a container for storage types. */
struct CowMatrix( ElementType_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor ) {
	alias ElementType_                                                 ElementType;
	alias ArrayData!ElementType                                        Data;
	alias storageOrder_                                                storageOrder;
	alias CowMatrix!(ElementType, transposeStorageOrder!storageOrder ) Transposed;
	alias CowArrayRef!ElementType                                      ArrayType;
	
	enum isRowMajor = (storageOrder == StorageOrder.RowMajor);
	
	/** Allocate a new matrix of given dimensions. Initialize with zero. */
	this( size_t newRows, size_t newCols )
	in {
		assert( newRows != 0 && newCols != 0, zeroDimMsg_( newRows, newCols ) );
	} body {
		this( newRows, newCols, null );
		scal( newRows * newCols, Zero!ElementType, ptr_, 1 );
	}
	
	/** Allocate a new uninitialized matrix of given dimensions. */
	this( size_t newRows, size_t newCols, void* )
	in {
		assert( newRows != 0 && newCols != 0, zeroDimMsg_( newRows, newCols ) );
	} body {
		data_.reset( newRows * newCols );
		ptr_     = data_.ptr;
		rows_    = newRows;
		cols_    = newCols;
		leading_ = minor_;
	}
	
	/** Create a new matrix with given data, data pointer and dimensions. */
	this(  Data data, ElementType* ptr, size_t numRows, size_t numColumns, size_t leadingDimension  )
	in {
		assert( data.owns( ptr ), "Pointer passed to ctor is not owned by data." );
	} body {
		data_    = data;
		ptr_     = ptr;
		rows_    = numRows;
		cols_    = numColumns;
		leading_ = leadingDimension;
		
		assert( leading_ >= minor_, "Leading dimension is smaller than minor dimension." );
		assert( data.owns( ptr_ + leading_ * major_ - 1 ), "Size passed to ctor exceeds size of data." );
	}
	
	/** Create a new matrix with a given major dimension (i.e number of columns for ColumnMajor matrices)
	    and an array with the elements in minor order.
	*/
	this( size_t newMajor, ElementType[] initializer )
	in {
		assert( newMajor != 0 && initializer.length % newMajor == 0, initMsg_( newMajor, initializer ) );
	} body {
		data_.reset( initializer );
		ptr_     = data_.ptr;
		major_   = newMajor;
		minor_   = initializer.length / major_;
		leading_ = minor_;
	}
	
	/** Create a matrix as a copy of another exisiting one. */
	this( typeof(this)* other ) {
		data_    = other.data_;
		ptr_     = other.ptr_;
		rows_    = other.rows_;
		cols_    = other.cols_;
		leading_ = other.leading_;
	}
	
	/** Create a matrix as a slice of an existing one. */
	this( typeof(this)* other, size_t rowStart, size_t numRows, size_t colStart, size_t numCols )
	in {
		assert( rowStart + numRows <= other.rows_ && colStart + numCols <= other.cols_,
			    sliceMsg_( rowStart, numRows, colStart, numCols ) );
	} body {
		data_    = other.data_;
		rows_    = numRows;
		cols_    = numCols;
		leading_ = other.leading_;
		ptr_     = other.ptr_ + mapIndex(rowStart, colStart);
	}
	
	/// ditto
	this( typeof(this)* other, size_t firstIndex, size_t numRows, size_t numCols )
	in {
		
	} body {
		data_    = other.data_;
		rows_    = numRows;
		cols_    = numCols;
		leading_ = other.leading_;
		ptr_     = other.ptr_ + firstIndex;
	}
	
	/** Resize the matrix and set all the elements to zero. */
	void resize( size_t newRows, size_t newCols ) {
		resize( newRows, newCols, null );
		generalMatrixScaling!storageOrder( rows_, cols_, Zero!ElementType, ptr_, leading_ );
	}
	
	/** Resize the matrix and leave the elements uninitialized. */
	void resize( size_t newRows, size_t newCols, void* ) {
		auto newLength = newRows * newCols;
		if( newLength != data_.length || data_.refCount() > 1 ) {
			data_.reset( newLength );
			ptr_     = data_.ptr;
			rows_    = newRows;
			cols_    = newCols;
			leading_ = minor_;
		}
	}
	
	/** Assignment has copy semantics. The actual copy is only performed on modification of the copy however. */
	ref typeof( this ) opAssign( CowMatrix rhs ) {
		data_    = move( rhs.data_ );
		ptr_     = rhs.ptr_;
		rows_    = rhs.rows_;
		cols_    = rhs.cols_;
		leading_ = rhs.leading_;
		
		return this;
	}
	
	/** Element access. */
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return ptr_[ mapIndex( i, j ) ];
	}
	
	/// ditto
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		unshareData_();
		mixin( "ptr_[ mapIndex( i, j ) ]" ~ op ~ "= rhs;" );
	}
	
	/** Remove the first major subvector (e.g. column for column major matrices). Part of the BidirectionalRange
	    concept.
	*/
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		-- major_;
		ptr_ += leading_;
	}
	
	/** Remove the last major subvector (e.g. column for column major matrices). Part of the BidirectionalRange
	    concept.
	*/
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack on empty." );
	} body {
		-- major_;
	}
	
	@property {
		/** Get a const pointer to the memory used by this storage. */
		const(ElementType*) cdata() const {
			return ptr_;
		}
		
		/** Get a mutable pointer to the memory used by this storage. */
		ElementType* data() {
			unshareData_();
			return ptr_;
		}
		
		/** Returh the length of the range (number of major subvectors). Part of the BidirectionalRange concept. */
		size_t length() const {
			return major_;
		}
		
		/** Is the array empty? Part of the BidirectionalRange concept. */
		bool empty() const {
			return major_ == 0;
		}
		
		/** Get the leading dimesnion of the matrix. */
		size_t leading() const {
			return leading_;
		}
		
		/** Get the number of rows. */
		size_t rows() const {
			return rows_;
		}
		
		/** Get the number of columns. */
		size_t columns() const {
			return cols_;
		}
		
		/** Get the number of major subvectors. */
		size_t major() const { 
			return major_;
		}
		
		/** Get the number of minor subvectors. */
		size_t minor() const {
			return minor_;
		}
		
		/** Get the address of this. Needed for a hack to avoid copying in certain cases. */
		typeof(this)* ptr() {
			return &this;
		}
	}
	
	size_t mapIndex( size_t i, size_t j ) const {
		if( isRowMajor )
			return i * leading_ + j;
		else
			return j * leading_ + i;
	}
	
private:
	mixin MatrixErrorMessages;

	static if( isRowMajor ) {
		alias rows_ major_;
		alias cols_ minor_;
	} else {
		alias rows_ minor_;
		alias cols_ major_;
	}
	
	void unshareData_() {
		if( data_.refCount() == 1 )
			return;
		
		if( leading_ == minor_ ) {
			auto len = rows_ * cols_;
			if( ptr_ == data_.ptr && len == data_.length )
				data_.unshare();
			else
				data_.reset( ptr_[ 0 .. len ] );
		} else {
			auto oldp = ptr_; // save the old ptr
			
			// NOTE: oldp won't be invalidated, because we know data_ is shared.
			data_.reset( rows_ * cols_ );
			auto newp = data_.ptr;
			
			generalMatrixCopy!( storageOrder, storageOrder )(
				rows_, cols_,
				oldp, leading_,
				newp, minor_
			);
			
			leading_ = minor_;
		}
		
		ptr_ = data_.ptr;
	}
	
	size_t        leading_, rows_, cols_;
	Data          data_;
	ElementType*  ptr_;
}

/** A simple alias for the preferred reference type for this container. */
template CowMatrixRef( T, StorageOrder order_ = StorageOrder.ColumnMajor ) {
	alias RefCounted!(CowMatrix!(T,order_), RefCountedAutoInitialize.yes )
			CowMatrixRef;
}

unittest {
	// TODO: CowArray tests.
}