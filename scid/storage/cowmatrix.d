module scid.storage.cowmatrix;

import scid.internal.assertmessages;
import scid.internal.hlblas;
import scid.bindings.blas.dblas;
import scid.storage.arraydata;
import scid.matrix, scid.common.meta;

import std.algorithm, std.array, std.typecons;


struct CowMatrix( ElementType_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor ) {
	alias ElementType_                                 ElementType;
	alias ArrayData!ElementType                        Data;
	alias storageOrder_                                storageOrder;	
	
	enum isRowMajor = (storageOrder == StorageOrder.RowMajor);
	
	this( size_t newRows, size_t newCols )
	in {
		assert( newRows != 0 && newCols != 0, zeroDimMsg_( newRows, newCols ) );
	} body {
		this( newRows, newCols, null );
		scal( newRows * newCols, Zero!ElementType, ptr_, 1 );
	}
	
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
	
	this( typeof(this)* other ) {
		data_    = other.data_;
		ptr_     = other.ptr_;
		rows_    = other.rows_;
		cols_    = other.cols_;
		leading_ = other.leading_;
	}
	
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
	
	this( typeof(this)* other, size_t firstIndex, size_t numRows, size_t numCols )
	in {
		
	} body {
		data_    = other.data_;
		rows_    = numRows;
		cols_    = numCols;
		leading_ = other.leading_;
		ptr_     = other.ptr_ + firstIndex;
	}
	
	void resizeOrClear( size_t newRows, size_t newCols ) {
		resizeOrClear( newRows, newCols, null );
		hlGeneralMatrixScal!storageOrder( rows_, cols_, Zero!ElementType, ptr_, leading_ );
	}
	
	void resizeOrClear( size_t newRows, size_t newCols, void* ) {
		auto newLength = newRows * newCols;
		if( newLength != data_.length ) {
			data_.reset( newLength );
			ptr_     = data_.ptr;
			rows_    = newRows;
			cols_    = newCols;
			leading_ = minor_;
		}
	}
	
	ref typeof( this ) opAssign( CowMatrix rhs ) {
		data_    = move( rhs.data_ );
		ptr_     = rhs.ptr_;
		rows_    = rhs.rows_;
		cols_    = rhs.cols_;
		leading_ = rhs.leading_;
		
		return this;
	}
	
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return ptr_[ mapIndex( i, j ) ];
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		unshareData_();
		mixin( "ptr_[ mapIndex( i, j ) ]" ~ op ~ "= rhs;" );
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		-- major_;
		ptr_ += leading_;
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack on empty." );
	} body {
		-- major_;
	}
	
	@property {
		ElementType*        data()             { unshareData_(); return ptr_; }
		const(ElementType*) cdata()      const { return ptr_; }
		size_t              leading()    const { return leading_; }
		size_t              rows()       const { return rows_; }
		size_t              columns()    const { return cols_; }
		size_t              major()      const { return major_; }
		size_t              minor()      const { return minor_; }
		bool                empty()      const { return major_ == 0; }
		typeof(this)*       ptr()              { return &this; }
	}
	
	size_t mapIndex( size_t i, size_t j ) const {
		if( isRowMajor )
			return i * leading_ + j;
		else
			return j * leading_ + i;
	}
	
private:
	mixin SliceIndex2dMessages;

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
			
			hlGeneralMatrixCopy!( storageOrder, storageOrder )(
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

template CowMatrixRef( T, StorageOrder order_ = StorageOrder.ColumnMajor ) {
	alias RefCounted!(CowMatrix!(T,order_), RefCountedAutoInitialize.no )
			CowMatrixRef;
}

unittest {
	// TODO: CowArray tests.
}