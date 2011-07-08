module scid.storage.cowmatrix;

import scid.internal.assertmessages;
import scid.storage.arraydata;
import scid.matrix;

import std.algorithm, std.array, std.typecons;


struct CowMatrix( ElementType_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor ) {
	alias ElementType_                                 ElementType;
	alias ArrayData!ElementType                        Data;
	alias storageOrder_                                storageOrder;	
	
	enum isRowMajor = (storageOrder == StorageOrder.RowMajor);
	
	this( size_t newRows, size_t newCols, ElementType initWith = ElementType.init )
	in {
		assert( newRows != 0 && newCols != 0, zeroDimMsg_( newRows, newCols ) );
	} body {
		this( newRows, newCols, null );
		slice_[] = initWith;
	}
	
	this( size_t newRows, size_t newCols, void* )
	in {
		assert( newRows != 0 && newCols != 0, zeroDimMsg_( newRows, newCols ) );
	} body {
		data_.reset( newRows * newCols );
		slice_   = data_.array;
		rows_    = newRows;
		cols_    = newCols;
		leading_ = minor_;
	}
	
	this( size_t newMajor, ElementType[] initializer )
	in {
		assert( newMajor != 0 && initializer.length % newMajor == 0, initMsg_( newMajor, initializer ) );
	} body {
		data_.reset( initializer );
		slice_   = data_.array;
		major_   = newMajor;
		minor_   = initializer.length / major_;
		leading_ = minor_;
	}
	
	this( typeof(this)* other ) {
		data_    = other.data_;
		slice_   = other.slice_;
		rows_    = other.rows_;
		cols_    = other.cols_;
		leading_ = other.leading_;
	}
	
	this( typeof(this)* other, size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd )
	in {
		assert( rowStart < rowEnd && rowEnd <= other.rows_ && colStart < colEnd && colEnd <= other.cols_,
			    sliceMsg_( rowStart, rowEnd, colStart, colEnd ) );
	} body {
		data_    = other.data_;
		rows_    = rowEnd - rowStart;
		cols_    = colEnd - colStart;
		leading_ = other.leading_;
		slice_   = other.slice_[ mapIndex(rowStart, colStart) .. mapIndex(rowEnd-1, colEnd-1) + 1 ];
	}
	
	ref typeof( this ) opAssign( CowMatrix rhs ) {
		data_    = move( rhs.data_ );
		slice_   = rhs.slice_;
		rows_    = rhs.rows_;
		cols_    = rhs.cols_;
		leading_ = rhs.leading_;
		
		return this;
	}
	
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return slice_[ mapIndex( i, j ) ];
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		unshareData_();
		mixin( "slice_[ mapIndex( i, j ) ]" ~ op ~ "= rhs;" );
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		-- major_;
		if( leading_ >= slice_.length )
			slice_ = null;
		else
			slice_ = slice_[ leading_ .. $ ];	
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack on empty." );
	} body {
		-- major_;
		if( leading_ >= slice_.length )
			slice_ = null;
		else
			slice_ = slice_[ leading_ .. $ ];	
	}
	
	@property {
		ElementType[]        data()             { unshareData_(); return slice_; }
		const(ElementType[]) cdata()      const { return slice_; }
		size_t               length()     const { return slice_.length; }
		size_t               leading()    const { return leading_; }
		size_t               rows()       const { return rows_; }
		size_t               columns()    const { return cols_; }
		size_t               major()      const { return major_; }
		size_t               minor()      const { return minor_; }
		bool                 empty()      const { return slice_.empty; }
		typeof(this)*        ptr()              { return &this; }
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
		
		if( leading_ == major_ ) {
			if( slice_ == data_.array )
				data_.unshare();
			else
				data_.reset( slice_ );
		} else {
			auto oldData = slice_; // save the old array
			
			// NOTE: oldData won't be invalidated, because we know data_ is shared.
			data_.reset( rows_ * cols_ );
			auto newData = data_.array;
			
			while( true ) { // ends when oldData would be empty
				newData[ 0 .. minor_ ] = oldData[ 0 .. minor_ ];
				newData = newData[ minor_   .. $ ];
				
				if( leading_ >= oldData.length )
					break;
				else
					oldData = oldData[ leading_ .. $ ];
			}
			
			leading_ = minor_;
		}
		
		slice_   = data_.array;
	}
	
	size_t        leading_, rows_, cols_;
	Data          data_;
	ElementType[] slice_;
}

template CowMatrixRef( T, StorageOrder order_ = StorageOrder.ColumnMajor ) {
	alias RefCounted!(CowMatrix!(T,order_), RefCountedAutoInitialize.no )
			CowMatrixRef;
}

unittest {
	// TODO: CowArray tests.
}