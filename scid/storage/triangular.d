module scid.storage.triangular;

import scid.internal.assertmessages;
import scid.storage.cowarray;
import scid.storage.packedmat;
import scid.matrix, scid.vector;
import scid.common.traits;
import std.math, std.algorithm;

template TriangularStorage( ElementOrArray, MatrixTriangle triangle = MatrixTriangle.Upper, StorageOrder storageOrder = StorageOrder.ColumnMajor )
	if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias PackedStorage!( TriangularArrayAdapter!(CowArrayRef!ElementOrArray, triangle, storageOrder) ) TriangularStorage;
	else
		alias PackedStorage!( TriangularArrayAdapter!(ElementOrArray, triangle, storageOrder) )             TriangularStorage;
}

struct TriangularArrayAdapter( ArrayRef_, MatrixTriangle tri_, StorageOrder storageOrder_ ) {
	alias ArrayRef_                ArrayRef;
	alias BaseElementType!ArrayRef ElementType;
	
	enum triangle     = tri_;
	enum storageOrder = storageOrder_;
	enum storageType  = MatrixStorageType.Triangular;
	enum isRowMajor   = storageOrder == StorageOrder.RowMajor;
	enum isUpper      = triangle     == MatrixTriangle.Upper;
	
	this( size_t newSize, ElementType initWith = ElementType.init ) {
		size_  = newSize;
		array_ = ArrayRef( newSize * (newSize + 1) / 2, initWith );
	}
	
	this( size_t newSize, void* ) {
		size_  = newSize;
		array_ = ArrayRef( newSize * (newSize + 1) / 2, null );
	}
	
	this( ElementType[] initializer ) {
		auto tri  = (sqrt( initializer.length * 8.0 + 1.0 ) - 1.0 ) / 2.0;
		
		assert( tri - cast(int) tri <= 0, msgPrefix_ ~ "Initializer list is not triangular." );
		
		size_  = cast(size_t) tri;
		array_ = ArrayRef( initializer );
	}
	
	this( typeof(this) *other ) {
		array_ = ArrayRef( other.array_.ptr );
		size_  = other.size_;
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		move( rhs.array_, array_ );
		size_  = rhs.size_;
		return this;
	}
	
	void resize( size_t size ) {
		array_ = ArrayRef( size*(size+1)/2, null );
		size_  = size;
	}
	
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < size_ && j < size_, boundsMsg_(i,j) );
	} body {
		static if( isUpper ) {
			if( i > j )
				return zero_;
		} else {
			if( i < j )
				return zero_;
		}
		
		return array_.index( map_( i, j ) );
	}
	
	private void indexAssign(string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < size_ && j < size_, boundsMsg_(i,j) );
		static if( isUpper )
			assert( i <= j, "Modification of zero element in triangle matrix." );
		else
			assert( i >= j, "Modification of zero element in triangle matrix." );
	} body {
		array_.indexAssign!op(rhs, map_( i, j ) );
	}
	
	@property {
		typeof(this)*       ptr()         { return &this; }
		ElementType*        data()        { return array_.data; }
		const(ElementType)* cdata() const { return array_.cdata; }
		size_t              size()  const { return size_; }
	}
	
	alias size rows;
	alias size columns;
	alias size major;
	alias size minor;
	
private:
	mixin SliceIndex2dMessages;

	enum zero_ = cast(ElementType)(0.0);
		
	size_t mapHelper_( bool colUpper )( size_t i, size_t j ) const {
		static if( colUpper ) return i + j * (j + 1) / 2;
		else                  return i + ( ( size_ + size_ - j - 1 ) * j ) / 2;
	}
	
	size_t map_( size_t i, size_t j ) const {
		static if( isRowMajor )
			return mapHelper_!( !isUpper )( j, i );
		else
			return mapHelper_!( isUpper )( i, j );
	}
	
	size_t   size_;
	ArrayRef array_;
}

