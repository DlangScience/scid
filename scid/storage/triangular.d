module scid.storage.triangular;

import scid.internal.assertmessages;
import scid.storage.cowarray;
import scid.storage.packedmat;
import scid.matrix, scid.vector;
import scid.common.storagetraits;
import scid.ops.expression;
import std.math, std.algorithm;
import std.exception;
import scid.storage.external;

template TriangularStorage( ElementOrArray, MatrixTriangle triangle = MatrixTriangle.Upper, StorageOrder storageOrder = StorageOrder.ColumnMajor )
	if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias PackedStorage!( TriangularArrayAdapter!(CowArrayRef!ElementOrArray, triangle, storageOrder) ) TriangularStorage;
	else
		alias PackedStorage!( TriangularArrayAdapter!(ElementOrArray, triangle, storageOrder) )             TriangularStorage;
}

struct TriangularArrayAdapter( ContainerRef_, MatrixTriangle tri_, StorageOrder storageOrder_ ) {
	alias ContainerRef_                ContainerRef;
	alias BaseElementType!ContainerRef ElementType;
	alias ContainerRef                 ArrayType;
	
	alias TriangularArrayAdapter!(
		ContainerRef,
		tri_ == MatrixTriangle.Upper ? MatrixTriangle.Lower : MatrixTriangle.Upper,
		storageOrder_ 
	) Transposed;
	
	alias TriangularArrayAdapter!( ExternalArray!(ElementType, ArrayTypeOf!ContainerRef), tri_, storageOrder_ )
		Temporary;
	
	enum triangle     = tri_;
	enum storageOrder = storageOrder_;
	enum storageType  = MatrixStorageType.Triangular;
	enum isRowMajor   = storageOrder == StorageOrder.RowMajor;
	enum isUpper      = triangle     == MatrixTriangle.Upper;
	
	
	this( A ... )( size_t newSize, A arrayArgs ) {
		size_  = newSize;
		containerRef_ = ContainerRef( newSize * (newSize + 1) / 2, arrayArgs );
	}
	
	this()( ElementType[] initializer ) {
		auto tri  = (sqrt( initializer.length * 8.0 + 1.0 ) - 1.0 ) / 2.0;
		
		enforce( tri - cast(int) tri <= 0, msgPrefix_ ~ "Initializer list is not triangular." );
		
		size_  = cast(size_t) tri;
		containerRef_ = ContainerRef( initializer );
	}
	
	this()( ElementType[][] initializer ) {
		if( !initializer.length )
			return;
		
		size_  = initializer.length;
		containerRef_ = ContainerRef( packedArrayLength(size_) , null );
		
		foreach( i ; 0 .. size_ ) {
			static if( isUpper ) {
				foreach( j ; i .. size_ )
					this.indexAssign( initializer[ i ][ j ], i, j );
			} else {
				foreach( j ; 0 .. (i+1) )
					this.indexAssign( initializer[ i ][ j ], i, j );
			}
		}
	}
	
	this()( TriangularArrayAdapter *other ) {
		containerRef_ = ContainerRef( other.containerRef_.ptr );
		size_  = other.size_;
	}
	
	void resize( A ... )( size_t newRows, size_t newCols, A arrayArgs ) {
		size_t arrlen = packedArrayLength(newRows);
		
		enforce( newRows == newCols );
		
		if( !isInitd_ )
			containerRef_ = ContainerRef( arrlen, arrayArgs );
		else
			containerRef_.resize( arrlen, arrayArgs );
		
		size_ = newRows;
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		move( rhs.containerRef_, containerRef_ );
		size_  = rhs.size_;
		return this;
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
		
		return containerRef_.index( map_( i, j ) );
	}
	
	private void indexAssign(string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < size_ && j < size_, boundsMsg_(i,j) );
		static if( isUpper )
			assert( i <= j, "Modification of zero element in triangle matrix." );
		else
			assert( i >= j, "Modification of zero element in triangle matrix." );
	} body {
		containerRef_.indexAssign!op(rhs, map_( i, j ) );
	}
	
	@property {
		typeof(this)*       ptr()         { return &this; }
		ElementType*        data()        { return containerRef_.data; }
		const(ElementType)* cdata() const { return containerRef_.cdata; }
		size_t              size()  const { return size_; }
	}
	
	alias size rows;
	alias size columns;
	alias size major;
	alias size minor;
	
	template Promote( Other ) {
		private import scid.storage.generalmat;
	}
	
private:
	bool isInitd_() const {
		// TODO: This assumes the reference type is RefCounted. Provide a more general impl.
		return containerRef_.RefCounted.isInitialized();	
	}
	
	mixin MatrixErrorMessages;

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
	ContainerRef containerRef_;
}

