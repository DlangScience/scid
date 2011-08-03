module scid.storage.symmetric;

import scid.internal.assertmessages;
import scid.storage.cowarray;
import scid.storage.packedmat;
import scid.matrix, scid.vector;
import scid.common.storagetraits;
import scid.ops.common;
import std.math, std.algorithm;
import std.complex, std.exception;
import scid.storage.external;


template SymmetricStorage( ElementOrArray, MatrixTriangle triangle = MatrixTriangle.Upper, StorageOrder storageOrder = StorageOrder.ColumnMajor )
	if( isScalar!(BaseElementType!ElementOrArray) ) {
	
	static if( isScalar!ElementOrArray )
		alias PackedStorage!( SymmetricArrayAdapter!(CowArrayRef!ElementOrArray, triangle, storageOrder) ) SymmetricStorage;
	else
		alias PackedStorage!( SymmetricArrayAdapter!(ElementOrArray, triangle, storageOrder) )             SymmetricStorage;
}

struct SymmetricArrayAdapter( ContainerRef_, MatrixTriangle tri_, StorageOrder storageOrder_ ) {
	alias ContainerRef_                ContainerRef;
	alias BaseElementType!ContainerRef ElementType;
	alias ContainerRef                 ArrayType;
	
	alias SymmetricArrayAdapter!(
		ContainerRef,
		tri_ == MatrixTriangle.Upper ? MatrixTriangle.Lower : MatrixTriangle.Upper,
		storageOrder_ 
	) Transposed;
	
	alias SymmetricArrayAdapter!( ExternalArray!(ElementType, ArrayTypeOf!ContainerRef), tri_, storageOrder_ )
		Temporary;
	
	enum triangle     = tri_;
	enum storageOrder = storageOrder_;
	enum isRowMajor   = storageOrder == StorageOrder.RowMajor;
	enum isUpper      = triangle     == MatrixTriangle.Upper;
	
	/** Is the matrix hermitian? */
	enum isHermitian = isComplexScalar!ElementType;
	
	static if( isHermitian )
		enum storageType  = MatrixStorageType.Hermitian;
	else
		enum storageType  = MatrixStorageType.Symmetric;
	
	this( A ... )( size_t newSize, A arrayArgs ) {
		size_  = newSize;
		static if( A.length == 0 || !is( A[ 0 ] : size_t ) ) {
			containerRef_ = ContainerRef( newSize * (newSize + 1) / 2, arrayArgs );
		} else {
			assert( arrayArgs[0] == newSize,
				format( "Non-square dimensions for triangular matrix (%d,%d).", newSize, arrayArgs[ 0 ] ) );
			containerRef_ = ContainerRef( newSize * (newSize + 1) / 2, arrayArgs[ 1 .. $ ] );
		}
	}
	
	this()( ElementType[] initializer ) {
		auto tri  = (sqrt( initializer.length * 8.0 + 1.0 ) - 1.0 ) / 2.0;
		
		assert( tri - cast(int) tri <= 0, msgPrefix_ ~ "Initializer list is not triangular." );
		
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
	
	this()( SymmetricArrayAdapter *other ) {
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
		swap( rhs.containerRef_, containerRef_ );
		size_  = rhs.size_;
		return this;
	}

	ElementType index( size_t row, size_t column ) const
	in {
		assert( row < size_ );
		assert( column < size_ );
	} body {
		if( needSwap_( row, column ) ) {
			static if( isHermitian ) {
				return gconj( containerRef_.index( map_( column, row ) ) );
			} else {
				return containerRef_.index( map_( column, row ) );
			}
		} else {
			return containerRef_.index( map_( row, column ) );
		}
	}	

	void indexAssign(string op="")( ElementType rhs, size_t row, size_t column )
	in {
		assert( row < size_ );
		assert( column < size_ );
	} body {
		if( needSwap_( row, column ) ) {
			static if( isHermitian ) {
				containerRef_.indexAssign!op( gconj(rhs), map_( column, row ) );
			} else {
				containerRef_.indexAssign!op( rhs, map_( column, row ) );
			}
		} else {
			containerRef_.indexAssign!op( rhs, map_( row, column ) );
		}
	}
	
	@property {
		typeof(this)*       ptr()         { return &this; }
		ElementType*        data()        { return isInitd_() ? containerRef_.data  : null; }
		const(ElementType)* cdata() const { return isInitd_() ? containerRef_.cdata : null; }
		size_t              size()  const { return size_; }
	}
	
	alias size rows;
	alias size columns;
	alias size major;
	alias size minor;
	
	template Promote( Other ) {
		static if( isScalar!Other ) {
			alias SymmetricArrayAdapter!( Promotion!(ContainerRef,Other), triangle, storageOrder )
				Promote;
		}
	}
	
private:
	mixin MatrixErrorMessages;

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

	bool needSwap_( size_t i, size_t j ) const {
		static if( isUpper ) {
			return i > j;
		} else {
			return i < j;
		}
	}
	
	bool isInitd_() const {
		// TODO: This assumes the reference type is RefCounted. Provide a more general impl.
		return containerRef_.RefCounted.isInitialized();	
	}
	
	size_t   size_;
	ContainerRef containerRef_;
}