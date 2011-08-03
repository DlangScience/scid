module scid.storage.external;

import scid.common.traits;
import scid.common.meta;
import scid.matrix;
import scid.ops.common;   
import std.string, std.array;
import std.exception;
import std.conv;
import scid.blas;

mixin template EmulateRefCounted() {
	static struct RefCounted {
		static bool isInitialized() pure { return true; }
	}
}

struct ExternalArray( T, DefaultType_ ) {
	alias T ElementType;
	alias DefaultType_ ArrayType;
	
	mixin EmulateRefCounted;
	
	this( Allocator )( size_t len, Allocator allocator, void* )
			if( isAllocator!Allocator ) {
		
		array_ = allocator.uninitializedArray!(T[])( len );
		enforce( array_, "ExternalArray: Allocator returned null." );
	}
	
	this( Allocator )( size_t len, Allocator allocator )
			if( isAllocator!Allocator ) {
		
		array_ = allocator.uninitializedArray!(T[])( len );
		enforce( array_, "ExternalArray: Allocator returned null." );
		blas.scal( array_.length, Zero!T, array_.ptr, 1 );
	}
	
	this( E, Allocator )( E[] initializer, Allocator allocator )
			if( isAllocator!Allocator && isConvertible!(E,T) ) {
		
		array_ = allocator.uninitializedArray!(T[])( initializer.length );
		
		enforce( array_, "ExternalArray: Allocator returned null." );
				
		foreach( i ; 0 .. array_.length )
			array_[ i ] = to!T( initializer[ i ] );
		
		
	}
	
	this()( T[] data )
	in {
		assert( !data.empty, "ExternalArray: Null/empty data." );
	} body {
		array_  = data;
	}
	
	this()( ExternalArray *other ) {
		array_ = other.array_;	
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		array_  = rhs.array_;
		return this;
	}
	
	ElementType index( size_t i ) const {
		return array_[ i ];
	}
	
	void indexAssign( string op="" )( ElementType rhs, size_t i ) {
		mixin("array_[ i ] " ~ op ~ "= rhs;");
	}
		
	@property {
		size_t        length() const { return array_.length; }
		const(T)*     cdata()  const { return array_.ptr; }
		T*            data()         { return array_.ptr; }
		typeof(this)* ptr()          { return &this; }
	}
	
private:
	T[]    array_;
}

struct ExternalMatrix( T, StorageOrder order_, DefaultType_ ) {
	alias T ElementType;
	alias DefaultType_ MatrixType;
	alias ArrayTypeOf!MatrixType ArrayType;
	alias order_ storageOrder;
	
	mixin EmulateRefCounted;
	
	enum isRowMajor = (order_ == StorageOrder.RowMajor);
	
	alias ExternalMatrix!(T, transposeStorageOrder!order_, TransposedOf!DefaultType_ )
		Transposed;
	
	this( Allocator )( size_t rows, size_t columns, Allocator allocator, void* )
			if( isAllocator!Allocator ) {
		array_ = allocator.uninitializedArray!(T[])( rows * columns );
		enforce( array_, "ExternalMatrix: Allocator returned null." );
		rows_ = rows;
		cols_ = columns;
	}
	
	this( Allocator )( size_t rows, size_t columns, Allocator allocator )
			if( isAllocator!Allocator ) {
		array_ = allocator.uninitializedArray!(T[])( rows * columns );
		blas.scal( array_.length, Zero!T, array_.ptr, 1 );
		enforce( array_, "ExternalMatrix: Allocator returned null." );
		rows_ = rows;
		cols_ = columns;
	}
	
	this( E, Allocator )( E[][] initializer, Allocator allocator ) if( isAllocator!Allocator && isConvertible!(E, ElementType) ) {
		if( initializer.length == 0 || initializer[0].length == 0)
			return;
		
		rows_ = initializer.length;
		cols_ = initializer[0].length;
		size_t len = rows_ * cols_;
		array_ = ( cast(T*) allocator.allocate( len * T.sizeof ) )[ 0 .. len ];
		enforce( array_, "ExternalMatrix: Allocator returned null." );
		
		foreach( i ; 0 .. rows )
			foreach( j ; 0.. columns )
				indexAssign( to!ElementType(initializer[i][j]), i, j );
	}
	
	this()( size_t major, T[] data )
	in {
		assert( major > 0, "ExternalMatrix: Zero dimension." );
		assert( data.length % major == 0 , format("ExternalMatrix: Invalid major dimension for given data (%d, %d).", major, data.length) );
	} body {
		array_ = data;
		major_ = major;
		minor_ = data.length / major;
	}
	
	this()( ExternalMatrix *other ) {
		array_ = other.array_;
		rows_  = other.rows_;
		cols_  = other.cols_;
	}
	
	ElementType index( size_t i, size_t j ) const {
		return array_[ mapIndex(i, j) ];
	}
	
	void indexAssign( string op="" )( ElementType rhs, size_t i, size_t j ) {
		mixin("array_[ mapIndex(i, j) ] " ~ op ~ "= rhs;");
	}
	
	@property {
		size_t        rows()    const { return rows_; }
		size_t        columns() const { return cols_; }
		size_t        minor()   const { return minor_; }
		size_t        major()   const { return major_; }
		const(T)*     cdata()   const { return array_.ptr; }
		T*            data()          { return array_.ptr; }
		typeof(this)* ptr()           { return &this; }
		
		alias minor leading;
	}
	
	size_t mapIndex( size_t i, size_t j ) const {
		static if( isRowMajor )
			return i * cols_ + j;
		else
			return j * rows_ + i;
	}
	
private:
	T[] array_;
	size_t rows_, cols_;
	static if( isRowMajor ) {
		alias rows_ major_;
		alias cols_ minor_;
	} else {
		alias cols_ major_;
		alias rows_ minor_;
	}
}