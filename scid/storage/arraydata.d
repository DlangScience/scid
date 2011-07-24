/** Implementation of ArrayData, a reference counted array implementation which is wrapped by the containers used for
    storage.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.storage.arraydata;

import core.stdc.stdlib;
import std.algorithm, std.conv;

import scid.bindings.blas.dblas;
alias scid.bindings.blas.dblas blas;

// debug = ArrayData;

debug( ArrayData ) {
	debug = ArrayDataAllocs;
}

debug( ArrayDataAllocs ) {
	import std.stdio;
}

/**
  * Defines a reference-counted array of a primitive type T. It is used by
  * CowArray to implement copy-on-write. It is unsafe and shouldn't really
  * be used for anything else.
  * 
  * Changing the length of the array can only be done by creating a new array.
  * Doesn't support slicing or range functionality.
  */

struct ArrayData( T ) {
	/** Allocate a new array of a given length. */
	this( size_t length ) {
		reset( length );
	}
	
	/** Allocate a new array whose elements are a copy of a given builtin array. */
	this( T[] array ) {
		reset( array );
	}
	
	/** Postblit constructor that increments the reference count accordingly. */
	this( this ) {
		if( !data_ )
			return;
		
		++ data_.refCount;
		debug( ArrayData )
			writeln( "Incremented refCount of ", data_, " to ", data_.refCount );
	}
	
	/** Allocate a new array of a given length. */
	void reset( size_t length ) {
		realloc_( length );
	}
	
	/** Allocate a new array whose elements are a copy of a given builtin array. */
	void reset( T[] array ) {
		if( array.length == 0 ) {
			reset();
			return;
		}
		
		// If the given array overlaps with the current one and this is the only
		// reference to it, then we need to duplicate the given array.
		if( data_ && data_.refCount == 1 ) {
			if( array.ptr < data_.ptr + data_.length && data_.ptr < array.ptr + array.length ) {
				debug( ArrayDataAllocs )
					writeln("ArrayData: Duplicated array in reset() because of overlap.");
				array = array.dup;
			}
		}
				
		if( !isInitialized() || array.length != data_.length || data_.refCount > 1  )
			realloc_( array.length );
		
		blas.copy( data_.length, array.ptr, 1, data_.ptr, 1 );
	}
	
	/** Just clears the reference - makes this reference an empty array. Does not affect other
	  * ArrayData objects pointing to the array.
	  */
	void reset() {
		deref_();
	}
	
	/** Destructor that decrements the ref count and frees the data if needed. */
	~this() {
		deref_();
	}
	
	/** Returns the number of references to the data, or zero if there is no data. */
	size_t refCount() const {
		if( data_ )
			return data_.refCount;
		else
			return 0;
	}
	
	/** Returns the length of the array. */
	@property size_t length() const {
		return data_ is null ? 0 : data_.length;
	}
	
	/** Returns the array. */
	@property T* ptr() {
		return data_ is null ? null : data_.ptr;
	}
	
	/** Returns the array. */
	@property const(T)* ptr() const {
		return data_ is null ? null : data_.ptr;
	}
	
	/** Returns true if the pointer is inside the memory owned by the array. */
	bool owns( T* address ) const {
		if( data_ is null )
			return false;
		
		return address >= ptr && address < (ptr + length);
	}
	
	/** Returns true if there is any data referenced, false otherwise. */
	bool isInitialized() const {
		return data_ !is null;
	}
	
	/** If the data is shared, it makes a copy of it so that this object
	  * holds the only reference to it. Does nothing otherwise.
	  */
	void unshare() {
		if( !data_ || data_.refCount == 1 )
			return;
		
		Data_ *oldData = data_;
		realloc_( data_.length );
		
		blas.copy( length, oldData.array.ptr, 1, data_.array.ptr, 1 );
	}
	
	/** Assignment updates the ref count accordingly. */
	void opAssign( typeof( this ) rhs ) {
		swap( data_, rhs.data_ );
	}
	
	/** Indexing operators. */
	T opIndex( size_t i ) const
	in {
		assert( isInitialized() );
		assert( i < data_.length, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(data_.length) );
	} body {
		return *(data_.ptr+i);
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t i )
	in {
		assert( isInitialized() );
		assert( i < data_.length, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(data_.length) );
	} body {
		*(data_.ptr+i) = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i )
	in {
		assert( isInitialized() );
		assert( i < data_.length, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(data_.length) );
	} body {
		mixin( "*(data_.ptr+i)" ~ op ~ "= rhs;" );
	}
	
	/** Equality of pointers, not of data. Not using == to avoid confusion */
	bool isSharingDataWith( ref typeof( this ) whom ) {
		return (data_ == whom.data_) && (data_ !is null);
	}
	
	/** String representation of the underlying array. */
	string toString() const {
		if( !isInitialized() )
			return "[]";
		else
			return to!string( data_.array );
	}
	
	// This is the data that ArrayData points to. Uses struct hack to
	// allocate everything in place.
	private struct Data_ {
		size_t refCount      = void;
		size_t length        = void;
		byte   ptrOffset     = void;
		
		T array_[1]  = void;
		
		@property T[] array() {
			return array_.ptr[ 0 .. this.length ];
		}
		
		@property const(T)[] array() const {
			return array_.ptr[ 0 .. this.length ];
		}
		
		@property T*        ptr()       { return array_.ptr; }
		@property const(T)* ptr() const { return array_.ptr; }
	}
	
	// Derferences the data, decrementing the ref count and freeing
	// the data if neccessary. It also sets data_ to null.
	private void deref_() {
		if( !data_ ) return;
		assert( data_.refCount > 0, "Zero refCount." );
		--data_.refCount;
		if( data_.refCount == 0 ) {
			debug(ArrayDataAllocs)
				writeln("Array ", data_, " freed.");
			free( cast(void*)data_ - data_.ptrOffset );
		} else {
			debug(ArrayData)
				writeln( "Decremented refCount of array ", data_, " to ", data_.refCount );	
		}
		
		data_ = null;
	}
	
	// Reallocates the data. Takes care of deref'ing the current one
	// first.
	private void realloc_( size_t length ) {
		deref_();
		
		if( length == 0 )
			return;
		
		// alignment fix shouldn't be required, I think
		const sz = sizeOfData_( length ) + Data_.alignof;
		auto p = malloc( sz )[ 0 .. sz ];
		auto offset = cast(size_t)(p.ptr) % Data_.alignof;
		if( offset )
			offset =  Data_.alignof - offset;
		
		// ASK: addRange()
		data_ = emplace!Data_(  p[offset .. sz] );
		data_.ptrOffset = cast(byte) offset;
		data_.refCount  = 1;
		data_.length    = length;
		
		debug(ArrayDataAllocs)
			writeln( "Allocated array ", data_, " of ", length, " elements." );
	}
	
	// Returns the number of bytes that need to be allocated for an
	// array of the given length.
	static private size_t sizeOfData_( size_t length ) {
		return (*data_).sizeof + T.sizeof*(length - 1);
	}
	
	// A pointer to the referenced data.
	private Data_ *data_ = null;
}


unittest {
    alias ArrayData!double Test;

    auto x = Test(10);
    
    assert( x.refCount() == 1 );
    
    auto y = x;
    assert( y.refCount() == 2 );
    assert( x.isSharingDataWith( y ) );
    
    Test z;
    assert( !z.isInitialized() );
    
    z = y;
    assert( z.refCount() == 3 );
    
    z.reset();
    y.reset();
    assert( x.refCount() == 1 );
    x.reset();
}