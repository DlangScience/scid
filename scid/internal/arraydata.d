module scid.internal.arraydata;

import std.algorithm;
import std.conv;
import std.stdio;
import core.stdc.stdlib;

/**
Defines a reference-counted array of a primitive type $(D T). It is used by
$(D Vector) and $(D Matrix) storage types to implement copy-on-write. It is
unsafe and shouldn't really be used for anything else.

Changing the length of the array can only be done by creating a new array.
Doesn't support slicing.
*/

struct ArrayData( T ) {
	/** Allocate a new array of a given length. */
	this( size_t length ) {
		realloc_( length );
	}
	
	/** Allocate a new array whose elements are a copy of a given builtin array. */
	this( T[] array ) {
		realloc_( array.length );
		data_.array[] = array[];
	}
	
	/** Postblit constructor that increments the reference count accordingly. */
	this( this ) {
		if( !data_ )
			return;
		
		++ data_.refCount;
		debug( ArrayData )
			writeln( "Incremented refCount of ", data_, " to ", data_.refCount );
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
	@property T[] array() {
		return data_ is null ? null : data_.array;
	}
	
	/** Clears the reference - makes it reference an empty array. Does not affect other
	  * $(D ArrayData) objects pointing to the array.
	  */
	void clear() {
		deref_();
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
		
		ArrayData_ *oldData = data_;
		realloc_( data_.length );
		
		data_.array[] = oldData.array[];
	}
	
	/** Assignment updates the ref count accordingly. */
	void opAssign( typeof( this ) rhs ) {
		swap( data_, rhs.data_ );
	}
	
	/** Indexing operators. */
	T opIndex( size_t i ) const
	in {
		assert( data_ );
		assert( i < data_.length );
	} body {
		return *(data_.array_.ptr+i);
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t i )
	in {
		assert( data_ );
		assert( i < data_.length, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(data_.length) );
	} body {
		*(data_.array_.ptr+i) = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i )
	in {
		assert( data_ );
		assert( i < data_.length );
	} body {
		auto ref e = *(data_.array_.ptr+i);
		e = mixin( "e" ~ op ~ "rhs" );
	}
	
	/** Equality of pointers, not of data. Not using == to avoid confusion */
	bool isSharingDataWith( ref typeof( this ) whom ) {
		return (data_ == whom.data_) && (data_ !is null);
	}
	
	/** String representation of the underlying array. */
	string toString() const {
		if( !data_ )
			return "[]";
		else
			return to!string( data_.array );
	}
	
	// This the data that ArrayData points to. Uses struct hack to
	// allocate everything in place.
	private struct ArrayData_ {
		size_t refCount;
		size_t length;
		
		private T array_[1];
		
		@property T[] array() {
			return array_.ptr[ 0 .. this.length ];
		}
		
		@property const(T)[] array() const {
			return array_.ptr[ 0 .. this.length ];
		}
	}
	
	// Derferences the data, decrementing the ref count and freeing
	// the data if neccessary. It also sets data_ to null.
	private void deref_() {
		if( !data_ ) return;
		assert( data_.refCount > 0, "Zero refCount." );
		if( --data_.refCount == 0 ) {
			debug(ArrayData)
				writeln("Array ", data_, " freed.");
			free( data_ );
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
		
		const sz = sizeOfData_( length );
		auto p = malloc( sz )[ 0 .. sz ];
		// ASK: addRange()
		emplace!ArrayData_(  p[0 .. ArrayData_.sizeof] );
		data_ = cast( typeof( data_ ) ) p.ptr;
		data_.refCount = 1;
		data_.length   = length;
		
		debug(ArrayData)
			writeln( "Allocated array ", data_, " of ", sz, " bytes." );
	}
	
	// Returns the number of bytes that need to be allocated for an
	// array of the given length.
	static private size_t sizeOfData_( size_t length ) {
		return (*data_).sizeof + T.sizeof*(length - 1);
	}
	
	// A pointer to the referenced data.
	private ArrayData_ *data_;
}
