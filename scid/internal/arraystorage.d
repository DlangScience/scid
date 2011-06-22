module scid.internal.arraystorage;

import scid.internal.arrayviewstorage;
import scid.internal.cowarray;
import scid.vector;
//import scid.internal.arrayviewstorag;

import std.typecons : RefCounted;

/**
A reference-counted wrapper for CowArray that implements the Storage
concept: assign, at, set, slice, sliceAssign.
*/
struct ArrayStorage( T ) {
	alias ArrayViewStorage!T View;
	alias T                  ElementType;
	
	/** Allocates a new array of a given size. */
	this( size_t size ) {
		array_ = typeof( array_ )( size );
	}
	
	/** Allocates a new array and initializes it with a copy of a given array. */
	this( T[] array ) {
		array_ = typeof( array_ )( array );
	}
	
	/** Allocates a new array as a copy of a given $(D ArrayView) */
	this( ref View view ) {
		assign( view );
	}
	
	/// ditto
	this( ref BaseVector!View vec ) {
		assign( vec.storage );
	}
	
	/** Length of the array. */
	@property size_t length() const { return array_.length; }
	
	/** Assignment copies the reference. Use assign for copying elements. */
	ref typeof(this) opAssign( typeof(this) rhs ) {
		array_ = rhs.array_;
		return this;
	}
	
	/** Storage concept. Calls opAssign on underlying array. */
	void assign( typeof(this) rhs ) {
		array_.refCountedPayload() = rhs.array_.refCountedPayload();
	}
	
	/** Storage concept. Calls opAssign on underlying array. */
	void assign( View rhs ) {
		if( rhs.stride == 1 ) {        // If the stride is 1 they can use the same memory.
			assign( rhs.array );
			size_t a = rhs.firstIndex;
			size_t b = a + rhs.length;
			array_ = array_.refCountedPayload()[ a .. b ];
		} else {                       // otherwise we duplicate the elements.
			if( rhs.length != length )
				array_ = typeof( array_ )( rhs.length );
				
			array_.refCountedPayload()[] = rhs;
		}
	}
	
	/** Storage concept. Calls opIndexAssign on underlying array. */
	void set( T rhs, size_t i ) {
		array_.refCountedPayload()[ i ] = rhs;
	}
	
	/** Storage concept. Calls opIndex on underlying array. */
	T at( size_t i ) const {
		return array_.refCountedPayload()[ i ];
	}
	
	/** Storage concept. Creates a $(D View) with the given slice parameters. */
	View slice( size_t i, size_t j, size_t stride = 1 ) {
		return View( this, i, j, stride );	
	}
	
	/** Storage concept. Calls opSlice assign on the underlying array. */
	void sliceAssign( S )( S rhs, size_t a, size_t b ) {
		array_[ a .. b ] = rhs;
	}
	
	/** String representation. Forwards the call to the array_. */
	string toString() const {
		return array_.toString();
	}
	
	/** InputRange method. Returns the first element in the array. */
	@property T    front() const { return array_.front; }
	
	/** InputRange method. Returns the last element in the array. */
	@property T    back() const  { return array_.back; }
	
	/** InputRange method. Returns whether the array is empty. */
	@property bool empty() const { return array_.empty; }
	
	/** InputRange method. Removes the first element. */
	void           popFront()    { array_.popFront(); }
	
	/** InputRange method. Removes the last element. */
	void           popBack()     { array_.popBack(); }
	
	// The wrapped ref-counted CowArray.
	private RefCounted!( CowArray!T ) array_;
}