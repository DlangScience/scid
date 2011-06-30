module scid.vector;

import scid.vectorview;
import scid.internal.cowarray;
import std.stdio;

import std.typecons : RefCounted;
import std.conv     : to;

enum VectorType {
	Column, Row
}

/** Vector struct with value semantics. Provides the storage used as a template parameter.
  * 
  * Internally it wraps a RefCounted!(ArrayTemplate!T) which it shares with views associated
  * with the vector. This way, when the vector is freed the views are not invalidated. 
  *
  * Supports slicing and implements the BidirectionalRange concept. Slices are semantic copies
  * of the original, to get a slice that propagates edits to the original use VectorView's.
  * The view() method returns a VectorView pointing to a slice of this vector (a stride can
  * also be specified.)
  *
  * To use with BLAS functions (and other functions taking a builtin array), use the data()
  * and cdata() methods. Whenever possible, use cdata() as data() might duplicate the memory.
 */
struct Vector( T, VectorType vtype = VectorType.Column, alias ArrayTemplate = CowArray ) {
	/** The type of elements in the vector. */
	alias T                ElementType;
	
	/** The type of the Array that's wrapped. */
	alias ArrayTemplate!T    Wrapped;
	
	/** The RefCounted array type. */
	alias RefCounted!Wrapped Array;
	
	/** The type of views associated with this vector type. */
	alias VectorView!Wrapped View;
	
	/** Compile time flag to deisgnate the type of vector. */
	enum bool isRowVector = vtype == VectorType.Row;
	
	/** Allocates a new vector of a given size. */
	this( size_t len, T initializeWith = T.init ) {
		array_.RefCounted.initialize( len, initializeWith );
	}
	
	/** Allocates a new vector of a given size. Do not initialize. */
	this( size_t len, void* ) {
		array_.RefCounted.initialize( len, null );
	}
	
	/** Allocates a new vector of a given size. Do not initialize the data. */
	//this( size_t len, void ) {
		//array_.RefCounted.initialize( len, [] );
	//}
	
	
	/** Allocates a new vector with elements equal to those of a given array. */
	this( T[] arr ) {
		array_.RefCounted.initialize( arr );
	}
	
	/** Postblit needs calls opAssign on the wrapped array. */
	this( this ) {
		if( array_.RefCounted.isInitialized() )
			array_.RefCounted.initialize( *&array_.refCountedPayload() );
	}
	
	// Creates a Vector that's a slice of another vector. Called by op*Slice*().
	private this( ref Wrapped arr, size_t a, size_t b ) {
		array_.RefCounted.initialize( arr, a, b );
	}
	
	/** Length of the vector. */
	@property size_t length() const { return array_.RefCounted.isInitialized() ? array_.length : 0; }
	
	/** Return the ref-counted array. */
	@property Array array() { return array_; }
	
	/** Return the wrapped memory block. Writable - duplicates the memory if it's shared. */
	@property T[] data() {
		assert(array_.RefCounted.isInitialized());
		return array_.data;
	}
	
	/** Return the wrapped memory block. Read-only. */
	@property const(T) [] cdata() const {
		assert(array_.RefCounted.isInitialized());
		return array_.cdata;
	}
	
	/** Get a view of the array. */
	View view( size_t start, size_t end, size_t stride = 1 )
	in {
		assert( stride > 0 );
		assert( start <= end );
		assert( end <= length );
	} body {
		return View( array_, start, end, stride );
	}
	
	/// ditto
	View view( size_t stride = 1 )
	in {
		assert( stride > 0 );
	} body {
		return View( array_, 0, length, stride );
	}
	
	/** Calls opAssign on the underlying array. */
	ref typeof(this) opAssign( typeof(this) rhs ) {
		// Workaround for Bug 6199.
		array_.refCountedPayload() =  *&rhs.array_.refCountedPayload();
		return this;
	}
	
	/** Create a new vector whose elements are copies of those of an array. */
	ref typeof(this) opAssign( T[] rhs ) {
		array_.refCountedPayload() = rhs;
		return this;
	}
	
	/** Assignment to a view. Copies the elements contained in the view.
	  * The new size will be that of the view. If the stride of the view
	  * is one, no actual copying is performed.
	  */
	ref typeof(this) opAssign( View rhs ) {
		if( rhs.stride == 1 ) {        // If the stride is 1 they can use the same memory.
			array_.refCountedPayload() =
				rhs.array.refCountedPayload()[
					rhs.firstIndex .. rhs.firstIndex + rhs.length
				]
			;
		} else {                       // otherwise we duplicate the elements.
			if( rhs.length != length )
				array_ = typeof( array_ )( rhs.length );
			
			// use range copying
			array_.refCountedPayload()[] = rhs;
		}
		
		return this;
	}
	
	/** Access to the ith element of the vector. */
	T opIndex( size_t i ) const {
		return array_.refCountedPayload()[ i ];
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t i ) {
		array_.refCountedPayload()[ i ]	= rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i ) {
		mixin( "array_.refCountedPayload()[ i ] " ~ op ~ "= rhs;" );
	}
	
	/** Returns a copy of this vector. */
	typeof( this ) opSlice() {
		// return this; - this doesn't work because of a bug whose number I'm forgetting right now, sorry.
		return opSlice( 0, length );
	}
	
	/** Returns a slice of this vector. No copying is actually performed. */
	typeof(this) opSlice( size_t i, size_t j ) {
		return typeof(this)( array_, i, j );
	}
	
	/** Slice assignment to a range or a value. Forwarded to the wrapped array. */
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) {
		array_.refCountedPayload()[ a .. b ] = rhs;
	}
	
	/// ditto
	void opSliceAssign( S )( S rhs ) {
		array_.refCountedPayload()[] = rhs;
	}
	
	/** BidrectionalRange method. Returns the first element in the array. */
	@property T    front() const { return array_.front; }
	
	/** BidrectionalRange method. Returns the last element in the array. */
	@property T    back() const  { return array_.back; }
	
	/** BidrectionalRange method. Returns whether the array is empty. */
	@property bool empty() const { return array_.empty; }
	
	/** BidrectionalRange method. Removes the first element. */
	void           popFront()    { array_.popFront(); }
	
	/** BidrectionalRange method. Removes the last element. */
	void           popBack()     { array_.popBack(); }
	
	// The wrapped ref-counted CowArray.
	private RefCounted!( CowArray!T ) array_;
}

version( unittest ) {
	import std.algorithm;
	import std.stdio;
}

unittest {
	alias Vector!int Vec;
	
	// static assert( isVector!Vec );
	
	// test ctors
	auto a = Vec(2);
	auto b = Vec( [1,2,3] );
	Vec  c;

	assert( a.length == 2, "Length ctor failed" );
	assert( b.length == 3 && b[0] == 1 && b[1] == 2 && b[2] == 3, "Array ctor failed" );
	assert( c.length == 0, "Default ctor failed" );

	// vector assignment
	a = b;
	assert( a.length == b.length );
	a[ 0 ] += 41;
	
	assert( a[ 0 ] == 42 && b[ 0 ] == 1, "Assignment failed" );

	// array assignment
	a = [10,9,8,7,6,5,4,3,2,1];
	for(int i = 10; i>0; --i)
		assert( a[ 10 - i ] == i, "Array assignment failed" );

	// slicing
	a[0 .. 3] = 42;
	assert( a[0] == 42 && a[1] == 42 && a[2] == 42 && a[3] == 7, "Slice assignment to value failed" );
	
	a = [1,2,3];
	a[] = map!("a+1")(a);
	assert( a[0] == 2 && a[1] == 3 && a[2] == 4, "Slice assignment to range failed" );
	
	a = [1,2,3,4,5,6];
	a[0 .. 3] = a[ 3 .. 6 ];

	for( int i=0; i<3; ++i )
		assert( a[ i ] == a[ i + 3 ], "Slice to slice assignment failed" );


	auto d = a[0 .. 3];
	d[ 0 ] = 42;
	assert( d[ 0 ] == 42 && a [ 0 ] == 4, "Slice duplication failed" );
	
	// foreach
	// static assert( !__traits(compiles, { Vec a; foreach( ref x ; a ) {} }() ) );
	a = [1,2,3,4];
	{
		int i = 0;
		foreach( x; a )
			assert( x == ++i, "Foreach failed." );
	}

	// toString()
	assert( to!string(a) == "[1, 2, 3, 4]", "toString() failed" );

	a = [];
	assert( to!string(a) == "[]", "toString() on empty array failed" );
	
	// views
	
	a = [0,1,2,3,4,5,6,7,8,9];
	
	auto v03 = a.view(0,3);
	v03[ 0 ] = 42;
	assert( v03.length == 3 && a[ 0 ] == 42, "View assignment propagation failed" );
	
	b = v03;
	assert( b.length == 3 && b[0] == 42 && b[1] == 1 && b[2] == 2, "Vector=View assignment failed" );
	
	v03[ 0 ] = 0;
	assert( b[ 0 ] == 42 && a[ 0 ] == 0, "Vector-View duplication failed" );
	
	b = [10,20,30];
	v03[] = b;
	assert( a[ 0 ] == 10 && a[ 1 ] == 20 && a[ 2 ] == 30, "View=Vector assignment failed" );
	
	v03 = b.view();
	v03[ 0 ] = 1;
	assert( b[ 0 ] == 1, "View reassignment failed" );
	
	foreach( i ; 0..10 )
		a[ i ] = i;
	
	auto evens = a.view(0,10,2);
	auto odds  = a.view(1,10,2);
	
	assert( to!string(evens) == "[0, 2, 4, 6, 8]", "Strided view failed" );
	assert( to!string(odds)  == "[1, 3, 5, 7, 9]", "Offset strided view failed" );
	
	odds[] = 2;
	evens[] = 1;
	
	assert( to!string(a) == "[1, 2, 1, 2, 1, 2, 1, 2, 1, 2]", "Strided view assignment failed" );
	
}