module scid.internal.cowarray;

import std.algorithm;
import std.range;
import std.conv;
import std.traits;
import std.stdio;

import scid.internal.arraydata;

version(unittest) {
    //import scid.common.testing;
	import std.stdio;
}

/** Copy-on-write array implementation. Wrapped by Vectors & Matrices for block storage. */
struct CowArray( T ) {
	alias T ElementType;
		
	/** Allocate a new array of a given size. */
	this( size_t length, T initializeWith = T.init ) {
		this( length, null );
		this[] = initializeWith;
	}
	
	/** Allocate a new array of a given size. Do not initialize.*/
	this( size_t length, void* ) {
		data_.reset( length );
		firstIndex_ = 0;
		length_     = length;
	}
	
	/** Allocate a new array and initialize it with a (semantic) copy of a slice of a given array. */
	this( typeof( this ) slicedArray, size_t i, size_t j ) {
		data_       = slicedArray.data_;
		firstIndex_ = i + slicedArray.firstIndex_;
		length_     = min(j - i, slicedArray.length - i);
	}
		
	/** Allocate a new array and initialize it with a copy of a given array. */
	this( T[] array ) {
		data_.reset( array );
		firstIndex_ = 0;
		length_     = data_.length;
	}
		
	/** Assign a different array to this one. */
	ref typeof(this) opAssign( typeof( this ) rhs ) {
		data_       = rhs.data_;
		firstIndex_ = rhs.firstIndex_;
		length_     = rhs.length_;
		return this;
	}
	
	/** Assign a built-in array to this one. */
	ref typeof(this) opAssign( T[] array ) {
		data_.reset( array );
		firstIndex_ = 0;
		length_     = array.length;
		return this;
	}
	
	/** Indexing operators. */
	T opIndex( size_t i ) const
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string( length ) );
	} body {
		return data_[ i + firstIndex_ ];
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t i )
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string( length ) );
	} body {
		unshareData_();
		data_[ i + firstIndex_ ] = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i )
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string( length ) );
	} body {
		unshareData_();
		mixin( "data_[ i + firstIndex_ ]" ~ op ~ "= rhs;" );
	}
	
	/** Slicing operators. */
	typeof(this) opSlice() {
		return opSlice( 0, length );
	}
	
	/// ditto
	typeof(this) opSlice( size_t i, size_t j ) {
		return typeof( this )( this, i, j );
	}
	
	/// ditto
	void opSliceAssign( S )( S rhs )  {
		static if( is( S : typeof(this) ) )
			this = rhs;
		else
			this.opSliceAssign( rhs, 0, length_ );
	}
	
	/// ditto
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) if( is( S : T ) ) {
		b = min(length_, b);
		unshareData_();
		size_t e = firstIndex_ + b;
		for( size_t i = firstIndex_ + a; i < e; ++ i )
			data_[ i ] = rhs;
		
	}
	
	/// ditto
	void opSliceAssign( Range )( Range rhs, size_t a , size_t b ) if( isInputRange!Range ) {
		b = min(length_, b);
		
		static if( __traits(compiles, rhs.length) )
			assert( b-a == rhs.length, "Length mismatch: " ~ to!string(b-a) ~ " != " ~ to!string(rhs.length) );
		
		unshareData_();
		size_t e = firstIndex_ + b;
		for( size_t i = firstIndex_ + a; i < e; ++ i ) {
			data_[ i ] = rhs.front;
			rhs.popFront();
		}
	}
	
	/** Returns the wrapped builtin array. Read-only. */
	@property const(T) [] cdata() const {
		assert(data_.isInitialized());
		return data_.array[ firstIndex_ .. firstIndex_ + length_ ];
	}
	
	/** Returns the wrapped builtin array. Writable. Do not save the array as the memory may be freed when you least expect it. */
	@property T[] data() {
		assert(data_.isInitialized());
		unshareData_();
		return data_.array[ firstIndex_ .. firstIndex_ + length_ ];
	}	
	
	/** Returns the length of the array. */
	@property size_t length() const { return length_; }
	
	/** BidirectionalRange method. Return the first element in the array. */
	@property T    front() const { return this[ 0 ]; }
	
	/** BidirectionalRange method. Return the last element in the array. */
	@property T    back() const  { return this[ length_ - 1]; }
	
	/** BidirectionalRange method. Return whether the array is empty. */
	@property bool empty() const { return length_ == 0; }
	
	/** BidirectionalRange method. Remove the first element. */
	void           popFront()    { -- length_; ++ firstIndex_; }
	
	/** BidirectionalRange method. Remove the last element. */
	void           popBack()     { -- length_; }
	
	// Duplicate the data if multiple CowArrays share it.
	private void unshareData_() {
		// if the data isn't shared do nothing
		if( data_.refCount() == 1 )
			return;
		
		// if the array encompasses all of the data
		// just duplicate the data
		if( firstIndex_ == 0 && length_ == data_.length ) {
			data_.unshare();
			return;
		}
		
		// otherwise duplicate only the slice that is needed
		data_.reset( data_.array[ firstIndex_ .. length_ + firstIndex_ ] );
		firstIndex_ = 0;
	}
	
	// TODO: replace firstIndex_ & length_ with a slice in data_.array
	private size_t      firstIndex_, length_;
	private ArrayData!T data_;
}

unittest {
	alias CowArray!int IntArr;
	
	// test constructors
	auto    a1 = IntArr( 3 );
	auto    a2 = IntArr( [1,2,3] );
	IntArr  a3;
	
	assert( a1.length == 3, "Length ctor test failed." );
	assert( a2.length == 3, "List ctor length test failed." );
	assert( a3.length == 0, "Default ctor test failed." );
	
	assert( a2[0] == 1 && a2[1] == 2 && a2[2] == 3, "List ctor value test failed." );
	
	
	// test copying
	a1 = a2;
	assert( a1[0] == 1 && a1[1] == 2 && a1[2] == 3, "Copy equality test failed." );
	
	a1[ 0 ] = 42;
	assert( a1[ 0 ] == 42 && a2[ 0 ] == 1, "Copy duplication test failed."  );
	
	// test range functionality
	int i = 1;
	foreach( x ; a2 )
		assert( x == (i++), "Foreach test failed." );
	
	foreach_reverse( x; a2 )
		assert( x == (--i), "Reverse foreach test failed." );
	
	// test string rep
	assert( to!string(a1) == "[42, 2, 3]", "toString test failed." );
	assert( to!string(a3)== "[]", "Empty toString test failed." );
	
	// test slicing
	IntArr primes = [2,3,5,7,11,13,17,19];
	a3 = primes[ 2 ..  6 ];
	assert( a3.length == 4, "Length is " ~ to!string( a3.length ) );
	assert( a3[ 0 ] == 5  &&
	        a3[ 1 ] == 7  &&
	        a3[ 2 ] == 11 &&
	        a3[ 3 ] == 13,
		   "Assignment to slice test failed." );
	
	a1[] = 42;
	foreach( x ; a1 )
		assert( x == 42, "Slice assignment to value test failed." );
	
	
	a3[ 2 .. a3.length ] = 42;
	assert( a3[ 0 ] == 5  &&
	        a3[ 1 ] == 7  &&
	        a3[ 2 ] == 42 &&
	        a3[ 3 ] == 42,
		   "Delimited slice assignment to value test failed." );
	
	primes[ 0 .. 4  ] = primes[ 4 .. primes.length ];
	for( i = 0; i < 4; ++ i )
		assert( primes[ i ] == primes [ i + 4 ],
			   "Slice to slice assignment test failed." );
}