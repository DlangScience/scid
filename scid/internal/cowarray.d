module scid.internal.cowarray;

import std.algorithm;
import std.range;
import std.conv;
import std.traits;

import scid.internal.arraydata;

version(unittest) {
    import scid.common.testing;
	import std.stdio;
}

/** Copy-on-write array implementation. Wrapped by $(D ArrayStorage) to provide a storage type for $(D Vector). */
struct CowArray( T ) {
	alias T ElementType;
		
	/** Allocate a new array of a given size. */
	this( size_t length ) {
		data_       = typeof( data_ )( length );
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
		data_       = typeof( data_ )( array );
		firstIndex_ = 0;
		length_     = data_.length;
	}
		
	/** Assign a different array to this one. */
	void opAssign( typeof(this) rhs ) {
		// Share the ArrayData. COW behaviour.
		data_       = rhs.data_;
		firstIndex_ = rhs.firstIndex_;
		length_     = rhs.length_;
	}
	
	/** Assign a built-in array to this one. */
	void opAssign( T[] array ) {
		data_       = typeof( data_ )( array );
		firstIndex_ = 0;
		length_     = array.length;
	}
	
	/** Indexing operators. */
	T opIndex( size_t i ) const
	in {
		assert( i < length_ );
	} body {
		return data_[ i + firstIndex_ ];
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t i )
	in {
		assert( i < length_ );
	} body {
		unshareData_();
		data_[ i + firstIndex_ ] = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i )
	in {
		assert( i < length_ );
	} body {
		unshareData_();
		mixin( "data_[ i + firstIndex_ ]" ~ op ~ "= rhs" );
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
	
	/// ditto
	size_t opDollar( int dim )() const { return length; }
	
	/** Returns the length of the array. */
	@property size_t length() const { return length_; }
	
	/** InputRange method. Return the first element in the array. */
	@property T    front() const { return this[ 0 ]; }
	
	/** InputRange method. Return the last element in the array. */
	@property T    back() const  { return this[ length_ - 1]; }
	
	/** InputRange method. Return whether the array is empty. */
	@property bool empty() const { return length_ == 0; }
	
	/** InputRange method. Remove the first element. */
	void           popFront()    { -- length_; ++ firstIndex_; }
	
	/** InputRange method. Remove the last element. */
	void           popBack()     { -- length_; }
	
	/** String representation of the Array */
	string toString() const {
		// TODO: Faster implementation.
		if( length == 0 )
			return "[]";
			
		string s = '[' ~ to!string( data_[ firstIndex_ ] );
		int i = firstIndex_ + 1;
		int e = firstIndex_ + length_;
		for( ; i < e; ++i)
			s ~= ", " ~ to!string( data_[ i ] );
		return s ~ ']';
	}
	
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
		data_       = typeof( data_ )( data_.array[ firstIndex_ .. length_ + firstIndex_ ] );
		firstIndex_ = 0;
	}
	
	private size_t      firstIndex_, length_;
	private ArrayData!T data_;
}

unittest {
	alias CowArray!int IntArr;
	
	// test constructors
	auto    a1 = IntArr( 3 );
	auto    a2 = IntArr( [1,2,3] );
	IntArr  a3;
	
	check( a1.length == 3, "Length ctor test failed." );
	check( a2.length == 3, "List ctor length test failed." );
	check( a3.length == 0, "Default ctor test failed." );
	
	check( a2[0] == 1 && a2[1] == 2 && a2[2] == 3, "List ctor value test failed." );
	
	
	// test copying
	a1 = a2;
	check( a1[0] == 1 && a1[1] == 2 && a1[2] == 3, "Copy equality test failed." );
	
	a1[ 0 ] = 42;
	check( a1[ 0 ] == 42 && a2[ 0 ] == 1, "Copy duplication test failed."  );
	
	// test range functionality
	int i = 1;
	foreach( x ; a2 )
		check( x == (i++), "Foreach test failed." );
	
	foreach_reverse( x; a2 )
		check( x == (--i), "Reverse foreach test failed." );
	
	// test string rep
	check( a1.toString() == "[42, 2, 3]", "toString test failed." );
	check( a3.toString() == "[]", "Empty toString test failed." );
	
	// test slicing
	IntArr primes = [2,3,5,7,11,13,17,19];
	a3 = primes[ 2 ..  6 ];
	check( a3.length == 4, "Length is " ~ to!string( a3.length ) );
	check( a3[ 0 ] == 5  &&
	        a3[ 1 ] == 7  &&
	        a3[ 2 ] == 11 &&
	        a3[ 3 ] == 13,
		   "Assignment to slice test failed." );
	
	a1[] = 42;
	foreach( x ; a1 )
		check( x == 42, "Slice assignment to value test failed." );
	
	
	a3[ 2 .. a3.length ] = 42;
	check( a3[ 0 ] == 5  &&
	        a3[ 1 ] == 7  &&
	        a3[ 2 ] == 42 &&
	        a3[ 3 ] == 42,
		   "Delimited slice assignment to value test failed." );
	
	primes[ 0 .. 4  ] = primes[ 4 .. primes.length ];
	for( i = 0; i < 4; ++ i )
		check( primes[ i ] == primes [ i + 4 ],
			   "Slice to slice assignment test failed." );
}