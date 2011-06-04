module vector;

import expressions;
import vecclosure;


VectorView!( T ) vector( T )( size_t size ) pure {
	return typeof( return )( new T[ size ] );
}

VectorView!( T ) vector( T )( size_t size, T init ) pure {
	auto arr = new T[ size ];
	arr[] = init;
	return typeof( return )( arr );
}

unittest {
	auto v = vector!int( 8 );
	auto w = vector!int( 16 );
	
	assert( v.size() == 8 );
	assert( w.size() == 16 );
	
	v[0] = 1;
	assert( v[1] == 1 );


struct VectorView( T, int stride_ = 1 ) {
	/** Ensure positive, non-zero stride. */
	static assert( stride_ > 0 );
	
	enum stride = stride_;
	
	/** To allow the vector to be used with expression templates */
	mixin LiteralExpression!( ExpressionKind.VectorLiteral, vectorClosure );
	
	/** Wrap a VectorView around the given array */
	this( T[] a ) pure nothrow
	in {
		static assert( a.length % stride_ == 0 );
	} body {
		array = a;
	}
	
	/** Number of elements in vector */
	size_t size() pure nothrow {
		return a.length / stride_;
	}
	
	/** Return the element at position i */
	T opIndex( size_t i ) pure nothrow
	in {
		assert( i * stride < array.length );
	} body {
		return a[ i * stride ];
	}
	
	/** Assign a value to the element at position i */
	T opIndexAssign( T value, size_t i ) pure nothrow
	in {
		assert( i * stride < array.length );
	} body {
		a[ i * stride ] = value;
	}
	
	/** Perform an assignment operation on the element at position i */
	T opIndexOpAssign( string op )( T value, size_t i ) pure nothrow
	in {
		assert( i * stride < array.length );
	} body {
		mixin("a[ i * stride ]" ~ op ~ "value");
	}
	
	/** The enclosed array */
	T[] array;
}
