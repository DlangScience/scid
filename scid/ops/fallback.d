/** Fallback routines for evaluating primitive operations.
    
    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ops.fallback;

import scid.matvec;
import scid.common.traits, scid.common.meta;
import std.conv  : to;
import scid.ops.expression;
import scid.ops.common;
import scid.ops.eval;

//debug = fallbackCalls;

debug( fallbackCalls ) {
	import std.stdio : write, writeln, stdout;
	
}

/** Matrix product. Similiar to the mm BLAS routines, performs matrix-vector products to compute the result. */
void fallbackMatrixProduct( Transpose transA, Transpose transB, A, B, E, Dest )
		( E alpha, auto ref A a, auto ref B b, E beta, auto ref Dest dest ) {
	size_t m = transA ? a.rows    : a.columns;
	size_t n = transB ? a.columns : a.rows;
	
	if( !alpha ) {
		if( !beta ) dest.resize( m, n );
		else        evalScaling( beta, dest );
		return;
	} 
	
	if( !beta )
		dest.resize( m, n, null );
	else
		assert( dest.rows == m && dest.columns == n, "Matrix size mismatch in matrix addition." );
		
		
	void doElement( size_t i, size_t j ) {
		E dot;
			
		static if( !transA && !transB )     dot = evalDot!( transA, transB )( a.row( i ),    b.column( j ) );
		else static if( !transA && transB ) dot = evalDot!( transA, transB )( a.row( i ),    b.row( j )    );
		else static if( transA && !transB ) dot = evalDot!( transA, transB )( a.column( i ), b.column( j ) );
		else static if( transA && transB )  dot = evalDot!( transA, transB )( a.column( i ), b.row( j )    );
			
		dot *= alpha; dot += beta * dest[i, j];
		dest[ i, j ] = dot;	
	}
		
	static if( Dest.storageOrder == StorageOrder.RowMajor ) {
		foreach( i ; 0 .. m ) {
			foreach( j ; 0 .. n ) doElement( i, j );
		}
	} else {
		foreach( j ; 0 .. n ) {
			foreach( i ; 0 .. m ) doElement( i, j );
		}
	}
}

/** Manually compute the dot product of two vectors. */
auto fallbackDot( Transpose transA, Transpose transB, A, B )( ref A a, ref B b ) {
	alias BaseElementType!A T;
	
	assert( a.length == b.length, "fbdot: Vector length mismatch in fallback dot product." );
	debug( fallbackCalls )  {
		string aStr = a.toString() ~ (transA ? ".t" : "");
		string bStr = b.toString() ~ (transB ? ".t" : "");
		write( "fb_dot( ", aStr , ", ", bStr, " ) => " );
	}
	int n = a.length;
	auto r = Zero!T;
		
	foreach( i ; 0 .. n ) {
		static if( !isComplex!T || (transA == transB) )
			r += a[ i ] * b[ i ];
		else static if(  transA && !transB )
			r += gconj( a[ i ] ) * b[ i ];
		else 
			r += a[ i ] * gconj( b[ i ] );
	}
		
	static if( transA && transB && isComplex!T )
		r = gconj( r );
		
	debug( fallbackCalls ) writeln( r );
	return r;
}

/** Computes a matrix-vector product by computing dot products. */
void fallbackMatrixVectorProduct( Transpose transM, Transpose transV, M, V, Dest, E )( E a, ref M m, ref V v, E b, ref Dest dest ) {
	auto rows = transM ? m.columns : m.rows;
	auto cols = transM ? m.rows    : m.columns;
	
	assert( cols == v.length, "fbmv: Dimensions mismatch in matrix-vector product - " ~ to!string( cols ) ~ " vs. " ~ to!string( v.length ) );
	
	if( b )
		assert( rows == dest.length, "fbmv: Dimensions mismatch in vector-vector addition - " ~ to!string( rows ) ~ " vs. " ~ to!string( dest.length ) );
	else
		dest.resize( rows );
	
	assert( rows == dest.length );
	
	foreach( i ; 0 .. rows ) {
		static if( transM )
			dest[ i ] = evalDot!(transM,transV)( m.column( i ), v ) + b * dest[ i ];
		else
			dest[ i ] = evalDot!(transM,transV)( m.row( i ), v ) + b * dest[ i ];
	}
}

/** Scales a vector or a matrix. Vector scaling is performed manually while matrix scaling falls back to vector scaling. */
void fallbackScaling( Scalar, Dest )( Scalar alpha, auto ref Dest dest ) {
	enum destClosure = closureOf!Dest;
	static if( isVectorClosure!destClosure ) {
		auto l = dest.length;
		debug( fallbackCalls ) write( "fbscal( ", alpha, ", ", to!string(dest), " ) => " );
		for( size_t i = 0; i < l ; ++ i )
			dest[ i ] *= alpha;
		debug( fallbackCalls ) writeln( to!string(dest) );	
	} else static if( closure == Closure.Matrix ) {
		static if( storageOrderOf!Dest == StorageOrder.RowMajor ) {
			size_t n = dest.rows;
			for( auto i = 0; i < n ; ++ i )
				evalScaling( alpha, dest.row( i ) );
		} else {
			size_t n = dest.columns;
			for( auto i = 0; i < n ; ++ i )
				evalScaling( alpha, dest.column( i ) );
		}
	} else static assert( false );
}

/** Copies a matrix or a vector. Vector copying is performed manually, while matrix copying falls back to vector
    copying.
*/
void fallbackCopy( Transpose srcTrans, Source, Dest )( auto ref Source source, auto ref Dest dest ) {
	alias srcTrans tr;
	enum destClosure = closureOf!Dest;
	
	static if( isVectorClosure!destClosure ) {
		// fallback vector copy
		debug( fallbackCalls ) { write( "fbcopy( ", to!string(source), ", ", to!string(dest), " ) => " ); }
		assert( source.length == dest.length, "fbcopy: Length mismatch in fallback vector assignment." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i )
			dest[ i ] = source[ i ];
		debug( fallbackCalls ) writeln( to!string(dest) );
	} else static if( destClosure == Closure.Matrix ) {
		// fallback matrix copy - copy major-by-major
		size_t m = tr ? source.columns : source.rows;
		size_t n = tr ? source.rows    : source.columns;
		assert( m == dest.rows && n == dest.columns, "fbcopy: Matrix size mismatch in fallback matrix assignment.");
		
		static if( storageOrderOf!Dest == StorageOrder.RowMajor ) {	
			for( auto i = 0; i < m ; ++ i )
				static if( tr )
					evalCopy!tr( source.column(i), dest.row( i ) );
				else
					evalCopy( source.row(i),    dest.row( i ) );
		} else {
			for( auto i = 0; i < n ; ++ i )
				static if( tr )
					evalCopy!tr( source.row(i), dest.column( i ) );
				else
					evalCopy( source.column(i), dest.column( i ) );
		}
	} else static assert( false );
}

/** Performs a scaled addition on a matrix or vector (i.e. axpy in BLAS terms, a*alpha + b). Vectors are processed,
    manually, while matrices fall back on vector routines.
*/
void fallbackScaledAddition( Transpose srcTrans, Scalar, Source, Dest )
		( Scalar scalar, auto ref Source src, auto ref Dest dest ) {
	enum destClosure = closureOf!Dest;		
			
	static if( isVectorClosure!destClosure ) {	
		// fallback vector axpy
		auto alphaValue = eval( alpha );
		debug( fallbackCalls ) {
			auto dummy = dest.data;
			write( "fbaxpy( ", alphaValue, ", ", to!string(source), ", ", to!string(dest), " ) => " );
		}
		assert( source.length == dest.length, "Length mismatch in fallback vector addition." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i )
			dest[ i ] += source[ i ] * alphaValue;	
		debug( fallbackCalls ) writeln( to!string(dest) );
	} else static if( destClosure == Closure.Matrix ) {
		// fallback matrix axpy - axpy major-by-major
		auto alphaValue = eval( alpha );
		size_t m = tr ? source.columns : source.rows;
		size_t n = tr ? source.rows    : source.columns;
		assert( m == dest.rows && n == dest.columns, "Matrix size mismatch in matrix addition.");
		
		static if( dest.storageOrder == StorageOrder.RowMajor ) {	
			for( auto i = 0; i < m ; ++ i )
				static if( tr ) evalScaledAddition!tr( alpha, source.column(i), dest.row( i ) );
				else            evalScaledAddition( alpha, source.row(i),    dest.row( i ) );
		} else {
			for( auto i = 0; i < n ; ++ i )
				static if( tr ) evalScaledAddition!tr( alpha, source.row(i),    dest.column( i ) );
				else            evalScaledAddition( alpha, source.column(i), dest.column( i ) );
		}
	} else static assert( false );		
}