/** Fallback routines for evaluating primitive operations.
    
    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ops.fallback;

import scid.matvec;
import scid.common.traits, scid.common.meta;
import std.conv;
import scid.ops.expression;
import scid.ops.common;
import scid.ops.eval;
import std.string, std.algorithm;
import scid.internal.regionallocator;

//debug = fallbackCalls;

debug( fallbackCalls ) {
	import std.stdio : write, writeln, stdout;
}

private {
	string toString(T)( T[]   v ) @property if(isScalar!T) { return to!string( v ); }
	string toString(T)( T[][] m ) @property if(isScalar!T) { return to!string( m ); }
}

/** Matrix product. Similiar to the mm BLAS routines, performs matrix-vector products to compute the result. */
void fallbackMatrixProduct( Transpose transA, Transpose transB, A, B, E, Dest )
		( E alpha, auto ref A a, auto ref B b, E beta, auto ref Dest dest ) {
	size_t m = transA ? a.columns : a.rows;
	size_t n = transB ? b.rows : b.columns;
	
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
		E dot = rowColumnDot!( transA, transB )(a, i, b, j);	
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
		debug( fallbackCalls ) write( "fbscal( ", alpha, ", ", dest.toString, " ) => " );
		for( size_t i = 0; i < l ; ++ i )
			dest[ i ] *= alpha;
		debug( fallbackCalls ) writeln( dest.toString() );	
	} else static if( destClosure == Closure.Matrix ) {
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
		debug( fallbackCalls ) { write( "fbcopy( ", source.toString(), ", ", dest.toString(), " ) => " ); }
		assert( source.length == dest.length, "fbcopy: Length mismatch in fallback vector assignment." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i )
			dest[ i ] = source[ i ];
		debug( fallbackCalls ) writeln( dest.toString() );
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

/** Performs a scaled addition on a matrix or vector (i.e. axpy in BLAS terms, a*alpha + b). Vectors are processed
    manually, while matrices fall back on vector routines.
*/
void fallbackScaledAddition( Transpose srcTrans, Scalar, Source, Dest )
		( Scalar alpha, auto ref Source source, auto ref Dest dest ) {
	enum destClosure = closureOf!Dest;
	alias srcTrans tr;
			
	static if( isVectorClosure!destClosure ) {	
		// fallback vector axpy
		auto alphaValue = eval( alpha );
		debug( fallbackCalls ) {
			auto dummy = dest.data;
			write( "fbaxpy",(tr?".t":""), "( ",  alphaValue, ", ", source.toString(), ", ", dest.toString(), " ) => " );
		}
		assert( source.length == dest.length, "Length mismatch in fallback vector addition." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i ) {
			static if( tr && isComplex!(BaseElementType!Source) ) {
				dest[ i ] += gconj(source[ i ]) * alphaValue;
			} else {
				dest[ i ] += source[ i ] * alphaValue;
			}
			
		}
		debug( fallbackCalls ) writeln( dest.toString() );
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

/** Solves a linear problem (i.e. sv/sm in BLAS terms). */
void fallbackSolve( Transpose transM, Side side, M, Dest )( auto ref M mat, auto ref Dest dest ) {
	auto alloc = newRegionAllocator();
	
	static if( transM )
		auto matInverse = eval( inv( mat.t ), alloc );
	else
		auto matInverse = eval( inv( mat ), alloc );
	
	static if( side == Side.Left )
		dest[] = eval( matInverse * dest, alloc );
	else
		dest[] = eval( dest * matInverse, alloc );
}

/** Matrix inversion. */
void fallbackInverse( M )( ref M mat ) {
	alias BaseElementType!M T;
	
	assert( mat.rows == mat.columns,
		   format( "Cannot invert non-square matrix with dimensions %d and %d.", mat.rows, mat.columns ) );
	
	size_t n = mat.rows;
	
	if( n == 1 ) {
		assert( mat[ 0, 0 ] != Zero!T, "Singular matrix in inversion" );
		mat[ 0, 0 ] = One!T / mat[ 0, 0 ];
		return;
	}
	
	T get( size_t i, size_t j ) {
		//static if( storageOrderOf!M == StorageOrder.ColumnMajor )
			return mat[ j, i ];
		//else
			//return mat[ i, j ];
	}
	
	void set(string op = "/")( T x, size_t i, size_t j ) {
		//static if( storageOrderOf!M == StorageOrder.ColumnMajor )
			mixin( "mat[ j, i ] " ~ op ~ "= x;" );
		//else
			//mixin( "mat[ i, j ] " ~ op ~ "= x;" );
	}
	
	// normalize first row
	auto normWith = get(0,0);
	foreach( i ; 1 .. n )
		set!"/"( normWith, i, n );
	/*
	// LU Factorization
	foreach( i ; 1 .. n ) {
		// Column of L
		foreach( j ; i .. n ) {
			T sum = Zero!T;
			foreach( k ; 0 .. i )
				sum += get( i, k ) * get( k, i );
			
			set!"-"( sum, j, i );
		}
	
		if( i == n - 1 )
			continue;
		
		// Row of U
		foreach( j ; i+1 .. n ) {
			T sum = Zero!T;
			foreach( k ; i .. j )
				sum += get( i, k ) * get( k, j );
			set( (get(i,j) - sum ) / get( i, i ), i, j );
		}
	}
	
	// Invert L
	foreach( i ; 0 .. n ) {
		foreach( j ; i .. n ) {
			T x = One!T;
			if( i != j ) {
				x = Zero!T;
				foreach( k ; i .. j )
					x -= get(j,k) * get(k,i);
			}
			set( x / get(j,j), j, i );
		}
	}
	
	// Invert U
	foreach( i ; 0 .. n ) {
		foreach( j ; i .. n ) {
			if ( i == j ) continue;
			T sum = Zero!T;
			foreach ( k ; i .. j )
				sum +=get(k,j)*( (i==k) ? One!T : get(i,k) );
			set( -sum, i, j );
		}
	}
	
	// final inversion
	foreach( i ; 0 .. n ) {
		foreach( j ; 0 .. n ) {
			T sum = Zero!T;
			foreach( k ; max(i,j) .. n )
				sum += ( (j==k) ? 1.0 : get( j, k ) ) * get( k, i );
			set( sum, j, i );
		}
	}*/
}