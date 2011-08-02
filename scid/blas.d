module scid.blas;

import std.math;
import scid.common.meta;
import scid.common.traits;
import std.traits;
import std.algorithm;
import std.conv, std.string;
import std.ascii;
// debug = blasCalls;
version = nodeps;

debug( blasCalls ) {
	import std.stdio;
	import scid.internal.assertmessages;
}


version( nodeps ) {
	private enum bool forceNaive = true;
} else {
	private enum bool forceNaive = false;
	static import scid.bindings.blas.dblas;
	private alias scid.bindings.blas.dblas blas_;
}

struct blas {
	static void swap( T )( int n, T* x, int incx, T* y, int incy ) {
		debug( blasCalls )
			writef( "swap( %s, %s ) ", stridedToString( x, n, incx ), stridedToString( y, n, incy ) );
		
		static if( isFortranType!T && !forceNaive )
			blas_.swap( n, x, incx, y, incy );
		else
			naive_.swap( n, x, incx, y, incy );
		
		debug( blasCalls )
			writefln( "=> ( %s, %s )", stridedToString( x, n, incx ), stridedToString( y, n, incy ) );
	}
	
	static void scal( T )( int n, T alpha, T* x, int incx ) {
		debug( blasCalls )
			write( "scal( ", alpha, ", ", stridedToString(x, n, incx), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			blas_.scal( n, alpha, x, incx );
		else
			naive_.scal( n, alpha, x, incx );
		
		debug( blasCalls )
			writeln( stridedToString(x, n, incx) );
	}
	
	static void copy( T )( uint n, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls )
			write( "copy( ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			blas_.copy( n, x, incx, y, incy );
		else
			naive_.copy( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString(y, n, incy) );
	}
	
	static void axpy( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls )
			write( "axpy( ", alpha, ", ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			blas_.axpy( n, alpha, x, incx, y, incy );
		else
			naive_.axpy( n, alpha, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString(y, n, incy) );
	}
	
	static T dot( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( !isComplex!T ) {
		debug( blasCalls )
			write( "dot( ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			auto r = blas_.dot( n, x, incx, y, incy );
		else
			auto r = naive_.dot( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T dotu( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "dotu( ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			auto r = blas_.dotu( n, x, incx, y, incy );
		else
			auto r = naive_.dot!true( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T dotc( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "dotc( ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			auto r = blas_.dotc( n, x, incx, y, incy );
		else
			auto r = naive_.dotc( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T nrm2( T )( int n, const(T)* x, int incx ) {
		debug( blasCalls )
			write( "nrm2( ", stridedToString(x, n, incx), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			auto r = blas_.nrm2( n, x, incx );
		else
			auto r = naive_.nrm2( n, x, incx );
		
		debug( blasCalls )
			writeln( r );
	}
	
	static void gemv( char trans, T )( int m, int n, T alpha, const(T)* a, int lda, const(T)* x, int incx, T beta, T *y, int incy ) {
		debug( blasCalls )
			writef( "gemv( %s, %s ) ", matrixToString(trans,m,n,a, lda), stridedToString( x, n, incx ) );
		
		static if( isFortranType!T && !forceNaive )
			blas_.gemv( trans, m, n, alpha, a, lda, x, incx, beta, y, incy );
		else
			naive_.gemv!trans( m, n, alpha, a, lda, x, incx, beta, y, incy );
		
		debug( blasCalls )
			writeln("=> ", stridedToString( x, n, incx ) );
	}
	
	static void trmv( char uplo, char trans, char diag, T )( int n, const(T)* a, int lda, T* x, int incx ) {
		debug( blasCalls )
			writef( "trmv( %s, %s, %s, %s ) ", uplo, diag, matrixToString(trans,n,n,a, lda), stridedToString( x, n, incx ) );
		
		static if( isFortranType!T && !forceNaive )
			blas_.trmv( uplo, trans, diag, n, a, lda, x, incx );
		else
			naive_.trmv!( uplo, trans, diag )( n, a, lda, x, incx );
		
		debug( blasCalls )
			writeln( "=> ", stridedToString( x, n, incx ) );
	}
	
	static void sbmv( char uplo, T )( int n, int k, T alpha, const(T)* a, int lda, const(T)* x, int incx, T beta, T *y, int incy ) {
		debug( blasCalls )
			write( "sbmv( ", uplo, ", ", n, ", ", k, ", ", alpha, ", ", stridedToString(a,n,1), ", ", stridedToString(x,n,incx), ", ", stridedToString(y,n,incy), " ) => " );
		
		static if( isFortranType!T && !forceNaive )
			blas_.sbmv( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy );
		else
			static assert( false );
			
		debug( blasCalls )
			writeln( stridedToString(y,n,incy) );
	}
	
	static void gemm( char transa, char transb, T )( int m, int n, int k, T alpha, const(T)* a, int lda, const(T) *b, int ldb, T beta, T *c, int ldc ) {
		debug( blasCalls )
			writef( "gemm( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s ) =>",
				  transa, transb, m, n, k, alpha, matrixToString( transa, m, k, a, lda ),
				  matrixToString( transb, k, n, b, ldb ), beta, matrixToString( 'N', m, n, c, ldc ) );
		
		static if( isFortranType!T && !forceNaive )
			blas_.gemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );
		else
			static assert( false, "There is no naive implementation of gemm available." );
		
		debug( blasCalls )
			writeln( matrixToString( 'N', m, n, c, ldc ) );
	}
	
	// Level 3
	
	static void trsm( char side, char uplo, char transa, char diag, T )( int m, int n, T alpha, const(T)* a, int lda, T* b, int ldb ) {
		debug( blasCalls )
			writef( "trsm( %s, %s, %s, %s, %s ) ", side, uplo, diag, matrixToString(transa, m, n, a, lda), matrixToString( 'n', m, n, b, ldb ) );
		
		static if( isFortranType!T && !forceNaive )
			blas_.trsm( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );
		else
			naive_.trsm!( side, uplo, transa, diag )( m, n, alpha, a, lda, b, ldb );
			//static assert( false, "There is no naive implementation of trsm available." );
		
		debug( blasCalls )
			writeln( matrixToString( 'n', m, n, b, ldb ) );
	}
	
	// Extended BLAS, stuff I needed and wasn't implemented by BLAS.
	
	// x := conj( x )
	static T xconj( T )( T x ) if( isComplex!T ) {
		return naive_.xconj( x );
	}
	
	// x := x.H
	static void xcopyc( T )( int n, T* x, int incx ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xcopyc( ", stridedToString(x, n, incx), " ) => " );
		
		naive_.xcopyc( n, x, incx );
		
		debug( blasCalls )
			writeln( stridedToString(x, n, incx) );
	}
	
	// y := x.H
	static void xcopyc( T )( int n, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xcopyc( ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		naive_.xcopyc( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString(y, n, incy) );
	}
	
	// y := alpha*x.H + y
	static void xaxpyc( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xaxpyc( ", alpha, ", ", stridedToString(x, n, incx), ", ", stridedToString(y, n, incy), " ) => " );
		
		naive_.xcopyc( n, alpha, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString(y, n, incy) );
	}
	
	// General matrix copy
	// B := A    or
	// B := A.T  or
	// B := A.H, for A and B mxn matrices
	static void xgecopy( char transA, T )( int m, int n, const(T)* a, int lda, T* b, int ldb ) {
		debug( blasCalls ) {
			writeln();
			writeln( "xgecopy( ", matrixToString(transA,m,n,a,lda), ", ", matrixToString('N',m,n,a,ldb), " ) => ..." );
		}
		
		naive_.xgecopy!transA( m, n, a, lda, b, ldb );
		
		debug( blasCalls ) {
			writeln( "/xgecopy()" );
			writeln();
		}
	}
	
	
	// General matrix copy
	// B := conj(A)   or
	// B := conj(A.T), for A and B mxn matrices
	static void xgecopyc( T, char transA )( int m, int n, const(T)* a, int lda, T* b, int ldb ) {
		debug( blasCalls ) {
			writeln();
			writeln( "xgecopyc( ", matrixToString(transA,m,n,a,lda), ", ", matrixToString('N',m,n,a,ldb), " ) => ..." );
		}	
		
		naive_.xgecopy!( transA, true )(  m, n, a, lda, b, ldb );
		
		debug( blasCalls ) {
			writeln( "/xgecopyc()" );
			writeln();
		}
	}
	
}

private struct naive_ {
	private static void reportNaive_() {
		debug( blasCalls )
			write( "<n> " );
	}
	
	private static bool checkMatrix( T )( char trans, int m, int n, const T* a, int lda ) {
		if( trans == 'T' || trans == 'C' )
			std.algorithm.swap( m, n );
		assert(trans == 'T' || (trans == 'C' && isComplex!(Unqual!T))  || trans == 'N',
			   "Invalid transposition character '" ~ trans ~ "' for '" ~ T.stringof ~ ".");
		assert( a != null, "Null matrix." );
		assert( lda >= m, format("Leading dimension less than minor dimension: %d vs %d",lda,m) );
		return true;
	}
	
	private static bool checkVector( T )( T* x, int incx ) {
		assert( incx != 0 );
		assert( x != null );
		return true;
	}
	
	// LEVEL 1
	static void swap( T )( int n, T* x, int incx, T* y, int incy ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n )
			return;
		
		assert( checkVector( x, incx )	);
		assert( checkVector( y, incy )	);
		
		if( incx == 1 && incy == 1 ) {
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; *x = *y; *y = aux;
				++ x; ++ y;
			}
		} else {
			n *= incx;
			auto xe = x + n * incx;
			T aux;
			while( x != xe ) {
				aux = *x; *x = *y; *y = aux;
				x += incx; y += incy;
			}
		}
	}
	
	static void scal( T )( int n, T alpha, T* x, int incx ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n || alpha == One!T )
			return;
		
		assert( checkVector( x, incx )	);
		if( incx == 1 ) {
			if( !alpha )
				x[ 0 .. n ] = Zero!T;
			else
				x[ 0 .. n ] *= alpha;
		} else if( alpha ) {
			n *= incx;
			auto xe	= x + n;
			while( x != xe ) {
				(*x) *= alpha;
				x += incx;
			}
		} else {
			n *= incx;
			auto xe	= x + n;
			while( x != xe ) {
				(*x) = Zero!T;
				x += incx;
			}
		}
	}
	
	static void copy( T )( int n, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n )
			return;
		
		assert( checkVector( x, incx ) );
		assert( checkVector( y, incy ) );
		
		if( incx == 1 && incy == 1 )
			y[ 0 .. n ] = x[ 0 .. n ];
		else {
			n *= incx;
			auto xe = x + n;
			while( x != xe ) {
				*y = *x;
				x += incx; y += incy;
			}
		}
	}
	
	static void axpy( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls ) reportNaive_();
		
		if( alpha == Zero!T || !n )
			return;
		
		assert( checkVector( y, incy )	);
		assert( checkVector( x, incx )	);
		
		if( incx == 1 && incy == 1 ) {
			y[ 0 .. n ] += x[ 0 .. n ] * alpha;
		} else {
			n *= incx;
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; aux *= alpha;
				*y += aux;
				x += incx; y += incy;
			}
		}
	}
	
	static T dot( bool forceComplex = false, T )( int n, const(T)* x, int incx, const(T)* y, int incy )
				if( !isComplex!T || forceComplex ) {
		debug( blasCalls ) reportNaive_();
				
		assert( checkVector( y, incy )	);
		assert( checkVector( x, incx )	);
		
		T r = Zero!T;
		if( incx == 1 && incy == 1 ) {
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; aux *= *y;
				r += aux;
				++x; ++y;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; aux *= *y;
				r += aux;
				x += incx; y += incy;
			}
		}
		
		return r;
	}
	
	template dotu( T ) {
		alias dot!(true,T) dotu;
	}
	
	static T dotc( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( isComplex!T ) {
		debug( blasCalls ) reportNaive_();
		
		assert( checkVector( y, incy )	);
		assert( checkVector( x, incx )	);
		
		T r = Zero!T;
		if( incx == 1 && incy == 1 ) {
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; aux = conj( aux );
				aux *= *y;
				r += aux;
				++x; ++y;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			T aux;
			while( x != xe ) {
				aux = *x; aux = conj( aux );
				aux *= *y;
				r += aux;
				x += incx; y += incy;
			}
		}
		
		return r;
	}
	
	static T nrm2( T )( int n, const(T)* x, int incx ) {
		debug( blasCalls ) reportNaive_();
		
		assert( checkVector( x, incx )	);
		T r = Zero!T;
		T aux;
		if( incx == 1 ) {
			auto xe = x + n;
			while( x != xe ) {
				aux = *x; aux *= aux;
				r += aux;
				++ x;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			while( x != xe ) {
				aux = *x; aux *= aux;
				r += aux;
				++ x;
			}
		}
		
		aux = sqrt( r );
		return aux;
	}
	
	// LEVEL 2
	static void gemv( char trans_, T )( int leny, int lenx, T alpha, const(T)* a, int lda, const(T)* x, int incx, T beta, T *y, int incy ) {
		enum trans = cast(char)toUpper( trans_ );
		
		debug( blasCalls ) reportNaive_();
		if( !lenx || !alpha ) {
			if( trans == 'N' )
				scal( leny, beta, y, incy );
			else
				scal( lenx, beta, y, incy );
			return;
		}
		
		assert( checkVector( y, incy ) );
		assert( checkVector( x, incx ) );
		assert( checkMatrix( trans, leny, lenx, a, lda ) );
		
		static if( trans == 'N' ) {
			foreach( i ; 0 .. leny ) {
				T temp = Zero!T;
				foreach( j ; 0 .. lenx )
					 temp += a[ j * lda + i ] * x[ j * incx ];
				temp *= alpha;
				
				y[ i * incy ] *= beta;
				y[ i * incy ] += temp;
			}
		} else {
			foreach( i ; 0 .. lenx ) {
				T temp = Zero!T;
				foreach( j ; 0 .. leny )
					 temp += a[ i * lda + j ] * x[ j * incx ];
				temp *= alpha;
				
				y[ i * incy ] *= beta;
				y[ i * incy ] += temp;
			}
		}
	}
	
	static void trmv( char uplo_, char trans_, char diag_, T)( int n, const(T)* a, int lda, T* x, int incx ) {
		enum uplo = cast(char)toUpper( uplo_ );
		enum trans = cast(char)toUpper( trans_ );
		enum diag = cast(char)toUpper( diag_ );
		
		static assert( uplo == 'U' || uplo == 'L' );
		static assert( trans == 'N' || trans == 'T' || ( trans == 'C' && isComplex!T ) );
		static assert( diag == 'U' || diag == 'N' );
		assert( n >= 0 );
		assert( lda >= max(1, n) );
		assert( incx != 0 );
		
		if( n == 0 )
			return;
		
		static if( (uplo == 'U') ^ (trans != 'N') ) {
			// Upper No-Transpose or Lower Transposed
			foreach( j ; 0 .. n ) {
				T xj = x[ j * incx ];
				static if( trans == 'N' ) {
					if( xj != Zero!T ) {
						foreach( i ; 0 .. j )
							x[ i * incx ] += xj * a[ j * lda + i ];
						
						static if( diag == 'N' )
							x[ j * incx ] *= a[ j * lda + j ];
					}
				} else {
					static if( diag == 'N' )
						xj *= a[ j * lda + j ];
						
					foreach( i ; j + 1 .. n ) {
						static if( trans == 'C' )
							xj += blas.xconj( a[ j * lda + i ] ) * x[ i * incx ];
						else
							xj += a[ j * lda + i ] * x[ i * incx ];
					}
					
					x[ j * incx ] = xj;
				}
			 }
		} else {
			// Lower No-Transpose or Upper Transposed
			for( int j = n - 1; j >= 0 ; -- j ) {
				T xj = x[ j * incx ];
				static if( trans == 'N' ) {
					if( xj != Zero!T ) {
						for( int i = n - 1 ; i > j ; -- i )
							x[ i * incx ] += xj * a[ j * lda + i ];
						
						static if( diag == 'N' )
							x[ j * incx ] *= a[ j * lda + j ];
					}
				} else {
					static if( diag == 'N' )
						xj *= a[ j * lda + j ];
					
					for( int i = j - 1 ; i >= 0 ; -- i )
						xj += a[ j * lda + i ] * x[ i * incx ];
					
					x[ j * incx ] = xj;
				}
			}
		}
	}
	
	// Level 3
	
	static void trsm( char side_, char uplo_, char trans_, char diag_,T )( int m, int n, T alpha, const(T)* a, int lda, T* b, int ldb ) {
		enum side = cast(char)toUpper(side_);
		enum uplo = cast(char)toUpper(uplo_);
		enum trans = cast(char)toUpper(trans_);
		enum diag = cast(char)toUpper(diag_);
		
		static assert( false, "No naive implementation of trsm available, sorry." );
	}

	// Extended
	// x := conj( x )
	static T xconj( T )( T x ) if( isComplex!(Unqual!T) ) {
		static if( is( T E : Complex!E ) )
			return x.conj;
		else
			return conj( x );
	}
	
	// x := x.H
	static void xcopyc( T )( int n, T* x, int incx ) if( isComplex!T ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n )
			return;
		
		assert( checkVector( x, incx ) );
		
		if( incx == 1 ) {
			auto xe = x + n;
			while( x != xe ) {
				*x = xconj( *x );
				++ x;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			while( x != xe ) {
				*x = xconj( *x );
				x += incx;
			}
		}
	}
	
	// y := x.H
	static void xcopyc( T )( int n, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n )
			return;
		
		assert( checkVector( x, incx ) );
		
		if( incx == 1 && incy == 1 ) {
			auto xe = x + n;
			while( x != xe ) {
				*y = xconj( *x );
				++ x; ++ y;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			while( x != xe ) {
				*y = xconj( *x );
				x += incx; y += incy;
			}
		}
	}
	
	// y := alpha*x.H + y
	static void xaxpyc( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls ) reportNaive_();
		
		if( !n || alpha == Zero!T )
			return;
		
		assert( checkVector( x, incx ) );
		
		T aux;
		if( incx == 1 && incy == 1 ) {
			auto xe = x + n;
			while( x != xe ) {
				aux = *x; aux = xconj( aux ) * alpha;
				*y += aux;
				++ x; ++ y;
			}
		} else {
			n *= incx;
			auto xe = x + n;
			while( x != xe ) {
				aux = *x; aux = xconj( aux ) * alpha;
				*y += aux;
				x += incx; y += incy;
			}
		}
	}
	
	// General matrix copy
	// B := A    or
	// B := A.T  or
	// B := A.H    , for A and B mxn matrices
	static void xgecopy( char transA_, bool forceConjugate=false, T )( int m, int n, const(T)* a, int lda, T* b, int ldb ) if( isComplex!T || !forceConjugate ) {
		//debug( blasCalls ) reportNaive_();
		
		enum transA = cast(char)toUpper( transA_ );
		
		if( !m || !n )
			return;
		
		assert( checkMatrix( transA, m, n, a, lda ) );
		assert( checkMatrix( 'N', m, n, b, ldb ) );
		
		if( transA == 'N' ) {
			// if a is not transposed
			if( lda == m && ldb == m ) {
				// if there are no gaps in their memory, just copy
				// eveything
				static if( forceConjugate )
					blas.xcopyc( m * n, a, 1, b, 1 );
				else
					blas.copy( m *n, a, 1, b, 1 );
				
			} else {
				// if there are gaps, copy column-by-column
				n *= ldb;
				auto be = b + n;
				while( b != be ) {
					static if( forceConjugate )
						blas.xcopyc( m, a, 1, b, 1 );
					else
						blas.copy( m, a, 1, b, 1 );
					a += lda; b += ldb;
				}
			}
		} else if( transA == 'T' ) {
			// if a is transposed, copy a row-by-row to b column-by-column
			n *= ldb;
			auto be = b + n;
			while( b != be ) {
				blas.copy( m, a, lda, b, 1 );
				++ a; b += ldb;
			}
		} else {
			static if( !isComplex!(Unqual!T) ) {
				assert( false,
					"'" ~ transA ~ "', invalid value for 'transA' in matrix of type '" ~ T.stringof ~ "' copy." );
			} else {
				// assume transA == 'C'
				n *= ldb;
				auto be = b + n;
				while( b != be ) {
					blas.xcopyc( m, a, lda, b, 1 );
					++ a; b += ldb;
				}
			}
		}
	}
	
	// General matrix axpy
	// B := alpha*A   + B  or
	// B := alpha*A.T + B  or
	// B := alpha*A.H + B, for A and B mxn matrices
	static void xgeaxpy( char transA_, bool forceConjugate=false,  T )( int m, int n, T alpha, const(T)* a, int lda, T* b, int ldb ) if( isComplex!T || !forceConjugate ) {
		//debug( blasCalls ) reportNaive_();
		
		enum transA = cast(char)toUpper( transA_ );
		
		if( !m || !n )
			return;
		
		assert( checkMatrix( transA, m, n, a, lda ) );
		assert( checkMatrix( 'N', m, n, b, ldb ) );
		
		if( transA == 'N' ) {
			// if a is not transposed
			if( lda == m && ldb == m ) {
				// if there are no gaps in their memory, just copy
				// eveything
				static if( forceConjugate )
					blas.xaxpyc( m * n, alpha, a, 1, b, 1 );
				else
					blas.axpy( m * n, alpha, a, 1, b, 1 );
				
			} else {
				// if there are gaps, copy column-by-column
				n *= ldb;
				auto be = b + n;
				while( b != be ) {
					static if( forceConjugate )
						blas.xaxpyc( m, alpha, a, 1, b, 1 );
					else
						blas.axpy( m, alpha, a, 1, b, 1 );
					a += lda; b += ldb;
				}
			}
		} else if( transA == 'T' ) {
			// if a is transposed, copy a row-by-row to b column-by-column
			n *= ldb;
			auto be = b + n;
			while( b != be ) {
				blas.axpy( m, alpha, a, lda, b, 1 );
				++ a; b += ldb;
			}
		} else {
			static if( !isComplex!(Unqual!T) ) {
				assert( false,
					"'" ~ transA ~ "', invalid value for 'transA' in matrix of type '" ~ T.stringof ~ "' copy." );
			} else {
				// assume transA == 'C'
				n *= ldb;
				auto be = b + n;
				while( b != be ) {
					blas.xaxpyc( m, alpha, a, lda, b, 1 );
					++ a; b += ldb;
				}
			}
		}
	}
	
}
