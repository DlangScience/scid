module scid.blas;

import std.math;
import scid.common.meta;
import scid.common.traits;
import std.traits;
import std.algorithm;
import std.conv;
import std.string;
// debug = blasCalls;
// version = nodeps;

debug( blasCalls ) {
	import std.stdio, std.range;
	
	private string stridedToString_( S )( const(S)* ptr, size_t len, size_t stride ) {
		if( len == 0 )
			return "[]";
	
		auto app = appender!string("[");
		app.put( to!string(*ptr) );
		auto e = ptr + len * stride;
		ptr += stride;
		for( ; ptr < e ; ptr += stride ) {
			app.put( ", " );
			app.put( to!string( *ptr ) );
		} 
		app.put(']');
		
		return app.data();
	}
	
	private string matrixToString_( S )( char trans, size_t m, size_t n, const(S)* a, size_t lda )
	in {
		assert( trans == 't' || trans == 'n' || trans == 'c' );
	} body {
		if( m == 0 || n == 0 )
			return "[]";
		
		auto app = appender!string("[");
		if( trans == 'n' ) {
			app.put( stridedToString_( a, n, lda ) );
			auto e = a + m;
			for( ; a < e ; ++a ) {
				app.put( ", " );
				app.put( stridedToString_( a, n, lda ) );
			}
			app.put(']');
		} else {
			app.put( stridedToString_(a, m, 1) );
			auto e = a + n * lda;
			for( ; a < e ; a += lda ) {
				app.put( ", " );
				app.put( stridedToString_( a, m, 1 ) );
			}
			app.put(']');
		}
		
		return app.data();
	}
}


version( nodeps ) {
	private alias naive_ blas_;
} else {
	static import scid.bindings.blas.dblas;
	private alias scid.bindings.blas.dblas blas_;
}

struct blas {
	static void swap( T )( int n, T* x, int incx, T* y, int incy ) {
		static if( isFortranType!T )
			blas_.swap( n, x, incx, y, incy );
		else
			naive_.swap( n, x, incx, y, incy );
	}
	
	static void scal( T )( int n, T alpha, T* x, int incx ) {
		debug( blasCalls )
			write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
		
		static if( isFortranType!T )
			blas_.scal( n, alpha, x, incx );
		else
			naive_.scal( n, alpha, x, incx );
		
		debug( blasCalls )
			writeln( stridedToString_(x, n, incx) );
	}
	
	static void copy( T )( uint n, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls )
			write( "copy( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		static if( isFortranType!T )
			blas_.copy( n, x, incx, y, incy );
		else
			naive_.copy( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString_(y, n, incy) );
	}
	
	static void axpy( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) {
		debug( blasCalls )
			write( "axpy( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		static if( isFortranType!T )
			blas_.axpy( n, alpha, x, incx, y, incy );
		else
			naive_.axpy( n, alpha, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString_(y, n, incy) );
	}
	
	static T dot( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( !isComplex!T ) {
		debug( blasCalls )
			write( "dot( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		static if( isFortranType!T )
			auto r = blas_.dot( n, x, incx, y, incy );
		else
			auto r = naive_.dot( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T dotu( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "dotu( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		static if( isFortranType!T )
			auto r = blas_.dotu( n, x, incx, y, incy );
		else
			auto r = naive_.dot!true( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T dotc( T )( int n, const(T)* x, int incx, const(T)* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "dotc( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		static if( isFortranType!T )
			auto r = blas_.dotc( n, x, incx, y, incy );
		else
			auto r = naive_.dotc( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( r );
		
		return r;
	}
	
	static T nrm2( T )( int n, const(T)* x, int incx ) {
		debug( blasCalls )
			write( "nrm2( ", stridedToString_(x, n, incx), " ) => " );
		
		static if( isFortranType!T )
			auto r = blas_.nrm2( n, x, incx );
		else
			auto r = naive_.nrm2( n, x, incx );
		
		debug( blasCalls )
			writeln( r );
	}
	
	static void gemv( T )( char trans, int m, int n, T alpha, const(T)* a, int lda, const(T)* x, int incx, T beta, T *y, int incy ) {
		static if( isFortranType!T )
			blas_.gemv( trans, m, n, a, lda, x, incx, beta, y, incy );
		else
			naive_.gemv( trans, m, n, a, lda, x, incx, beta, y, incy );
	}
	
	static void sbmv( T )( char uplo, int n, int k, T alpha, const(T)* a, int lda, const(T)* x, int incx, T beta, T *y, int incy ) {
		debug( blasCalls )
			write( "sbmv( ", uplo, ", ", n, ", ", k, ", ", alpha, ", ", stridedToString_(a,n,1), ", ", stridedToString_(x,n,incx), ", ", stridedToString_(y,n,incy), " )" );
		
		static if( isFortranType!T )
			blas_.sbmv( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy );
		else
			static assert( false );
			//naive_.gemv( trans, m, n, a, lda, x, incx, beta, y, incy );
		debug( blasCalls )
			writeln( stridedToString_(y,n,incy) );
	}
	
	// Extended BLAS, stuff I needed and wasn't implemented by BLAS.
	
	// x := conj( x )
	static T xconj( T )( T x ) if( isComplex!T ) {
		return naive_.xconj( x );
	}
	
	// x := x.H
	static void xcopyc( T )( int n, T* x, int incx ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xcopyc( ", stridedToString_(x, n, incx), " ) => " );
		
		naive_.xcopyc( n, x, incx );
		
		debug( blasCalls )
			writeln( stridedToString_(x, n, incx) );
	}
	
	// y := x.H
	static void xcopyc( T )( int n, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xcopyc( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		naive_.xcopyc( n, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString_(y, n, incy) );
	}
	
	// y := alpha*x.H + y
	static void xaxpyc( T )( int n, T alpha, const(T)* x, int incx, T* y, int incy ) if( isComplex!T ) {
		debug( blasCalls )
			write( "xaxpyc( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
		
		naive_.xcopyc( n, alpha, x, incx, y, incy );
		
		debug( blasCalls )
			writeln( stridedToString_(y, n, incy) );
	}
	
	// General matrix copy
	// B := A    or
	// B := A.T  or
	// B := A.H, for A and B mxn matrices
	static void xgecopy( T )( char transA, int m, int n, const(T)* a, int lda, T* b, int ldb ) {
		debug( blasCalls ) {
			writeln();
			writeln( "xgecopy( ", matrixToString_(transA,m,n,a,lda), ", ", matrixToString_('n',m,n,a,ldb), " ) => ..." );
		}
		
		naive_.xgecopy( transA, m, n, a, lda, b, ldb );
		
		debug( blasCalls ) {
			writeln( "/xgecopy()" );
			writeln();
		}
	}
	
	
	// General matrix copy
	// B := conj(A)   or
	// B := conj(A.T), for A and B mxn matrices
	static void xgecopyc( T )( char transA, int m, int n, const(T)* a, int lda, T* b, int ldb ) {
		debug( blasCalls ) {
			writeln();
			writeln( "xgecopyc( ", matrixToString_(transA,m,n,a,lda), ", ", matrixToString_('n',m,n,a,ldb), " ) => ..." );
		}
			
		
		naive_.xgecopy!true( transA, m, n, a, lda, b, ldb );
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
		if( trans == 't' || trans == 'c' )
			std.algorithm.swap( m, n );
		assert(trans == 't' || (trans == 'c' && isComplex!(Unqual!T))  || trans == 'n',
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
		
		if( incx == 1 && incy == 1 )
			std.algorithm.swap( x[0 .. n], y[0 ..n] );
		else {
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
	static void gemv( T )( char trans, int leny, int lenx, T alpha, const(T)* a, int lda, T* x, int incx, T beta, T *y, int incy ) {
		debug( blasCalls ) reportNaive_();
		
		if( !lenx || !alpha ) {
			if( trans == 'n' )
				scal( leny, beta, y, incy );
			else
				scal( lenx, beta, y, incy );
			return;
		}
		
		assert( checkVector( y, incy ) );
		assert( checkVector( x, incx ) );
		assert( checkMatrix( trans, leny, lenx, a, lda ) );
		
		if( trans != 'n' )
			std.algorithm.swap( lenx, leny );
		
		leny *= lda;
		auto ae = a + leny;
		T aux;
		if( beta == Zero!T ) {
			if( trans != 'c' ) {
				while( a != ae ) {
					aux = dot( lenx, ae, 1, x, incx ); aux *= alpha;
					*y = aux;
					ae += lda; y  += incy;
				}
			} else {
				while( a != ae ) {
					aux = dotc( lenx, ae, 1, x, incx ); aux *= alpha;
					*y  = aux;
					ae += lda; y  += incy;
				}
			}
		} else if( beta == One!T ) {
			if( trans != 'c' ) {
				while( a != ae ) {
					aux = dot( lenx, ae, 1, x, incx ); aux *= alpha;
					*y += aux;
					ae += lda; y  += incy;
				}
			} else {
				while( a != ae ) {
					aux = dotc( lenx, ae, 1, x, incx ); aux *= alpha;
					*y += aux;
					ae += lda; y  += incy;
				}
			}
		} else {
			T aux2;
			if( trans != 'c' ) {
				while( a != ae ) {
					aux = dot( lenx, ae, 1, x, incx ); aux *= alpha;
					aux2 = *y; aux2 *= beta; aux += aux2;
					*y = aux;
					ae += lda; y  += incy;
				}
			} else {
				while( a != ae ) {
					aux = dot( lenx, ae, 1, x, incx ); aux *= alpha;
					aux2 = *y; aux2 *= beta; aux += aux2;
					*y  = aux;
					ae += lda; y  += incy;
				}
			}
		}
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
	static void xgecopy( bool forceConjugate=false, T )( char transA, int m, int n, const(T)* a, int lda, T* b, int ldb ) if( isComplex!T || !forceConjugate ) {
		//debug( blasCalls ) reportNaive_();
		
		if( !m || !n )
			return;
		
		assert( checkMatrix( transA, m, n, a, lda ) );
		assert( checkMatrix( 'n', m, n, b, ldb ) );
		
		if( transA == 'n' ) {
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
		} else if( transA == 't' ) {
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
				// assume transA == 'c'
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
	static void xgeaxpy( bool forceConjugate=false, T )( char transA, int m, int n, T alpha, const(T)* a, int lda, T* b, int ldb ) if( isComplex!T || !forceConjugate ) {
		//debug( blasCalls ) reportNaive_();
		
		if( !m || !n )
			return;
		
		assert( checkMatrix( transA, m, n, a, lda ) );
		assert( checkMatrix( 'n', m, n, b, ldb ) );
		
		if( transA == 'n' ) {
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
		} else if( transA == 't' ) {
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
				// assume transA == 'c'
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
