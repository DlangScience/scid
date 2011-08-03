module scid.lapack;
import scid.blas;
import scid.common.traits, scid.common.meta;
import std.algorithm, std.math;
import std.ascii, std.exception;

//debug = lapackCalls;
//version = nodeps;

debug( lapackCalls ) {
	import std.stdio;
	import scid.internal.assertmessages;
}



version( nodeps ) {
	private enum forceNaive = true;
} else {
	private enum forceNaive = false;
	static import scid.bindings.lapack.dlapack;
	alias scid.bindings.lapack.dlapack lapack_;
	
}

struct lapack {
	static void laswp( T )( int n, T *a, int lda, int k1, int k2, int *ipiv, int incx ) {
		debug( lapackCalls )
			writef( "laswp( %s, %s, %s, %s ) ", matrixToString( 'n', n, n, a, lda ), k1, k2, stridedToString( ipiv, n, incx ) );
		
		naive_.laswp( n, a, lda, k1, k2, ipiv, incx );
		
		debug( lapackCalls ) {
			writeln( "=> ", matrixToString( 'n', n, n, a, lda ) );
		}
																										 
	}
	
	static void getrs( char trans, T )( int n, int nrhs, T *a, int lda, int *ipiv, T *b, int ldb, ref int info ) {
		debug( lapackCalls )
			writef( "getrs( %s, %s, %s ) ", matrixToString( trans, n, n, a, lda ), ipiv[ 0 .. n ], matrixToString( trans, n, nrhs, b, ldb ) );
		
		static if( isFortranType!T && !forceNaive ) 
			lapack_.getrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info );
		else
			naive_.xgetrs!(trans, 'L')( n, nrhs, a, lda, ipiv, b, ldb, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}
	
	static void gesv( T )( int n, int nrhs, T *a, int lda, int* ipiv, T *b, int ldb, ref int info ) {
		debug( lapackCalls )
			writef( "gesv( %s, %s, %s ) ", matrixToString( trans, n, n, a, lda ), matrixToString( trans, n, nrhs, b, ldb ) );
		
		static if( isFortranType!T && !forceNaive ) 
			lapack_.gesv( n, nrhs, a, lda, ipiv, b, ldb, info );
		else
			naive_.gesv( n, nrhs, a, lda, ipiv, b, ldb, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}
	
	static void getrf( T )( int m, int n, T* a, int lda, int* ipiv, ref int info ) {
		debug( lapackCalls )
			writef( "getrf( %s ) ", matrixToString( 'N', m, n, a, lda ) );
		
		static if( isFortranType!T && !forceNaive )
			lapack_.getrf( m, n, a, lda, ipiv, info );
		else
			naive_.getrf( m, n, a, lda, ipiv, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', m, n, a, lda ) );
	}
	
	static void trtri( char uplo, char diag, T )( int n, T *a, int lda, ref int info ) {
		debug( lapackCalls )
			writef( "trtri( %s, %s ) ", uplo, matrixToString( 'N', n, n, a, lda ) );
		
		static if( isFortranType!T && !forceNaive )
			lapack_.trtri( uplo, diag, n, a, lda, info );
		else
			naive_.trtri!( uplo, diag )( n, a, lda, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', n, n, a, lda ) );
	}
	
	static void getri( T )( int n, T* a, int lda, int* ipiv, T* work, int lwork, ref int info ) {
		debug( lapackCalls )
			writef( "getri( %s ) ", matrixToString( 'N', n, n, a, lda ) );
		
		static if( isFortranType!T && !forceNaive )
			lapack_.getri( n, a, lda, ipiv, work, lwork, info );
		else
			naive_.getri( n, a, lda, ipiv, work, lwork, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', n, n, a, lda ) );
	}
	
	// Extended LAPACK
	static void xgetrs( char trans, char side, T )( int n, int nrhs, T *a, int lda, int *ipiv, T *b, int ldb, ref int info ) {
		debug( lapackCalls )
			writef( "xgetrs( %s, %s, %s, %s ) ", side, matrixToString( trans, n, n, a, lda ), ipiv[ 0 .. n ], matrixToString( trans, n, nrhs, b, ldb ) );
		
		naive_.xgetrs!( trans, side )( n, nrhs, a, lda, ipiv, b, ldb, info );
		
		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}
}

private struct naive_ {
	private static void reportNaive_() {
		debug( lapackCalls )
			write( "<n> " );
	}
	
	private static void reportNaiveln_() {
		debug( lapackCalls )
			writeln( "<n> ..." );
	}
	
	static void laswp( T )( int n, T *a, int lda, int k1, int k2, int *ipiv, int incx ) {
		reportNaiveln_();
		
		// convert FORTRAN indices
		--k1; --k2;
		
		enforce( n >= 0 );
		enforce( incx != 0 );
		enforce( a );
		enforce( ipiv );
		
		if( n == 0 )
		    return;
		
		if( incx > 0 ) {
		    for( int i = k1; i <= k2 ; ++ i ) {
		        int pivot = ipiv[ i ] - 1; // convert FORTRAN index
		        if( pivot != i )
		            blas.swap( n, a + i, lda, a + pivot, lda );
		    }
		} else {
		    for( int i = k2; i >= k1 ; -- i ) {
		        int pivot = ipiv[ i ] - 1; // convert FORTRAN index
		        if( pivot != i )
		            blas.swap( n, a + i, lda, a + pivot, lda );
		    }	
		}
	}
	
	static void xgetrs( char trans_, char side_, T )( int n, int nrhs, T *a, int lda, int *ipiv, T *b, int ldb, ref int info ) {
		enum trans = cast(char) toUpper( trans_ );
		enum side = cast(char) toUpper( side_ );
		
		reportNaiveln_();
		
		enforce( n >= 0 );
		enforce( nrhs >= 0 );
		enforce( lda >= max( 1, n ) );
		enforce( ldb >= max( 1, n ) );
		
		if( n == 0 || nrhs == 0 )
			return;
		
		static if( trans == 'N' ) {
			static if( side == 'R' ) {
				blas.trsm!( 'R', 'U', 'N', 'N' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( 'R', 'L', 'N', 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				for( int i = n - 1; i >= 0 ; -- i ) {
					int pivot = ipiv[ i ] - 1;
					if( pivot != i )
						blas.swap( n, b + i * ldb, 1, b + pivot * ldb, 1 );
				}
			} else if( side == 'L' ) {
				lapack.laswp( nrhs, b, ldb, 1, n, ipiv, 1 );
				blas.trsm!( 'L', 'L', 'N', 'U' )( n, nrhs, One!T, a, lda, b, ldb );	
				blas.trsm!( 'L', 'U', 'N', 'N' )( n, nrhs, One!T, a, lda, b, ldb );
				
			}
		} else {
			static if( side == 'R' ) {
				for( int i = 0; i < n ; ++ i ) {
					int pivot = ipiv[ i ] - 1;
					if( pivot != i )
						blas.swap( n, b + i * ldb, 1, b + pivot * ldb, 1 );
				}
				blas.trsm!( 'R', 'L', trans, 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( 'R', 'U', trans, 'N' )( n, nrhs, One!T, a, lda, b, ldb );
			} else static if( side == 'L' ) {
				blas.trsm!( side, 'U', trans, 'N' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( side, 'L', trans, 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				lapack.laswp( nrhs, b, ldb, 1, n, ipiv, -1 );
			}
		}
	}
	
	static void gesv( T )( int n, int nrhs, T *a, int lda, int* ipiv, T *b, int ldb, ref int info ) {
		reportNaiveln_();
		
		enforce( n >= 0 );
		enforce( nrhs >= 0 );
		enforce( lda >= max( 1, n ) );
		enforce( ldb >= max( 1, n ) );
		
		lapack.getrf( n,n, a, lda, ipiv, info );
		if( info == 0 )
			lapack.getrs!'N'( n, nrhs, a, lda, ipiv, b, ldb, info );
	}
	
	static void getri( T )( int n, T* a, int lda, int* ipiv, T* work, int lwork, ref int info ) {		
		reportNaiveln_();
		info = 0;
		
		work[ 0 ] = n * 2;
		if( lwork == -1 ) {
			info = 0;
			return;
		}
		
		if( n == 0 )
			return;
		
		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}
	
		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}
		
		lapack.trtri!( 'U', 'N' )( n, a, lda, info );
		if( info > 0 )
			return ;
		
		for( int j = n - 1; j >= 0 ; -- j ) {
			foreach( i ; j + 1 .. n ) {
				work[ i ] = get( i, j );
				set( Zero!T, i, j );
			}
			
			
			if( j < n-1 ) {
				blas.gemv!'N'( n, n - j - 1, MinusOne!T, a + (j+1)*lda,
					lda, work + j + 1, 1, One!T, a + j * lda, 1 );
			}
		}
		
		for( int j = n - 1 ; j >= 0 ; -- j ) {
			int pivot = ipiv[ j ] - 1; // convert from FORTRAN index
			if( pivot != j )
				blas.swap( n, a + j * lda, 1, a + pivot * lda, 1 );
		}
		
	}
	
	static void trtri( char uplo_, char diag_, T )( int n, T *a, int lda, ref int info ) {
		reportNaiveln_();
		
		enum char uplo = toUpper(uplo_);
		enum char diag = toUpper(diag_);
		
		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}
	
		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}
		
		T ajj;
		if( uplo == 'U' ) {
			foreach( j ; 0 .. n ) {
				if( diag == 'N' ) {
					// assert( get( j, j ) != Zero!T, "fbti: Singular matrix in inverse." );
					if( get( j, j ) == Zero!T ) {
						info = j;
						return;
					}
					
					set( One!T / get(j,j), j, j );
					ajj = -get( j, j );
				} else {
					ajj = MinusOne!T;
				}
				
				blas.trmv!( 'U', 'N', diag )( j, a, lda, a + j*lda, 1 );
				blas.scal( j, ajj, a + j*lda, 1 );	
			}
		} else {
			for( int j = n - 1 ; j >= 0 ; -- j ) {
				if( diag == 'N' ) {
					// assert( get( j, j ) != Zero!T, "fbti: Singular matrix in inverse." );
					set( One!T / get(j,j), j, j );
					ajj = -get( j, j );
				} else {
					ajj = MinusOne!T;
				}
				
				if( j < n-1 ) {
					blas.trmv!( 'L', 'N', diag )( n - 1 - j, a + (j+1)*lda + j+1, lda, a + j*lda + j+1, 1 );
					blas.scal( n - j, ajj, a + (j+1) * lda + j, 1 );
				}
			}
		}
	}
	
	static void getrf( T )( int m, int n, T* a, int lda, int* pivot, ref int info ) {
		reportNaiveln_();
		
		n = min( m, n );
		
		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}
	
		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}
		
		foreach( k ; 0 .. n ) {
			pivot[ k ] = k;
			T maxSoFar = abs( get( k, k ) );
			foreach( j ; k + 1 .. n ) {
				T cur = abs( get(j, k) );
				if( maxSoFar <= cur ) {
					maxSoFar = cur;
					pivot[ k ] = j;
				}
			}
		
			if( pivot[ k ] != k ) {
				foreach( j ; 0 .. n ) {
					T aux = get(k, j);
					set( get(pivot[k], j), k, j );
					set( aux, pivot[k], j );
				}
			}
		
			if( get(k,k) == Zero!T )
				info = k;
		
			foreach( i ; k + 1 .. n )
				set!"/"( get(k,k), i, k );
		
			foreach( i ; k + 1 .. n ) {
				foreach( j ; k + 1 .. n ) {
					set!"-"( get(i,k) * get(k,j), i, j );
				}
			}
			
			++ pivot[ k ]; // convert to FORTRAN index
		}
	}
}