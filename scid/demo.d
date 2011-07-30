module scid.demo;

import std.range, std.exception;

import std.typetuple;
import std.complex;
import scid.common.traits, scid.common.meta;
import scid.internal.regionallocator, scid.storages;

//version = demo;

version( demo ) {
	import scid.matvec;
	import std.stdio, std.conv;	
	
	void basicExpressions()() {
		writeln();
		writeln( "======================= Basic Expressions =======================" );
		// The simplest constructors take a built-in array.
		auto mat = Matrix!double( [ [1., 2., 3.], [2., 3., 4.] ] );
		auto vec = Vector!double( [1.,2.,3.,4.,5.,6.] );
		
		// Printing matrices and vectors
		writeln( "vec = ", vec.toString() ); // Using toString
		writeln( "mat = " );
		writeln( mat.pretty );               // pretty prints rows on multiple lines
	
		writeln();
		
		// Slicing vectors and matrices
		auto vecSlice = vec[ 1 .. 3 ];
		auto matSlice = mat[ 0 .. 2 ][ 1 .. 3 ];
	
		// RowVector times Matrix. The t property transposes vectors and matrices.
		auto w = eval( vecSlice.t * mat ); 
		writeln( "Expr 1: ", vecSlice.toString(), " * ", mat.toString(), " = ", w.toString() );
		
		enforce( w == [8.0, 13.0, 18.0] ); 
	
		// More complicated expression. mat[0][0..2] gets a view of a slice of the first row. The .t is neccessary since
		// vec is a column vector while the matrix slice is a row vector.
		vecSlice[] = vec[ 0 .. 2 ] * 5.0 - mat[0][0..2].t;
		writeln( "Expr 2: ", vec[0 .. 2].toString(), " * 5 - ", mat[0][0..2].toString(), " = ", vecSlice.toString() );
		
		enforce( vecSlice == [4.0, 8.0] );
	
		// One can use array literals as vector literals most of the time:
		double x = eval( [2.0, -1.0].t * vecSlice );
		writeln( "Expr 3: [2.0, -1.0].t * ", vecSlice.toString, " = ", x );
		
		enforce( x == 0.0 );
	}
	
	void rangeInterface()() {
		writeln();
		writeln( "======================== Range Interface ========================" );
		
		// InputRange. Using foreach with vectors
		auto v = Vector!double([1.0, 2.0, 3.0]);
		auto sum = 0.;
		foreach( e ; v ) sum += e;
		writeln( "The sum of ", v.toString(), "'s elements is ", sum );
		
		// Matrices are iterated by major subvector (e.g. columns for column-major matrices). Importantly, the elements
		// in the iterations are views. Therefore changing an element affects the matrix, no matter if ref is used or not.
		auto rowMat = Matrix!(double, StorageOrder.RowMajor)([ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ]);
		auto colMat = Matrix!double( [ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ] );
		
		uint i = 0;
		writeln( "Row major matrix: " );
		foreach( r ; rowMat ) {
			writeln( "Row ", i, ": ", r.toString() );
			
			if( i == 0 )      enforce( r == [1.,2.,3.] );
			else if( i == 1 ) enforce( r == [4.,5.,6.] );
			else if( i == 2 ) enforce( r == [7.,8.,9.] );
			
			i ++;
		}
		writeln();
		
		i = 0;
		writeln( "Column major matrix: " );
		foreach( c ; colMat ) {
			writeln( "Column ", i, ": ", c.toString() );
			
			if( i == 0 )      enforce( c == [1.,4.,7.] );
			else if( i == 1 ) enforce( c == [2.,5.,8.] );
			else if( i == 2 ) enforce( c == [3.,6.,9.] );
			
			i ++;
		}
		writeln();
	}
	
	void dataInterface()() {
		writeln();
		writeln( "========================= Data Interface =========================" );
		
		// Most storages provide data & cdata methods which allow access to the raw memory that can be passed to BLAS
		// or used by custom function.
		// cdata() - returns a const pointer to a memory block which, due to copy-on-write, might be shared.
		// data()  - returns a mutable pointer to the memory block. It first ensures that the memory is not shared.
		
		// By printing the memory of the two following matrices we can see the difference between the storage
		// orders:
		auto rowMat = Matrix!(double, StorageOrder.RowMajor)([ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ]);
		auto colMat = Matrix!double( [ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ] );
		
		writeln( "Row major data   : ", rowMat.cdata[ 0 .. 9 ] );
		writeln( "Column major data: ", colMat.cdata[ 0 .. 9 ] );
		
		// Assigning colMat to a new matrix will cause the two matrices to share the data
		auto otherMat = colMat;
		enforce( otherMat.cdata == colMat.cdata );
		
		// Calling data will cause the memory to be copied though:
		enforce( otherMat.data  != colMat.data );
		enforce( otherMat.cdata != colMat.cdata );
	}
	
	/** Some of these types have to be disabled otherwise the compiler runs out of memory. */
	template MatrixTypes( T ) {
		alias TypeTuple!(
			Matrix!T,
			//Matrix!(T,StorageOrder.RowMajor),
			//TriangularMatrix!T,
			//TriangularMatrix!(T, MatrixTriangle.Lower ),
			//TriangularMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			//TriangularMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
			//SymmetricMatrix!T,
			// SymmetricMatrix!(T, MatrixTriangle.Lower ),
			// SymmetricMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			// SymmetricMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor)
		) MatrixTypes;
	}
	
	/** Syntactically test all the operations on all the matrix types. */
	void opTest()() {
		alias TypeTuple!(double) ElementTypes;
		foreach(T; ElementTypes) {
			enum z = One!T;
			T[][] minit = [[z, z, z], [z, z, z], [z, z, z]];
			foreach(LhsType; MatrixTypes!T) {
				auto lhs = LhsType( minit );
				foreach(RhsType; MatrixTypes!T) {
					auto rhs = RhsType( minit );
					
					eval(lhs + rhs*z);
					eval(lhs*z - rhs);
					eval(lhs * rhs*z);
					eval(lhs.column(0)*z + rhs.column(1));
					eval(lhs.row(0) - rhs.row(1)*z);
					eval(lhs.row(0) * rhs.column(0)*z);
					
					eval( lhs.t*(lhs + rhs*z) );
					eval( (lhs.column(0) - rhs.column(0)*z).t*lhs );
					eval( (lhs.column(0)*z + rhs.column(0)).t*lhs.column(1) );
					
					lhs[] += rhs;
					lhs[] -= rhs;
					lhs[] *= rhs;
					lhs[] *= z;
					lhs[] /= z;
					
					lhs[0][] *= z;
					lhs[0][] /= z;
					lhs[0][] += lhs[][0].t;
					lhs[][0] -= lhs[1][].t;
					
					lhs[1..3][1..3] *= rhs[0..2][0..2];
					
				}
			}
		}
	}
	import std.string;
	import scid.matvec;
	import std.math;
	
	
	void testMat( M, E )( auto ref M m, size_t r, size_t c, E[] expected ) {
		debug {
			enum epsilon = 1e-3;
			assert( m.rows == r, format("Wrong no. of rows %d vs %d", m.rows, r) );
			assert( m.columns == c, format("Wrong no. of rows %d vs %d", m.columns, c ) );
			auto a = m.cdata[ 0 .. expected.length ];
			auto b = to!(BaseElementType!M[])(expected.dup);
			b[] -= a[];
			foreach( i, x ; b ) {
				assert( abs(x) <= epsilon,
				   "Expected " ~ to!string(expected) ~ ", got "~ to!string(a) ~ " (" ~ to!string(i) ~ ", " ~ to!string(abs(x)) ~ ")"  );
			}
		}
	}
	
	void testVec( V, E )( auto ref V v, E[] expected ) {
		debug {
			assert( v.length == expected.length, format("Wrong vector length: %d vs. %d", v.length, expected.length) );
			assert( v.cdata[ 0 .. expected.length ] == to!(BaseElementType!V)(expected) );
		}
	}
	
	import scid.ops.expression;
	void dMatOpsTest()() {
		alias Matrix!double            dGeMat;
		alias SymmetricMatrix!double   dSyMat;
		
		auto a = dGeMat( 3, [1.,2,3,4,5,6,7,8,9] );
		auto b = dGeMat( 3, [1.,2,3,4,5,6] );
		
		dGeMat c = b * a;
		testMat( c, 2, 3, [22,28,49,64,76,100] );
		c[] = c[0 .. 2][ 0 .. 2 ].t * ( (b[][0] - a[1..3][0]).t * eval(c[][0]) ) / 50.;
		testMat( c, 2, 2, [-22,-49,-28,-64] );
		
		dSyMat s = c.t*c;
		
		testMat( s, 2, 2, [2885, 3752, 4880] );
		
		auto d = eval( s - dSyMat([2800.,3700,4800]) );
		static assert( is( typeof(d) : dSyMat ) );
		testMat( d, 2, 2, [85,52,80] );
		assert( d[1][0] == 52 );
		
		auto e = eval( d - b[0..2][1..3]*10 );
		static assert( is( typeof(e) : dGeMat ) );
		testMat( e, 2, 2, [ 55, 12, 2, 20 ] );
	}
	
	void zMatOpsTest()() {
		alias Matrix!cdouble            zGeMat;
		alias SymmetricMatrix!cdouble   zSyMat;
		
		auto a = zGeMat( 3, [1.+4.i,2+3.i,3+2.i,4+1.i,5+0.i,6-1.i,7-2.i,8-3.i,9-4.i] );
		auto b = zGeMat( 3, [1.+2.i,2.+1.i,3+0.i,4-1.i,5-2.i,6-3.i] );
				   
		zGeMat c = b * a;
		testMat( c, 2, 3, [ 18 + 19.i, 33 + 22i, 45 - 8i, 60 - 23i, 72 -35i, 87 - 68i ] );
		
		
		c[] = c[0 .. 2][ 0 .. 2 ].t * ( (b[][0] - a[1..3][0]).t * eval(c[][0]) ) / (10.+0.i);
		testMat( c, 2, 2, [-146.60 + 192.80i, -422.00 -  28.60i, -281.60 + 235.40i, -575.00 - 151.60i] );
		
		c[] = c + zGeMat([[150-190i,280-230i], [430+28i,570+150i]]);
		zSyMat s = c.t*c;
		
		testMat( s, 2, 2, [ 83.760 +  0.000i,  
		                   -29.360 +  7.040i,
		                    59.280 +  0.000i ]);
		assert( abs(s[1][0] + 29.360 + 7.040i) <= 1e-3 );
		
		auto d = eval( s - zSyMat([80.+0.i,-28,59]) );
		static assert( is( typeof(d) : zSyMat ) );
		
		testMat( d, 2, 2, [ 3.76 + 0.0i, -1.36 + 7.04i, 0.28 + 0.0i ] );
		assert( abs(d[1][0] + 1.36 + 7.04i) <= 1e-3 );
		
		
		auto e = eval( d + b[0..2][1..3]*(10.+0.i) );
		static assert( is( typeof(e) : zGeMat ) );
		
		testMat( e, 2, 2, [ 33.760 +  0.000i, 38.640 - 17.040i, 48.640 - 12.960i, 60.280 - 30.000i] );
	}
	import scid.ops.eval, scid.common.meta;
	void main() {
		auto dm = DiagonalMatrix!double([1,2,3,4,5]);
		//writeln( dm.column(0).pretty );
		dm[] =dm * dm;
		
		writeln( dm.pretty );
		//static assert( isScalar!(BaseElementType!double) );
		//writeln( MinusOne!double );
		//auto x = Matrix!double([[1.0, 4, 3], [4.0, 5, 6], [7.0, 8, 9]]);
		//auto vec = Vector!double([8.0, 2, 3]);
		//evalSolve(x, vec);
	    //writeln(vec.pretty);
		//dMatOpsTest();
		//zMatOpsTest();
		readln();
		//writeln( y );
		// alias Matrix!cdouble Mat;
		//alias Matrix!double Mat;
		//auto x = Mat([[1,2,3],[4,5,6],[7,8,9]]);
		//x[] += 2*x.t;
		//writeln(x.pretty);
		//writeln( x.pretty );
		//Mat x2 = (x.t+x);
		//writeln( x2.pretty );
		
		//opTest();
		// basicExpressions();
		// rangeInterface();
		// dataInterface();
		//auto v = Vector!double([1.,2.,3.]);
		//writeln( eval((v-v*2.0).t * (v+v)));
	}
}
