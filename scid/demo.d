module scid.demo;
import std.range, std.exception;

import std.typetuple, std.demangle;
import std.complex;
import scid.common.traits, scid.common.meta;
import scid.internal.regionallocator;

version = demo;

version( demo ) {
	import scid.matvec;
	import std.stdio, std.conv;	
	
	void basicExpressions() {
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
	
	void rangeInterface() {
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
	
	void dataInterface() {
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
			TriangularMatrix!T,
			//TriangularMatrix!(T, MatrixTriangle.Lower ),
			//TriangularMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			//TriangularMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
			SymmetricMatrix!T,
			// SymmetricMatrix!(T, MatrixTriangle.Lower ),
			// SymmetricMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			// SymmetricMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor)
		) MatrixTypes;
	}
	
	/** Syntactically test all the operations on all the matrix types. */
	void opTest()() {
		alias TypeTuple!(double) ElementTypes;
		foreach(Element; ElementTypes) {
			enum z = Zero!Element;
			Element[][] minit = [[z, z, z], [z, z, z], [z, z, z]];
			foreach(LhsType; MatrixTypes!Element) {
				auto lhs = LhsType( minit );
				foreach(RhsType; MatrixTypes!Element) {
					auto rhs = RhsType( minit );
					
					eval(lhs + rhs*2.0);
					eval(lhs*2.0 - rhs);
					eval(lhs * rhs*2.0);
					eval(lhs.column(0)*2.0 + rhs.column(1));
					eval(lhs.row(0) - rhs.row(1)*2.0);
					eval(lhs.row(0) * rhs.column(0)*2.0);
					
					writeln( "<Temp (3 allocs)>" );
					// require temporaries
					eval( lhs.t*(lhs + rhs*2.0) );
					eval( (lhs.column(0) - rhs.column(0)*2.0).t*lhs );
					eval( (lhs.column(0)*2.0 + rhs.column(0)).t*lhs.column(1) );
					writeln( "</Temp>" );
				}
			}
		}
	}
	
	void main() {
	
		opTest();
		// basicExpressions();
		// rangeInterface();
		// dataInterface();
		//auto v = Vector!double([1.,2.,3.]);
		//writeln( eval((v-v*2.0).t * (v+v)));
		
		readln();
	}
}