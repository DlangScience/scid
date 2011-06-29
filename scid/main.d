module main;

import scid.vector;
import scid.matrix;
import scid.symmatrix;
import scid.trimatrix;

import std.conv;
import std.stdio;
import std.complex;

int main() {
	// Vector/Matrix playground.
	//    Vector( T, VectorType vtype = VectorType.Column, alias ArrayTemplate = CowArray )
	//    Matrix( T, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	//    SymmetricMatrix( T, MatrixTriangle triangle_ = MatrixTriangle.Upper, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	//    TriangularMatrix( T, MatrixTriangle triangle_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	
	// Vector assignment & Matrix slicing
	auto v = Vector!double( 4 );
	auto A = Matrix!double( 3, 3 );
	
	writeln( "===== Matrix slicing =====" );
	A[] = 0.0;
	
	foreach( i  ; 0 .. 4 )
		v[ i ] = cast(double)( i + 1 );
	
	writeln( "v = ", to!string(v) );
	
	foreach( i ; 0 .. 3 ) {
		auto off = i % 2;
		A[ i ][] = v[ off .. off + 3 ];
	}
	
	writeln( "Checker board matrix:" );
	writeln( to!string( A ) );
	
	writeln( "Lower right corner matrix slice:" );
	writeln( to!string( A[1..3][1..3] ) );
	
	// Matrix storages
	writeln();
	writeln();
	writeln( "===== Matrix storage types =====" );
	{
		uint k = 0;
		foreach( i ; 0 .. 3 )
			foreach( j ; 0..3 )
				A[ i, j ] = ++k;
	}
	writeln( "Row major numbered matrix:" );
	writeln( to!string( A ) );
	
	writeln( "Column major storage: ", A.cdata() );
	{
		Matrix!(double, StorageOrder.RowMajor) B;
		B.assign( A ); // annoying bug with template opAssign() & postblit...
		writeln( "Row major storage: ", B.cdata() );
	}
	//
	writeln();
	{
		uint k = 0;
		foreach( i ; 0 .. 3 )
			foreach( j ; 0..3 )
				A[ j, i ] = ++k;
	}
	writeln( "Column major numbered matrix:" );
	writeln( to!string( A ) );
	
	writeln( "Column major storage: ", A.cdata() );
	{
		Matrix!(double, StorageOrder.RowMajor) B;
		B.assign( A ); // annoying bug with template opAssign() & postblit...
		writeln( "Row major storage: ", B.cdata() );
	}
	
	// Hermitian matrix;
	writeln();
	writeln();
	writeln( "===== Hermitian Matrix =====" );
	auto B = SymmetricMatrix!cdouble( 2 );	
	B[0, 0] = 1.0 + 1.0i;
	B[1, 1] = 1.0 + 1.0i;
	B[0, 1] = 2.0 + 2.0i;
	
	writeln( "Hermitian matrix: " );
	writeln( to!string(B) );
	writeln( "Storage: ", B.cdata() );
	
	{
		Matrix!cdouble M;
		M.assign( B );
		writeln( "General matrix storage: ", M.cdata() );
	}
	
	
	return 0;
}
	