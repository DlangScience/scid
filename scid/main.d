module main;

import scid.vector;
import scid.matrix;
import scid.symmatrix;
import scid.trimatrix;
import scid.vectorview;

import std.conv;
import std.stdio;
import std.complex;

void putMat( M )( string title, M m ) {
	writeln( title,':' );
	writeln( to!string(m) );
	writeln();
}

int main() {
	// Vector/Matrix playground.
	//    Vector( T, VectorType vtype = VectorType.Column, alias ArrayTemplate = CowArray )
	//    Matrix( T, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	//    SymmetricMatrix( T, MatrixTriangle triangle_ = MatrixTriangle.Upper, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	//    TriangularMatrix( T, MatrixTriangle triangle_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )
	
	// uninitialized matrix
	auto A = Matrix!double(5,7, null);
	
	putMat( "Uninitialized matrix", A );
	
	// initialized matrices
	{
		auto B = Matrix!double(2,2);
		auto C = Matrix!double(2,2, 42.0);
		
		putMat( "Initialized matrix (default)", B );
		putMat( "Initialized matrix (with 42)", C );
	}
		
	uint k = 0;
	foreach( i ; 0..5 )
		foreach( j ; 0..7 )
			A[ i, j ] = ++k;
	
	putMat( "Numbered matrix", A );
	putMat( "Transposed", A.t );
	
	writeln( "Foreach - by column: " );
	foreach( column ; A ) {
		foreach( value ; column )
			write( value, ' ' );
		writeln();
	}
	writeln();
	writeln( "Transposed foreach - by row: " );
	foreach( column ; A.t ) {
		foreach( value ; column )
			write( value, ' ' );
		writeln();
	}
	
	readln();
	return 0;
}
	