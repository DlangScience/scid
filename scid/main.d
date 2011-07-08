module main;

import scid.vector;
import scid.matrix;
import std.stdio, std.conv, std.complex;

int main() {
	
	auto gen = Matrix!double(4, 4, null);
	auto tri = TriangularMatrix!double( [1.0,2.0,3.0,4.0,5.0,6.0] );
	auto sym = SymmetricMatrix!double( [1.0,2.0,3.0,4.0,5.0,6.0] );
	auto her = SymmetricMatrix!cdouble( [1.0+6.0i, 2.0+5.0i, 3.0+4.0i, 4.0+3.0i, 5.0+2.0i, 6.0+1.0i] );
	
	{
		double k = 0.0;
		foreach( i ; 0 .. gen.rows )
			foreach( j ; 0 .. gen.columns )
				gen[i, j] = ++ k;
	}
	
	writeln( "General matrix: " );
	writeln( gen.pretty );
	writeln( "Column-major data: ", gen.cdata );
	
	writeln();
	writeln( "Triangular matrix: " );
	writeln( tri.pretty );
	writeln( "Column Major data: ", tri.cdata );
	
	writeln();
	writeln( "Symmetric matrix: " );
	writeln( sym.pretty );
	writeln( "Column-major data: ", sym.cdata );
	
	writeln();
	writeln( "Hermitian matrix: " );
	writeln( her.pretty );
	writeln( "Column-major data: ", her.cdata );
	
	return 0;
}
	
