module main;

import std.stdio;
import scid.matrix;
import matclosure;
import scalclosure;
import vecclosure;
import scid.linalg;
import scid.bindings.blas.dblas;
import std.traits;
import scid.common.traits;
import std.complex;

import std.conv;

unittest {
	auto A = matrix!double(2,2);
	
	static assert( is( typeof( -A ) : MatrixClosure!( MatrixUnary!( MatrixView!double, ScalarClosure!double, false, false ) ) ) );
	static assert( is( typeof( transpose(A*2)*3 ) : MatrixClosure!( MatrixUnary!( MatrixView!double, ScalarClosure!double, true, false ) ) ) );
}

int main(string[] argv) {
	auto M = matrix!double(2,2);
	auto N = M;
	
	M[0,0] = 1;
	writeln(to!string(N));

	readln();
	return 0;
}
