module scid.main;

import scid.vector;
import scid.matrix;
import scid.common.traits;
import scid.internal.expression, scid.internal.hlblas;
import std.stdio, std.conv, std.complex, std.traits, std.range;
import std.math, scid.util;


int main() {
	alias Vector!double                  dVec;
	alias Matrix!double                  dMat;
	alias Vector!cdouble                 cdVec;
	alias Matrix!cdouble                 cdMat;

	dVec v = [1.0,2.0,3.0];
	dMat m = dMat( 3, [ 1.0, 4.0, 7.0,
		                2.0, 5.0, 8.0,
		                6.0, 6.0, 9.0] );
	
	auto x = eval( v.t * v );
	static assert( is( typeof(x) : double ) );
	
	assert( x == 14.0 );
	
	dVec u = m * v;
	auto w = v.t * m;
	
	//assert( u == [14.0,32.0,50.0] );
	//assert( w == [30.0,36.0,42.0] );
	
	//x = eval(w * u);
	//assert( x == 3672.0 );
		
	
	readln();
	return 0;
}
	
