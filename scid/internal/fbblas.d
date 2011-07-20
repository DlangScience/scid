module scid.internal.fbblas;

import scid.vector;
import scid.common.traits, scid.common.meta;
import scid.internal.hlblas;
import scid.internal.expression;
import scid.matrix, scid.vector;

debug = blasCalls;

debug( blasCalls ) {
	import std.stdio, std.conv;
}

void fbmm( Transpose transA, Transpose transB, A, B, E, Dest )( E alpha, auto ref A a, auto ref B b, E beta, auto ref Dest dest )
if( is( BaseElementType!A : E ) && is( BaseElementType!B : E ) && is( BaseElementType!Dest : E ) ) {
	
	size_t m = transA ? a.rows    : a.columns;
	size_t n = transB ? a.columns : a.rows;
	
	if( !alpha ) {
		if( !beta ) dest.resizeOrClear( m, n );
		else        hlScal( beta, dest );
		return;
	} 
	
	if( !beta )
		dest.resizeOrClear( m, n, null );
	else
		assert( dest.rows == m && dest.columns == n, "Matrix size mismatch in matrix addition." );
		
		
	void doElement( size_t i, size_t j ) {
		E dot;
			
		static if( !transA && !transB )     dot = hlDot!( transA, transB )( a.row( i ),    b.column( j ) );
		else static if( !transA && transB ) dot = hlDot!( transA, transB )( a.row( i ),    b.row( j )    );
		else static if( transA && !transB ) dot = hlDot!( transA, transB )( a.column( i ), b.column( j ) );
		else static if( transA && transB )  dot = hlDot!( transA, transB )( a.column( i ), b.row( j )    );
			
		dot *= alpha; dot += beta * dest[i, j];
		dest[ i, j ] = dot;	
	}
		
	static if( Dest.storageOrder == StorageOrder.RowMajor ) {
		foreach( i ; 0 .. m ) {
			foreach( j ; 0 .. n ) doElement( i, j );
		}
	} else {
		foreach( j ; 0 .. n ) {
			foreach( i ; 0 .. m ) doElement( i, j );
		}
	}
}

auto fbdot( Transpose transA, Transpose transB, A, B )( auto ref A a, auto ref B b ) if( is(BaseElementType!A : BaseElementType!B ) ){
	alias BaseElementType!A T;
	
	assert( a.length == b.length, "Vector length mismatch in fallback dot product." );
	debug( blasCalls )  {
		string aStr = to!string(a) ~ (transA ? ".t" : "");
		string bStr = to!string(b) ~ (transB ? ".t" : "");
		write( "fb_dot( ", aStr , ", ", bStr, " ) => " );
	}
	int n = a.length;
	auto r = Zero!T;
		
	foreach( i ; 0 .. n ) {
		static if( !isComplex!T || (transA == transB) )   r += a[ i ]          * b[ i ];
		else static if(  transA && !transB )              r += gconj( a[ i ] ) * b[ i ];
		else static if( !transA &&  transB )              r += a[ i ]          * gconj( b[ i ] );
		else static assert( false, "I've proven p && !p." );
	}
		
	static if( transA && transB )
		r = r.conj;
		
	debug( blasCalls ) writeln( r );
	return r;
}

void fbconj( Transpose trans )( auto ref V vec ) {
	static if( trans && BaseElementType!V ) {
		auto n = vec.length;
		foreach( i ; 0 .. n )
			vec[ i ] = gconj( vec[ i ] );
	}
}

void fbmv( Transpose transM, Transpose transV, M, V, Dest, E )( E a, auto ref M m, auto ref V v, E b, auto ref Dest dest )
if( is( BaseElementType!M : E ) && is( BaseElementType!V : E ) && is( BaseElementType!Dest : E ) &&
(closureOf!V == Closure.RowVector || closureOf!V == Closure.ColumnVector) && closureOf!M == Closure.Matrix ) {
	
	auto rows = transM ? m.columns : m.rows;
	auto cols = transM ? m.rows    : m.columns;
	
	assert( cols == v.length, "fbmv: Dimensions mismatch in matrix-vector product - " ~ to!string( cols ) ~ " vs. " ~ to!string( v.length ) );
	
	if( b )
		assert( rows == dest.length, "fbmv: Dimensions mismatch in vector-vector addition - " ~ to!string( rows ) ~ " vs. " ~ to!string( dest.length ) );
	else
		dest.resizeOrClear( rows );
	
	assert( rows == dest.length );
	
	
	foreach( i ; 0 .. rows ) {
		static if( transM )
			dest[ i ] = hlDot!(transM,transV)( m.column( i ), v ) + b * dest[ i ];
		else
			dest[ i ] = hlDot!(transM,transV)( m.row( i ), v ) + b * dest[ i ];
	}
}




