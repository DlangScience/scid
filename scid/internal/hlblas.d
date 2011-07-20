module scid.internal.hlblas;

import std.traits, std.range, std.array, std.conv;
import scid.vector, scid.matrix;
import scid.internal.expression;
import scid.internal.fbblas;
import scid.bindings.blas.dblas;
import scid.common.traits, scid.common.meta;
import std.math : conj;
import std.complex;

template isStridedStorage( T, E = BaseElementType!T ) {
	enum isStridedStorage =
		(is( typeof(T.stride) : size_t ) &&
		 is( typeof(T.cdata)  : const(E)* ) &&
		 is( typeof(T.data)   : E* ) &&
		 is( typeof(T.length) : size_t )) ||
		 is( T : E[] );
}

template isGeneralMatrix( T, E = BaseElementType!T ) {
	enum isGeneralMatrix =
		(is( typeof(T.leading) : size_t ) &&
		 is( typeof(T.cdata)   : const(E)* ) &&
		 is( typeof(T.data)    : E* ) &&
		 is( typeof(T.rows)    : size_t ) &&
		 is( typeof(T.columns) : size_t )) ||
		 is( T : E[] );
}

auto gconj( cdouble z        ) { return conj( z ); }
auto gconj( cfloat  z        ) { return conj( z ); }
auto gconj( Complex!double z ) { return z.conj; }
auto gconj( Complex!float  z ) { return z.conj; }

template transNot( Transpose tr )             { enum transNot = cast(Transpose) !tr; }
template transXor( Transpose a, Transpose b ) { enum transXor = cast(Transpose) (a ^ b); } 

mixin template StridedAxpyScalDot() {
	void axpy( S )( ElementType alpha, ref S rhs ) if( isStridedStorage!(S, ElementType) ) {
		hlStridedAxpy( alpha, rhs, this );
	}
	
	void scal()( ElementType alpha ) {
		hlStridedScal( alpha, this );	
	}
	
	auto dot( Transpose tr, Other )( ref Other other ) if( isStridedStorage!(Other, ElementType) ) {
		return hlStridedDot!tr( this, other );	
	}
}

mixin template GeneralMatrixAxpyScal() {
	void axpy( S )( ElementType alpha, ref S rhs ) if( isGeneralMatrix!(S, ElementType) ) {
		assert( rhs.rows == this.rows && rhs.columns == this.columns, "Matrix size mismatch in addition." );
		static assert( isGeneralMatrix!(typeof(this), ElementType) );
		hlGeneralMatrixAxpy!(S.storageOrder, storageOrder)(
			this.rows, this.columns,
			alpha,
			rhs.cdata, rhs.leading,
			this.data, this.leading
		);
	}
	
	void scal()( ElementType alpha ) {
		static assert( isGeneralMatrix!(typeof(this), ElementType) );
		hlGeneralMatrixScal!storageOrder( this.rows, this.columns, alpha, this.data, this.leading );
	}
}

template transposeClosure( Closure A, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( A == Closure.RowVector )         enum transposeClosure = Closure.ColumnVector;
		else static if( A == Closure.ColumnVector ) enum transposeClosure = Closure.RowVector;
		else                                        enum transposeClosure = A;
	} else
		enum transposeClosure = A;
}

template transposeVectorType( VectorType v, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( v == VectorType.Row )         enum transposeVectorType = VectorType.Column;
		else static if( v == VectorType.Column ) enum transposeVectorType = VectorType.Row;
		else                                     enum transposeVectorType = v;
	} else
		enum transposeVectorType = v;
}

template transposeStorageOrder( StorageOrder S, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( S == StorageOrder.RowMajor )         enum transposeStorageOrder = StorageOrder.ColumnMajor;
		else static if( S == StorageOrder.ColumnMajor ) enum transposeStorageOrder = StorageOrder.RowMajor;
		else                                            enum transposeStorageOrder = S;
	} else
		enum transposeStorageOrder = S;
}

private {
	@property auto ref storage( S )( ref S[] x ) if( isPrimitiveLiteral!S && closureOf!S != Closure.Matrix ) {
		return x;
	}

	@property auto data( S )( S[] x ) if( isPrimitiveLiteral!S && closureOf!S != Closure.Matrix ) {
		return x.ptr;
	}

	@property const(S)* cdata( S )( S[] x ) if( isPrimitiveLiteral!S && closureOf!S != Closure.Matrix ) {
		return x.ptr;
	}

	@property auto stride( S )( const(S[]) x ) if( isPrimitiveLiteral!S && closureOf!S != Closure.Matrix ) {
		return 1;
	}
}

void hlStridedCopy( Source, Dest )( auto ref Source source, ref Dest dest ) if( isStridedStorage!(Source,BaseElementType!Dest) && isStridedStorage!Dest ) {
	dest.resizeOrClear( source.length, null );
	copy( dest.length, source.cdata, source.stride, dest.data, dest.stride );
}

void hlStridedAxpy( Scalar, Source, Dest )( Scalar alpha, auto ref Source source, ref Dest dest ) if( isStridedStorage!(Source,Scalar) && is( Scalar : BaseElementType!Dest ) ) {
	assert( source.length == dest.length, "Length mismatch in strided addition." );
	axpy( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
}

void hlStridedScal( Scalar, Dest )( Scalar alpha, ref Dest dest ) if( isStridedStorage!Dest && is( Scalar : BaseElementType!Dest ) ) {
	scal( dest.length, alpha, dest.data, dest.stride );
}

auto hlStridedDot( Transpose transA, A, B )( auto ref A a, auto ref B b ) if( isStridedStorage!(A,BaseElementType!B) && isStridedStorage!B ) {
	assert( a.length == b.length, "Length mismatch in strided dot." );
	static if( !isComplex!(BaseElementType!A) )
		auto r = dot( a.length, a.cdata, a.stride, b.cdata, b.stride );	
	else static if( transA )
		auto r = dotc( a.length, a.cdata, a.stride, b.cdata, b.stride );
	else
		auto r = dotu( a.length, a.cdata, a.stride, b.cdata, b.stride );	
	return r;
}

auto hlGeneralMatrixScal( StorageOrder order, E )( size_t rows, size_t cols, E alpha,  E* data, size_t leading ) {
	static if( order == StorageOrder.RowMajor ) {
		alias rows major;
		alias cols minor;
	} else {
		alias rows minor;
		alias cols major;
	}	 
	
	if( leading == minor ) {
		// zero everything
		scal( rows * cols, alpha, data, 1 );
	} else {
		// zero by-major
		auto e = data + major * leading;
		for( ; data < e ; data += leading )
			scal( minor, alpha, data, 1 );
	}
}

auto hlGeneralMatrixCopy( StorageOrder srcOrder, StorageOrder dstOrder, E )
	( size_t rows, size_t cols, const(E)* source, size_t srcLeading, E* dest, size_t dstLeading ) {
	
	static if( srcOrder == StorageOrder.RowMajor ) { alias rows srcMajor; alias cols srcMinor; }
	else                                           { alias rows srcMinor; alias cols srcMajor; }
	
	static if( dstOrder == StorageOrder.RowMajor ) { alias rows dstMajor; alias cols dstMinor; }
	else                                           { alias rows dstMinor; alias cols dstMajor; }
	
	static if( srcOrder == dstOrder ) {
		// if they have they have the same storage order
		if( srcLeading == srcMinor && dstLeading == dstMinor ) {
			// if both matrices own all of the data, just vector-copy the data
			copy( rows * cols, source, 1, dest, 1 );
		} else {
			// if they don't all the data copy major-by-major
			auto e = dest + dstMajor * dstLeading;
			while( dest < e ) {
				copy( dstMinor, source, 1, dest, 1 );
				dest   += dstLeading;
				source += srcLeading;
			}
		}
	} else {
		// if they don't have the same order, copy by-minor from source to dest by-major
		auto e = dest + dstMajor * dstLeading;
		while( dest < e ) {
			copy( dstMinor, source, srcLeading, dest, 1 );
			dest += dstLeading;
			++ source;
		}
	}
}

auto hlGeneralMatrixAxpy( StorageOrder srcOrder, StorageOrder dstOrder, E )
	( size_t rows, size_t cols, E alpha, const(E)* source, size_t srcLeading, E* dest, size_t dstLeading ) {
	
	static if( srcOrder == StorageOrder.RowMajor ) { alias rows srcMajor; alias cols srcMinor; }
	else                                           { alias rows srcMinor; alias cols srcMajor; }
	
	static if( dstOrder == StorageOrder.RowMajor ) { alias rows dstMajor; alias cols dstMinor; }
	else                                           { alias rows dstMinor; alias cols dstMajor; }
	
	static if( srcOrder == dstOrder ) {
		// if they have they have the same storage order
		if( srcLeading == srcMinor && dstLeading == dstMinor ) {
			// if both matrices own all of the data, just vector-axp the data
			axpy( rows * cols, alpha, source, 1, dest, 1 );
		} else {
			// if they don't all the data axpy major-by-major
			auto e = dest + dstMajor * dstLeading;
			while( dest < e ) {
				axpy( dstMinor, alpha, source, 1, dest, 1 );
				dest   += dstLeading;
				source += srcLeading;
			}
		}
	} else {
		// if they don't have the same order, axpy by-minor from source to dest by-major
		auto e = dest + dstMajor * dstLeading;
		while( dest < e ) {
			axpy( dstMinor, alpha, source, srcLeading, dest, 1 );
			dest += dstLeading;
			++ source;
		}
	}
}


void hlScal( Scalar, Dest )( auto ref Scalar scalar, auto ref Dest dest ) {
	alias BaseElementType!Dest T;
	alias closureOf!Dest         closure;
	
	static assert( is( T : BaseElementType!Scalar ),
		"Types in " ~ to!string(closure) ~ "-Scalar multiplication do not match (" ~
		T.stringof ~ " and " ~ (BaseElementType!Scalar).stringof ~ ")."
	);
	
	static assert( closureOf!Scalar == Closure.Scalar,
		"Cannot perform " ~ to!string(closure) ~ "-" ~ to!string(closureOf!Scalar) ~
		"multiplication in place."
	);
	
	static if( __traits( compiles, dest.storage.scal( scalar ) ) ) {
		// custom expression scaling
		dest.storage.scal( scalar );
	} else static if( __traits( compiles, dest.storage.scal( eval(scalar) ) ) ) {
		// custom literal scaling
		dest.storage.scal( eval(scalar) );
	} else static if( closure == Closure.RowVector || closure == Closure.ColumnVector ) {
		// vector fallback scaling, anything that supports indexing works
		auto l = dest.length;
		auto scalarValue = eval(scalar);
		debug( blasCalls ) write( "fb_scal( ", scalarValue, ", ", to!string(dest), " ) => " );
		for( size_t i = 0; i < l ; ++ i )
			dest[ i ] *= scalarValue;
		debug( blasCalls ) writeln( to!string(dest) );
	} else static if( closure == Closure.Matrix ) {
		auto scalarValue = eval(scalar);
		static if( dest.storageOrder == StorageOrder.RowMajor ) {
			size_t n = dest.rows;
			for( auto i = 0; i < n ; ++ i )
				hlScal( alpha, dest.row( i ) );
		} else {
			size_t n = dest.columns;
			for( auto i = 0; i < n ; ++ i )
				hlScal( alpha, dest.column( i ) );
		}
	} else static assert( false );
		
}

void hlCopy( Transpose tr = Transpose.no, Source, Dest )( auto ref Source source, auto ref Dest dest ) {
	alias BaseElementType!Source                  T;
	alias transposeClosure!(closureOf!Source, tr) srcClosure;
	alias closureOf!Dest                          dstClosure;
	
	// check correct operation
	static assert( isLiteral!Dest, Dest.stringof ~ " is not a lvalue." );
	
	static assert( is( T : BaseElementType!Dest ),
		"Types in assignment do not have the same element type (" ~
		T.stringof ~ " and " ~ (BaseElementType!Dest).stringof ~ ")."
	);
	
	static assert( srcClosure == dstClosure,
		"Incompatible types for assignment '" ~ to!string(srcClosure) ~
		"' and '" ~ to!string(dstClosure) ~ "'."
	);

	static if( __traits( compiles, dest.storage.copy!tr( source.storage ) ) ) {
		// call custom copy routine on destination
		dest.storage.copy!tr( source.storage );
	} else static if( __traits( compiles, source.storage.copyRight!tr( dest.storage ) ) ) {
		// call custom copy routine on source
		source.storage.copyRight!tr( dest.storage );
	} else static if( isExpression!Source ) {
		enum op = Source.operation;
		static if( isScalingOperation!op ) {
			// copy from a scaling operation - clear the destination & axpy with 1*source
			static if( srcClosure == Closure.Matrix )
				// matrix resizing
				dest.storage.resizeOrClear( source.rows, source.columns );
			else
				// vector resizing
				dest.storage.resizeOrClear( source.length );
				
			hlAxpy!tr( source.rhs, source.lhs, dest );
		} else static if( isAddition!op ) {
			// copy from an addition - copy the first operand and axpy with the second
			hlCopy!tr( source.lhs, dest );
			hlAxpy!tr( One!T, source.rhs, dest );
		} else static if( isTransposition!op ) {
			// copy from a transposition operation - recurse with the negated tranposition flag
			hlCopy!( transNot!tr )( source.lhs, dest );
		} else static if( op == Operation.MatMatProd ) {
			// copy from a matrix multiplication - fwd to hlMatMult
			static if( tr )
				hlMatMat!( Transpose.yes, Transpose.yes )( One!T, source.rhs, source.lhs, Zero!T, dest );
			else
				hlMatMat!( Transpose.no, Transpose.no )( One!T, source.lhs, source.rhs, Zero!T, dest );
		} else static if( op == Operation.RowMatProd ) {
			hlMatCol!( transNot!tr, transNot!tr )( One!T, source.rhs, source.lhs, Zero!T, dest );
		} else static if( op == Operation.MatColProd ) {
			hlMatCol!( tr, tr )( One!T, source.lhs, source.rhs, Zero!T, dest );
		} else { 
			// can't solve the expression - use a temporary by calling eval
			hlCopy!tr( eval(source), dest );
		}
	} else static if( srcClosure == Closure.RowVector || srcClosure == Closure.ColumnVector ) {
		// fallback vector copy
		debug( blasCalls ) { write( "fb_copy( ", to!string(source), ", ", to!string(dest), " ) => " ); }
		assert( source.length == dest.length, "Length mismatch in fallback vector assignment." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i )
			dest[ i ] = source[ i ];
		debug( blasCalls ) writeln( to!string(dest) );
	} else static if( srcClosure == Closure.Matrix ) {
		// fallback matrix copy - copy major-by-major
		static if( dest.storageOrder == StorageOrder.RowMajor ) {
			size_t n = dest.rows;
			for( auto i = 0; i < n ; ++ i )
				hlCopy( source.row(i), dest.row( i ) );
		} else {
			size_t n = dest.columns;
			for( auto i = 0; i < n ; ++ i )
				hlCopy( source.column(i), dest.column( i ) );
		}
	}
		
}

void hlAxpy( Transpose tr = Transpose.no, Scalar, Source, Dest )( auto ref Scalar alpha, auto ref Source source, auto ref Dest dest ) {
	alias BaseElementType!Source                  T;
	alias transposeClosure!(closureOf!Source, tr) srcClosure;
	alias closureOf!Dest                          dstClosure;
	enum multClosure = operationClosures[ operationOf!("*",srcClosure,closureOf!Scalar) ];
	
	// check correct operation
	static assert( isLiteral!Dest, Dest.stringof ~ " is not an lvalue." );
	static assert( is( T : BaseElementType!Scalar ),
		"Types in " ~ to!string(srcClosure) ~ "-Scalar multiplication do not match (" ~
		T.stringof ~ " and " ~ (BaseElementType!Scalar).stringof ~ ")."
	);
	
	static assert( is( T : BaseElementType!Dest ),
		"Types in addition do not have the same element type (" ~
		T.stringof ~ " and " ~ (BaseElementType!Dest).stringof ~ ")."
	);
	
	static assert( multClosure == dstClosure,
		"Incompatible types for addition '" ~ to!string(multClosure) ~
		"' and '" ~ to!string(dstClosure) ~ "'."
	);
	
	static assert( closureOf!Scalar == Closure.Scalar,
		"Cannot perform " ~ to!string(closure) ~ "-" ~ to!string(closureOf!Scalar) ~
		"multiplication in place."
	);
	
	static if( __traits( compiles, dest.storage.axpy( alpha, source.storage ) ) ) {
		dest.storage.axpy( alpha, source.storage );
	} else static if( __traits( compiles, source.storage.axpyRight( alpha, dest.storage ) ) ) {
		source.storage.axpyRight( alpha, dest.storage );
	} else static if( __traits( compiles, dest.storage.axpy( eval( alpha ), source.storage ) ) ) {
		dest.storage.axpy( eval( alpha ), source.storage );
	} else static if( __traits( compiles, source.storage.axpyRight( eval( alpha ), dest.storage ) ) ) {
		source.storage.axpyRight( eval( alpha ), dest.storage );
	} else static if( isExpression!Source ) {
		enum op = Source.operation;
		static if( isScalingOperation!op ) {
			// axpy with a scaling operation - multiply alpha by the scalar
			hlAxpy!tr( alpha * source.rhs, source.lhs, dest );
		} else static if( isAddition!op ) {
			// axpy with an addition - distribute the scalar between the operands
			hlAxpy!tr( alpha, source.lhs, dest );
			hlAxpy!tr( alpha, source.rhs, dest );
		} else static if( isTransposition!op ) {
			// axpy with transposition - recurse with the transposition flag negated
			hlAxpy!( transNot!tr )( alpha, source.lhs, dest );
		} else static if( op == Operation.MatMatProd ) {
			// axpy from a matrix multiplication - fwd to hlMatMult
			static if( tr )
				hlMatMat!( Transpose.yes, Transpose.yes )( One!T, source.rhs, source.lhs, One!T, dest );
			else
				hlMatMat!( Transpose.no, Transpose.no )( One!T, source.lhs, source.rhs, One!T, dest );
		} else static if( op == Operation.RowMatProd ) {
			hlMatCol!( transNot!tr, transNot!tr )( One!T, source.rhs, source.lhs, One!T, dest );
		} else static if( op == Operation.MatColProd ) {
			hlMatCol!( tr, tr )( One!T, source.lhs, source.rhs, One!T, dest );
		} else {
			// if we can't solve the expression then use a temporary by calling eval
			hlAxpy!tr( alpha, eval(source), dest );
		}
	} else static if( srcClosure == Closure.RowVector || srcClosure == Closure.ColumnVector ) {	
		// fallback vector axpy
		auto alphaValue = eval( alpha );
		debug( blasCalls ) {
			auto dummy = dest.data;
			write( "fb_axpy( ", alphaValue, ", ", to!string(source), ", ", to!string(dest), " ) => " );
		}
		assert( source.length == dest.length, "Length mismatch in fallback vector addition." );
		auto n = dest.length;
		for( size_t i = 0 ; i < n ; ++ i )
			dest[ i ] += source[ i ] * alphaValue;	
		debug( blasCalls ) writeln( to!string(dest) );
	} else static if( srcClosure == Closure.Matrix ) {
		// fallback matrix axpy - axpy major-by-major
		auto alphaValue = eval( alpha );
		static if( dest.storageOrder == StorageOrder.RowMajor ) {
			size_t n = dest.rows;
			for( auto i = 0; i < n ; ++ i )
				hlAxpy( alpha, source.row(i), dest.row( i ) );
		} else {
			size_t n = dest.columns;
			for( auto i = 0; i < n ; ++ i )
				hlAxpy( alpha, source.column(i), dest.column( i ) );
		}
	} else static assert( false );
}

void hlMatMat( Transpose transA = Transpose.no, Transpose transB = Transpose.no, Alpha, Beta, A, B, Dest )( auto ref Alpha alpha, auto ref A a, auto ref B b, auto ref Beta beta, auto ref Dest dest ) {
	alias BaseElementType!Dest T;
	
	static assert( isLiteral!Dest, Dest.stringof ! " is not an lvalue." );
	static assert( is( T : BaseElementType!Alpha ) && is( T : BaseElementType!Beta ) && is( T : BaseElementType!A ) && is( T : BaseElementType!B ) );
	static assert( closureOf!A == Closure.Matrix && closureOf!B == Closure.Matrix && closureOf!Alpha == Closure.Scalar && closureOf!Beta == Closure.Scalar );
	
	static if( __traits( compiles, dest.storage.mm( alpha, a.storage, b.storage, beta ) ) ) {
		dest.storage.mm( alpha, a.storage, b.storage, beta );
	} else static if( __traits( compiles, dest.storage.mm( eval(alpha), a.storage, b.storage, eval(beta) ) ) ) {
		auto alphaValue = eval( alpha );
		auto betaValue  = eval( beta  );
		dest.storage.mm( alphaValue, a.storage, b.storage, betaValue );
	} else static if( isExpression!A ) {
		static if( A.operation == Operation.MatScalProd ) {
			hlMatMat!(transA,transB)( alpha * a.rhs, a.lhs, b, beta, dest ); 
		} else static if( A.operation == Operation.MatTrans ) {
			hlMatMat!( transNot!transA, transB )( alpha, a.lhs, b, beta, dest );
		} else {
			hlMatMat!( transA, transB )( alpha, eval(a), b, beta, dest );
		}
	} else static if( isExpression!B ) {
		static if( B.operation == Operation.MatScalProd ) {
			hlMatMat!( transA, transB )( alpha * b.rhs, a, b.lhs, beta, dest ); 
		} else static if( B.operation == Operation.MatTrans ) {
			hlMatMat!( transA, transNot!transB )( alpha, a, b.lhs, beta, dest );
		} else {
			hlMatMat!( transA, transB )( alpha, a, eval(b), beta, dest );
		}
	} else {
		auto vAlpha = eval( alpha ), vBeta  = eval( beta  );
		fbmm!(transA,transB)( vAlpha, a, b, vBeta, dest );
	}
}

void hlMatCol( Transpose transM, Transpose transV, Alpha, Beta, Mat, Vec, Dest )( auto ref Alpha alpha, auto ref Mat mat, auto ref Vec vec, auto ref Beta beta, auto ref Dest dest ) {
	
	enum vecClosure = transposeClosure!( Vec.closure, transV );
	pragma( msg, vecClosure );
	
	alias BaseElementType!Alpha T;
	
	static assert( vecClosure       == Closure.ColumnVector );
	static assert( closureOf!Alpha  == Closure.Scalar );
	static assert( closureOf!Beta   == Closure.Scalar );
	static assert( closureOf!Mat    == Closure.Matrix );
	static assert( is( BaseElementType!Dest : T ) );
	static assert( is( BaseElementType!Mat  : T ) );
	static assert( is( BaseElementType!Vec  : T ) );
	static assert( is( BaseElementType!Beta : T ) );
	static assert( isLiteral!Dest );
	
	static if( __traits( compiles, mat.storage.mv!(transM,transV)( alpha, vec.storage, beta, dest ) ) ) {
		mat.storage.mv!(transM,transV)( alpha, vec.storage, beta, dest );
	} else static if( __traits( compiles, mat.storage.mv!(transM,transV)( eval( alpha ), vec.storage, eval(beta), dest ) ) ) {
		mat.storage.mv!(transM,transV)( eval( alpha ), vec.storage, eval(beta), dest );
	} else static if( __traits( compiles, vec.storage.mvRight!(transM,transV)( alpha, mat.storage, beta, dest ) ) ) {
		mat.storage.mvRight!(transM,transV)( alpha, mat.storage, beta, dest );
	} else static if( __traits( compiles, vec.storage.mvRight!(transM,transV)( eval( alpha ), mat.storage, eval( beta ), dest ) ) ) {
		mat.storage.mvRight!(transM,transV)( eval( alpha ), vec.storage, eval( beta ), dest );
	} else static if( isExpression!Mat ) {
		static if( Mat.operation == Operation.MatScalProd ) {
			hlMatCol!( transM, transV )( alpha * mat.rhs, mat.lhs, vec, beta, dest );
		} else static if( Mat.operation == Operation.MatTrans ) {
			hlMatCol!( transNot!transM, transV )( alpha, mat.lhs, vec, beta, dest );
		} else {
			hlMatCol!( transM, transV )( alpha, eval(mat), vec, beta, dest );	
		}
	} else static if( isExpression!Vec ) {
		static if( Vec.operation == Operation.RowScalProd || Vec.operation == Operation.ColScalProd ) {
			hlMatCol!( transM, transV )( alpha * vec.rhs, mat, vec.lhs, beta, dest );
		} else static if ( Vec.operation == Operation.RowTrans || Vec.operation == Operation.ColTrans ) {
			hlMatCol!( transM, transNot!transV )( alpha, mat, vec.lhs, beta, dest );
		} else {
			hlMatCol!( transM, transV )( alpha, mat, eval(vec), beta, dest );
		}
	} else {
		auto vAlpha = eval( alpha );
		auto vBeta  = eval( beta );
		fbmv!(transM, transV)( vAlpha, mat, vec, vBeta, dest );
	}
}

auto hlDot( Transpose transA = Transpose.no, Transpose transB = Transpose.no, A, B )( auto ref A a, auto ref B b ) {
	alias BaseElementType!A T;
	
	static if( __traits( compiles, a.storage.dot!( transXor!(transA,transB) )( b.storage ) ) ) {
		auto r = a.storage.dot!( transXor!(transA,transB) )( b.storage );
		static if( transB ) return gconj( r );
		else                return r;
	} else static if( __traits( compiles, b.storage.dot!( transXor!(transA,transB) )( a.storage ) ) ) {
		auto r = b.storage.dot!( transXor!(transA,transB) )( a.storage );
		static if( transA ) return gconj( r );
		else                return r;
	} else static if( isExpression!A ) {
		static if( A.operation == Operation.RowScalProd || A.operation == Operation.ColScalProd ) {
			return hlDot!(transA,transB)(a.lhs.lhs, b) * eval(expr.lhs.rhs);
		} else static if( A.operation == Operation.RowTrans || A.operation == Operation.ColTrans ) {
			return hlDot!( transNot!transA, transB )( a.lhs, b );
		}  else {
			return hlDot!(transA,transB)( eval( a ), b );
		}
	} else static if( isExpression!B ) {
		return hlDot!(transB,transA)( b, a );
	} else {
		return fbdot!(transA, transB)( a, b );
	}
}

auto eval( T )( auto ref T value ) if( isLiteral!T ) {
	return value;	
}

auto eval(string op_, Lhs, Rhs)( auto ref Expression!(op_, Lhs, Rhs) expr )
	if( expr.closure == Closure.Scalar && (
		expr.operation == Operation.ScalScalSum  ||
		expr.operation == Operation.ScalScalProd ||
		expr.operation == Operation.ScalScalDiv  ||
		expr.operation == Operation.ScalScalSub  ||  
		expr.operation == Operation.ScalScalPow
	)) {
	return mixin("eval( expr.lhs )" ~ op_ ~ "eval( expr.rhs )");
}


auto eval( E )( auto ref E expr ) if( isExpression!E && E.closure != Closure.Scalar ) {
	pragma( msg, ExpressionResult!E );
	auto r = ExpressionResult!E( expr );
	return r;
}

auto eval(string op_, Lhs, Rhs)( auto ref Expression!(op_, Lhs, Rhs) expr )
	if( expr.operation == Operation.DotProd ) {
	return hlDot(expr.lhs, expr.rhs);
}

template isPrimitiveLiteral( T ) {
    enum isPrimitiveLiteral = (is( T : BaseElementType!T[] ) || is( T : BaseElementType!T[][] ) || is( T : BaseElementType!T )) && isFortranType!(BaseElementType!T);
}

//template StorageOf( T ) {
//    static if( isPrimitiveLiteral!T )
//        alias T StorageOf;
//    else
//        alias T.Storage StorageOf;
//}
//