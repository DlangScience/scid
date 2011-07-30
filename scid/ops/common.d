/** Common functions used to compute matrix and vector operations.
    
    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/

module scid.ops.common;

public import scid.ops.expression;

import std.complex;
import std.math;
import std.typecons;

import scid.matvec;
import scid.ops.eval;
import scid.common.traits;

import scid.blas;
import scid.lapack;

/** Generic complex conjugate that works on both builtin complex numbers and on std.Complex */
auto gconj( T )( T z ) if( is( T : cdouble ) || is( T : cfloat ) ) {
	return conj( z );
}
/// ditto
auto gconj( T )( T z ) if( is( T : Complex!double ) || is( T : Complex!float ) ) {
	return z.conj;
}

/** Inverses the value of a Transpose value. */
template transNot( Transpose tr )             { enum transNot = cast(Transpose) !tr; }

/** Performs exclusive or on Transpose values. */
template transXor( Transpose a, Transpose b ) { enum transXor = cast(Transpose) (a ^ b); }

/** Gets the Transposed type of a given Matrix/Vector storage. If the given type is ref-counted then the result will
    be ref-counted as well.
*/
template TransposedOf( T ) {
	static if( is( T.Transposed ) )
		alias T.Transposed TransposedOf;
	else static if( is( T E : RefCounted!(E,x), uint x ) ) {
		alias RefCounted!(TransposedOf!E,cast(RefCountedAutoInitialize)x) TransposedOf;
	} else static assert( false, T.stringof ~ " has no transpose." );
}

/** Transposes a closure (only affects RowVector <-> ColumnVector) with a given Transpose value. */
template transposeClosure( Closure A, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( A == Closure.RowVector )         enum transposeClosure = Closure.ColumnVector;
		else static if( A == Closure.ColumnVector ) enum transposeClosure = Closure.RowVector;
		else                                        enum transposeClosure = A;
	} else
		enum transposeClosure = A;
}

/** Transposes a vector type (Row <-> Column) based on the given Transpose value. */
template transposeVectorType( VectorType v, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( v == VectorType.Row )         enum transposeVectorType = VectorType.Column;
		else static if( v == VectorType.Column ) enum transposeVectorType = VectorType.Row;
		else                                     enum transposeVectorType = v;
	} else
		enum transposeVectorType = v;
}

/** Transposes a storage order (RowMajor <-> ColumnMajor) based on the given Transpose value. */
template transposeStorageOrder( StorageOrder S, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( S == StorageOrder.RowMajor )         enum transposeStorageOrder = StorageOrder.ColumnMajor;
		else static if( S == StorageOrder.ColumnMajor ) enum transposeStorageOrder = StorageOrder.RowMajor;
		else                                            enum transposeStorageOrder = S;
	} else
		enum transposeStorageOrder = S;
}

/** Transposes a vector if the given parameter is Transpose.yes. This simply means conjugate its elements if the
    element type is complex.
*/
void vectorTranspose( Transpose trans, V )( ref V vec ) {
	static if( trans && isComplex!(BaseElementType!V) ) {
		auto n = vec.length;
		foreach( i ; 0 .. n )
			vec[ i ] = gconj( vec[ i ] );
	}
}

/** Mixin template that provides scaledAddition(), scale() and dot() specializations for strided vector storage
    types.
*/
mixin template StridedScalingAdditionDot() {
	import scid.common.traits : isStridedVectorStorage;
	void scaledAddition( Transpose tr = Transpose.no, S )( ElementType alpha, auto ref S rhs ) if( isStridedVectorStorage!(S, ElementType) ) {
		stridedScaledAddition!tr( alpha, rhs, this );
	}
	
	void scale( ElementType alpha ) {
		stridedScaling( alpha, this );	
	}
	
	auto dot( Transpose tr = Transpose.no, Other )( auto ref Other other ) if( isStridedVectorStorage!(Other, ElementType) ) {
		return stridedDot!tr( this, other );	
	}
}

/** Mixin template that provides scaledAddition() and scale() for general matrix storage types. */
mixin template GeneralMatrixScalingAndAddition() {
	private import scid.lapack;
	private import scid.internal.regionallocator;
	
	void scaledAddition( Transpose tr = Transpose.no, S )( ElementType alpha, auto ref S rhs ) if( isGeneralMatrixStorage!(S, ElementType) ) {
		assert( rhs.rows == this.rows && rhs.columns == this.columns, "Matrix size mismatch in addition." );
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaledAddition!tr( alpha, rhs, this );
	}
	
	void scale( ElementType alpha ) {
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaling!storageOrder( this.rows, this.columns, alpha, this.data, this.leading );
	}
	
	void solveRight( Transpose trans, Dest )( auto ref Dest dest ) if( isStridedVectorStorage!Dest ) {
		static if( transposeStorageOrder!(storageOrder, trans) != StorageOrder.ColumnMajor ) {
			static assert( false, "Non column-major matrix or row vector for SolveVector." );
		}
		auto alloc = newRegionAllocator();
		auto ipiv = alloc.uninitializedArray!(int[])( this.rows );
		auto cpy = this;
		int info;
		if( dest.stride == 1 )
			lapack.gesv( cpy.rows, dest.length, cpy.data, cpy.leading, ipiv.ptr, dest.data, dest.length, info );
		else {
			static if( is( typeof(eval(dest)) : Dest ) )
				assert( false, "sv: Stride !=1 in non-view type." );
			else {
				auto destValue = eval( dest );
				assert( destValue.stride == 1 );
				lapack.gesv( cpy.rows, dest.length, cpy.data, cpy.leading, ipiv.ptr, destValue.data, dest.length, info );
			}
		}
		assert( info <= 0, "sv: Singular matrix in inversion." );
			
	}
	
	void solveRight( Transpose trans, Dest )( auto ref Dest dest ) if( isGeneralMatrixStorage!Dest ) {
		static if( transposeStorageOrder!(storageOrder, trans) != StorageOrder.ColumnMajor &&
				   storageOrderOf!Dest == StorageOrder.RowMajor ) {
			static assert( false, "Non column-major matrix or row vector for SolveVector." );
		}
		
		auto ipiv = new int[ this.rows ];
		auto cpy = this;
		int info;
		// gesv( cpy.rows, dest.length, thisCopy.data, cpy.leading, ipiv.ptr, dest.data, dest.leading, info );
		assert( info <= 0, "sm: Singular matrix in inversion." );
	}
}

/** Compute the dot product of a row and a column of possibly transposed matrices. */
auto rowColumnDot( Transpose transA = Transpose.no, Transpose transB = Transpose.no, A, B )( auto ref A a, size_t i, auto ref B b, size_t j ) {
	static if( !transA && !transB ) {
		return evalDot!( transA, transB )( a.row( i ), b.column( j ) );
	} else static if( !transA && transB ) {
		return evalDot!( transA, transB )( a.row( i ), b.row( j ) );
	} else static if( transA && !transB ) {
		return evalDot!( transA, transB )( a.column( i ), b.column( j ) );
	} else {
		return evalDot!( transA, transB )( b.column( i ), b.row( j ) );
	}
}

/** Specialized copy for strided vector storage types. */
void stridedCopy( Transpose tr, Source, Dest )( auto ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,BaseElementType!Dest) && isStridedVectorStorage!Dest ) {
	dest.resize( source.length, null );
	static if( isComplex!(BaseElementType!Source) && tr ) {
		blas.xcopyc( dest.length, source.cdata, source.stride, dest.data, dest.stride );
	} else
		blas.copy( dest.length, source.cdata, source.stride, dest.data, dest.stride );
	
}

/** Specialized scaled addition (axpy) for strided vector storage types. */
void stridedScaledAddition( Transpose tr, Scalar, Source, Dest )( Scalar alpha, ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,Scalar) && is( Scalar : BaseElementType!Dest ) ) {
	assert( source.length == dest.length, "Length mismatch in strided addition." );
	static if( isComplex!(BaseElementType!Source) && tr )
		blas.xaxpyc( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
	else
		blas.axpy( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
}

/** Specialized scaling for strided vector storage types. */
void stridedScaling( Scalar, Dest )( Scalar alpha, auto ref Dest dest ) if( isStridedVectorStorage!Dest && is( Scalar : BaseElementType!Dest ) ) {
	blas.scal( dest.length, alpha, dest.data, dest.stride );
}

/** Specialized dot for strided vector storage types. */
auto stridedDot( Transpose transA, A, B )( auto ref A a, auto ref B b ) if( isStridedVectorStorage!(A,BaseElementType!B) && isStridedVectorStorage!B ) {
	assert( a.length == b.length, "Length mismatch in strided dot." );
	static if( !isComplex!(BaseElementType!A) )
		auto r = blas.dot( a.length, a.cdata, a.stride, b.cdata, b.stride );	
	else static if( transA )
		auto r = blas.dotc( a.length, a.cdata, a.stride, b.cdata, b.stride );
	else
		auto r = blas.dotu( a.length, a.cdata, a.stride, b.cdata, b.stride );	
	return r;
}

/** Specialized scaling for general matrix storage types. */
void generalMatrixScaling( StorageOrder order, E )( size_t rows, size_t cols, E alpha,  E* data, size_t leading ) {
	static if( order == StorageOrder.RowMajor ) {
		alias rows major;
		alias cols minor;
	} else {
		alias rows minor;
		alias cols major;
	}	 
	
	if( leading == minor ) {
		// scale everything
		blas.scal( rows * cols, alpha, data, 1 );
	} else {
		// scale by-major
		auto e = data + major * leading;
		for( ; data < e ; data += leading )
			blas.scal( minor, alpha, data, 1 );
	}
}

/** Specialized copying for general matrix storage types. */
void generalMatrixCopy( Transpose tr, S, D )( auto ref S source, auto ref D dest ) {
	alias BaseElementType!S T;
	alias transposeStorageOrder!(storageOrderOf!S, tr) srcOrder;
	alias storageOrderOf!D dstOrder;
	
	static if( !isComplex!T ) {
		blas.xgecopy( (srcOrder == dstOrder) ? 't' : 'n',
			dest.rows,
			dest.columns,
			source.cdata,
			source.leading,
			dest.data,
			dest.leading
		);
	} else static if( isComplex!T ) {
		static if( srcOrder == dstOrder ) {
			static if( tr ) {
				blas.xgecopyc( 'n', dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgecopy( 'n', dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		} else {
			static if( tr ) {
				blas.xgecopy( 'c', dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgecopy( 't', dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		}
	}
}

/** Specialized scaled addition (axpy) for general matrix storage types. */
void generalMatrixScaledAddition( Transpose tr, S, D, E )( T alpha, auto ref S source, auto ref D dest ) {
	alias transposeStorageOrder!(storageOrderOf!S, tr) srcOrder;
	alias storageOrderOf!D dstOrder;
	
	static if( !isComplex!T ) {
		blas.xgeaxpy( (srcOrder == dstOrder) ? 't' : 'n',
			dest.rows, dest.columns, alpha,
			source.cdata, source.leading,
			dest.data, dest.leading
		);
	} else static if( isComplex!T ) {
		static if( srcOrder == dstOrder ) {
			static if( tr ) {
				blas.xgeaxpyc( 'n', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgeaxpy( 'n', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		} else {
			static if( tr ) {
				blas.xgeaxpy( 'c', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgeaxpy( 't', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		}
	}
}