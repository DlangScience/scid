/** Common functions used to compute matrix and vector operations.
    
    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/

module scid.ops.common;

import std.complex;
import std.math     : conj;
import std.typecons;

import scid.matvec;
import scid.ops.eval;
import scid.common.traits;

import scid.bindings.blas.dblas;
import scid.bindings.lapack.dlapack;

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
	void scaledAddition( S )( ElementType alpha, auto ref S rhs ) if( isStridedVectorStorage!(S, ElementType) ) {
		stridedScaledAddition( alpha, rhs, this );
	}
	
	void scale( ElementType alpha ) {
		stridedScaling( alpha, this );	
	}
	
	auto dot( Transpose tr, Other )( auto ref Other other ) if( isStridedVectorStorage!(Other, ElementType) ) {
		return stridedDot!tr( this, other );	
	}
}

/** Mixin template that provides scaledAddition() and scale() for general matrix storage types. */
mixin template GeneralMatrixScalingAndAddition() {
	import scid.bindings.lapack.dlapack;
	
	void scaledAddition( S )( ElementType alpha, auto ref S rhs ) if( isGeneralMatrixStorage!(S, ElementType) ) {
		assert( rhs.rows == this.rows && rhs.columns == this.columns, "Matrix size mismatch in addition." );
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaledAddition!(S.storageOrder, storageOrder)(
			this.rows, this.columns,
			alpha,
			rhs.cdata, rhs.leading,
			this.data, this.leading
		);
	}
	
	void scale( ElementType alpha ) {
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaling!storageOrder( this.rows, this.columns, alpha, this.data, this.leading );
	}
	
	void solveRight( Transpose trans, Dest )( auto ref Dest dest ) if( isStridedVectorStorage!Dest ) {
		static if( transposeStorageOrder!(storageOrder, trans) != StorageOrder.ColumnMajor ) {
			static assert( false, "Non column-major matrix or row vector for SolveVector." );
		}
		auto ipiv = new int[ this.rows ];
		auto cpy = this;
		int info;
		if( dest.stride == 1 )
			gesv( cpy.rows, dest.length, cpy.data, cpy.leading, ipiv.ptr, dest.data, dest.length, info );
		else {
			static assert( is(Dest.Referenced), "sv: Stride != 1 in non-view type." );
			Dest.Referenced destValue = dest;
			assert( destValue.stride == 1 );
			gesv( cpy.rows, dest.length, cpy.data, cpy.leading, ipiv.ptr, destValue.data, dest.length, info );
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

/** Specialized copy for strided vector storage types. */
void stridedCopy( Source, Dest )( auto ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,BaseElementType!Dest) && isStridedVectorStorage!Dest ) {
	dest.resize( source.length, null );
	copy( dest.length, source.cdata, source.stride, dest.data, dest.stride );
}

/** Specialized scaled addition (axpy) for strided vector storage types. */
void stridedScaledAddition( Scalar, Source, Dest )( Scalar alpha, ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,Scalar) && is( Scalar : BaseElementType!Dest ) ) {
	assert( source.length == dest.length, "Length mismatch in strided addition." );
	axpy( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
}

/** Specialized scaling for strided vector storage types. */
void stridedScaling( Scalar, Dest )( Scalar alpha, auto ref Dest dest ) if( isStridedVectorStorage!Dest && is( Scalar : BaseElementType!Dest ) ) {
	scal( dest.length, alpha, dest.data, dest.stride );
}

/** Specialized dot for strided vector storage types. */
auto stridedDot( Transpose transA, A, B )( auto ref A a, auto ref B b ) if( isStridedVectorStorage!(A,BaseElementType!B) && isStridedVectorStorage!B ) {
	assert( a.length == b.length, "Length mismatch in strided dot." );
	static if( !isComplex!(BaseElementType!A) )
		auto r = dot( a.length, a.cdata, a.stride, b.cdata, b.stride );	
	else static if( transA )
		auto r = dotc( a.length, a.cdata, a.stride, b.cdata, b.stride );
	else
		auto r = dotu( a.length, a.cdata, a.stride, b.cdata, b.stride );	
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
		// zero everything
		scal( rows * cols, alpha, data, 1 );
	} else {
		// zero by-major
		auto e = data + major * leading;
		for( ; data < e ; data += leading )
			scal( minor, alpha, data, 1 );
	}
}

/** Specialized copying for general matrix storage types. */
void generalMatrixCopy( StorageOrder srcOrder, StorageOrder dstOrder, E )
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

/** Specialized scaled addition (axpy) for general matrix storage types. */
void generalMatrixScaledAddition( StorageOrder srcOrder, StorageOrder dstOrder, E )
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