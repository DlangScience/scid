module scid.storage.constant;

import scid.matvec;
import scid.common.storagetraits;
import scid.ops.common;
import scid.ops.eval;

import std.conv;

template RelatedConstantType( T ) {
	static if( closureOf!T == Closure.Scalar )
		alias T RelatedConstantType;
	else static if( closureOf!T == Closure.RowVector ) {
		alias BasicVector!( ConstantStorage!( BaseElementType!T, VectorType.Row ) ) RelatedConstantType;
	} else static if( closureOf!T == Closure.ColumnVector ) {
		alias BasicVector!( ConstantStorage!( BaseElementType!T, VectorType.Column ) ) RelatedConstantType;
	} else static if( closureOf!T == Closure.Matrix ) {
		static if( storageOrderOf!T == StorageOrder.RowMajor ) {
			alias BasicMatrix!( ConstantStorage!( BaseElementType!T, VectorType.Row ) ) RelatedConstantType;
		} else {
			alias BasicMatrix!( ConstantStorage!( BaseElementType!T, VectorType.Column ) ) RelatedConstantType;
		}
	}
}

auto relatedConstant( T, V )( V value_, auto ref T toWhom ) if(closureOf!V == Closure.Scalar) {
	alias BaseElementType!T E;
	
	static assert( isConvertible!( BaseElementType!V, E ),
		"Incompatible constant type (" ~ typeof(BaseElementType!V).stringof ~ " ) for " ~ T.stringof );
	
	auto value = to!E( eval(value_) );
	static if( closureOf!T == Closure.Scalar )
		return to!T( value );
	else static if( closureOf!T == Closure.RowVector ) {
		return BasicVector!( ConstantStorage!( E, VectorType.Row ) )( value, toWhom.length );
	} else static if( closureOf!T == Closure.ColumnVector ) {
		return BasicVector!( ConstantStorage!( E, VectorType.Column ) )( value, toWhom.length );
	} else static if( closureOf!T == Closure.Matrix ) {
		static if( storageOrderOf!T == StorageOrder.RowMajor ) {
			alias BasicMatrix!( ConstantStorage!( E, VectorType.Row ) ) R;
		} else {
			alias BasicMatrix!( ConstantStorage!( E, VectorType.Column ) ) R;
		}
		
		return R( value, toWhom.rows, toWhom.columns );
	}
}

struct ConstantStorage( T, alias OrderOrVectorType ) {
	alias T     ElementType;
	alias T     ContainerRef;
	alias typeof(this) Transposed;
	
	static if( is( typeof( OrderOrVectorType ) : VectorType ) ) {
		alias OrderOrVectorType vectorType;
		static if( vectorType == VectorType.Row )
			enum storageOrder = StorageOrder.RowMajor;
		else
			enum storageOrder = StorageOrder.ColumnMajor;
	} else static if( is( typeof( OrderOrVectorType ) : StorageOrder ) ) {
		alias OrderOrVectorType storageOrder;
		static if( storageOrder == StorageOrder.RowMajor )
			enum vectorType = VectorType.Row;
		else
			enum vectorType = VectorType.Column;
	}
	
	static if( vectorType == VectorType.Column ) {
		alias rows minor;
		alias columns major;
	} else {
		alias columns minor;
		alias rows major;
	}
	
	
	this( T value, size_t length ) {
		this.value  = value;
		this.length = length;
	}
	
	this( T value, size_t r, size_t c ) {
		this.value   = value;
		this.rows    = r;
		this.columns = c;
	}
	
	T index( size_t i, size_t j = 0 ) const {
		return value;
	}
	
	void indexAssign( string op = "" )( T rhs, size_t i, size_t j = 0 ) {
		value = rhs;
	}
	
	void forceRefAssign( ref typeof( this ) rhs ) {
		this = rhs;
	}
	
	ref typeof( this ) opAssign( ref typeof( this ) rhs ) {
		value   = rhs.value;
		rows    = rhs.rows;
		columns = rhs.columns;
		return this;
	}
	
	typeof( this ) view( size_t start, size_t end, size_t stride = 1 ) {
		size_t len = end - start;
		len = (len / stride) + (len % stride == 0);
		return typeof( return )( value, len );
	}
	
	typeof( this ) view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( value, rowEnd - rowStart, colEnd - colStart );
	}
	
	typeof( this ) slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return view( rowStart, rowEnd, colStart, colEnd );
	}
	
	typeof( this ) slice( size_t start, size_t end ) {
		return view( start, end );
	}
	
	void resize( size_t newRows, void* arg = null ) {
		rows = newRows;
	}
	
	void resize( size_t newRows, size_t newCols, void *arg = null ) {
		rows = newRows;
		columns = newCols;
	}
	
	auto row( size_t ) {
		return BasicVector!(ConstantStorage!(T, VectorType.Row))( value, columns );
	}
	
	auto column( size_t ) {
		return BasicVector!(ConstantStorage!(T, VectorType.Column))( value, rows );
	}
	
	private void commonOp_( string op, Dest )( ElementType v, auto ref Dest lhs ) {
		static if( isStridedVectorStorage!Dest ) {
			if(	lhs.stride == 1 ) {
				mixin( "lhs.data[ 0 .. lhs.length ] " ~ op ~ " v;");
			} else {
				auto lhsPtr = lhs.data;
				auto lhsInc = lhs.stride;
				auto lhsEnd = lhsPtr + lhs.length * lhsInc;
				for( ; lhsPtr != lhsEnd ; lhsPtr += lhsInc )
					mixin( "(*lhsPtr) " ~ op ~ " v;");
			}
		} else static if( isGeneralMatrixStorage!Dest ) {
			enum rowMajor = storageOrderOf!Dest == StorageOrder.RowMajor;
			
			static if( rowMajor ) auto lhsMaj = lhs.rows,    lhsMin = lhs.columns;
			else                  auto lhsMaj = lhs.columns, lhsMin = lhs.rows;
			
			auto lhsPtr = lhs.data;
			auto lhsInc = lhs.leading;
			
			if( lhsInc == lhsMin )
				mixin( "lhsPtr[ 0 .. lhsMaj * lhsInc ] " ~ op ~ " v;");
			else {
				auto lhsEnd = lhsPtr + lhsInc * lhsMaj;
			
				for( ; lhsPtr != lhsEnd ; lhsPtr += lhsInc )
					mixin( "lhsPtr[ 0 .. lhsMin ] " ~ op ~ " v;");
			}
		} else static if( is( Dest : typeof( this ) ) ) {
			mixin( "lhs.value " ~ op ~ " v;" );
		} else
			static assert( false );
	}

	void copyRight( Transpose tr = Transpose.no, Dest )( auto ref Dest lhs ) {
		commonOp_!"="( value, lhs );
	}
	
	void scaledAdditionRight( Transpose tr = Transpose.no, Dest )( ElementType alpha, auto ref Dest lhs ) {
		commonOp_!"+="( value * alpha, lhs );
	}
	
	void scale( ElementType alpha ) {
		value *= alpha;	
	}
	
	void invert() {
		assert( false, "inv: Constant matrices are always singular." );
	}
	
	void solveRight() {
		invert();
	}
	
	template Promote( Other ) {
		static if( is( typeof( this ) : Other ) ) {
			alias typeof( this ) Other;
		}
	}
	
	
	@property {
		bool   empty()  const { return major == 0; }
		
		alias major length;
		alias value front;
		alias value back;
	}
	
	void popFront() { -- major; }
	void popBack()  { -- major; }
	
	ElementType value;
	size_t rows, columns;
	

}
