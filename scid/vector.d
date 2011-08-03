module scid.vector;

import scid.storage.array;
import scid.storage.arrayview;
import scid.storage.external;
import scid.storage.cowarray;
import scid.common.traits;
import scid.common.meta;
import scid.ops.eval;

import std.traits, std.range, std.algorithm, std.conv;

import scid.internal.assertmessages;

enum VectorType {
	Row, Column
}

template Vector( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
	
	static if( isScalar!ElementOrStorage )
		alias BasicVector!( ArrayStorage!( ElementOrStorage, vectorType ) ) Vector;
	else
		alias BasicVector!( ElementOrStorage )              Vector;
}

template VectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
			
	alias BasicVector!( ArrayViewStorage!( ElementOrStorage, vectorType ) ) VectorView;
}

template StridedVectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
	
	alias BasicVector!( StridedArrayViewStorage!( ElementOrStorage, vectorType ) ).View StridedVectorView;
}

template ExternalVectorView( ElementOrContainer, VectorType vectorType = VectorType.Column )
		if( isScalar!( BaseElementType!ElementOrContainer ) ) {
	
	static if( isScalar!ElementOrContainer ) {
		alias BasicVector!(
			ArrayViewStorage!(
				ExternalArray!( ElementOrContainer, CowArrayRef!ElementOrContainer ),
				vectorType
			)
		) ExternalVectorView;
	} else {
		alias BasicVector!(
			ArrayViewStorage!(
				ExternalArray!( BaseElementType!ElementOrContainer, ElementOrContainer ),
				vectorType
			)
		) ExternalVectorView;
	}
}

template signOfOp( string op, T ) {
	static if( op == "+" )
		enum T signOfOp = One!T;
	else static if( op == "-" )
		enum T signOfOp = MinusOne!T;
}

struct BasicVector( Storage_ ) {
	alias BaseElementType!Storage                          ElementType;
	alias Storage_                                         Storage;
	alias BasicVector!( typeof(Storage.init.view(0,0)) )   View;
	alias BasicVector!( Storage.Transposed )               Transposed;
	alias storage                                          this;
	
	static if( is( Storage.Temporary ) )
		alias BasicVector!( Storage.Temporary ) Temporary;
	else
		alias typeof( this ) Temporary;
	
	static if( is( typeof(Storage.vectorType) ) )
		alias Storage.vectorType vectorType;
	else
		alias VectorType.Column vectorType;
	
	/** Whether the storage can be resized. */
	enum isResizable = is( typeof( Storage.init.resize(0) ) );
	
	static if( is( typeof(Storage.init.view(0,0,0)) R ) )
		alias BasicVector!R StridedView;
	
	//static assert( isVectorStorage!Storage );

	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isVector!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage(args);
	}
	
	this( Expr )( Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	this( A )( BasicVector!A other ) {
		static if( is( A : Storage ) ) storage = other.storage;
		else                           this[] = other;
	}
	
	this()( auto ref Storage stor ) {
		storage = stor;
	}
	
	ElementType opIndex( size_t i ) const {
		return storage.index( i );
	}
	
	void opIndexAssign( ElementType rhs, size_t i ) {
		storage.indexAssign( rhs, i );
	}
	
	void opIndexOpAssign( string op )( ElementType rhs, size_t i ) {
		storage.indexAssign!op( rhs, i );
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		rhs.storage = storage;
		return this;
	}
	
	typeof( this ) opSlice() {
		return typeof(this)( storage );
	}
	
	typeof( this ) opSlice( size_t start, size_t end ) {
		return typeof(this)( storage.slice( start, end ) );	
	}
	
	bool opEquals( Rhs )( Rhs rhs ) const if( isInputRange!Rhs ) {
		size_t i = 0;
		foreach( x ; rhs ) {
			if( i >= length || this[ i ++ ] != x )
				return false;
		}
		
		return true;
	}
	
	/** Resize the vector and leave the memory uninitialized. If not resizeable simply check that the length is
	    correct.
	*/
	void resize( size_t newLength, void* ) {
		static if( isResizable ) {
			storage.resize( newLength, null );
		} else {
			assert( length == newLength,
				lengthMismatch_( newLength ) );
		}
	}
	
	/** Resize the vectors and set all the elements to zero. If not resizeable, check that the length is correct
	    and just set the elements to zero.
	*/
	void resize( size_t newLength ) {
		static if( isResizable ) {
			storage.resize( newLength );
		} else {
			this.resize( newLength, null );
			evalScaling( Zero!ElementType, this );
		}
	}

	void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
		static if( is( Rhs E : E[] ) && isConvertible( E, ElementType  ) )
			evalCopy( BasicVector(rhs), this );
		else
			evalCopy( rhs, this );
		
	}
	
	void opSliceAssign( Rhs )( Rhs rhs, size_t start, size_t end ) {
		auto v = view( start, end );
		v[] = rhs;
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "-" ) {
		evalScaledAddition( signOfOp!(op,ElementType), rhs, this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( (op == "*" || op == "/") && isConvertible!(Rhs,ElementType) ) {
		static if( op == "/" )
			rhs = One!ElementType / rhs;
		evalScaling( to!ElementType(rhs), this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if( op == "+" || op == "-" ) {
		auto v = view( start, end );
		evalScaledAddition( signOfOp!(op,ElementType), rhs, v );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if( (op == "*" || op == "/") && isConvertible!(Rhs,ElementType) ) {
		auto v = view( start, end );
		auto s = to!ElementType(rhs);
		static if( op == "/" )
			s = One!ElementType / s;
		evalScaling( s, v );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( storage.view( start, end ) );
	}
	
	static if( is( StridedView ) ) {
		StridedView view( size_t start, size_t end, size_t stride ) {
			return typeof( return )( storage.view( start, end, stride ) );	
		}
	}
	
	static if( isInputRange!Storage ) {
		void popFront() { storage.popFront(); }
		void popBack()  { storage.popBack(); }
	}
	
	@property {
		bool empty() const {
			static if( is(typeof(Storage.init.empty)) )
				return storage.empty;
			else
				return storage.length == 0;
		}
		
		size_t length() const {
			return storage.length;
		}
		
		ElementType front() const
		in {
			assert( !empty, emptyMsg_("front") );	
		} body {
			static if( is( typeof(Storage.init.front) ) )
				return storage.front;
			else
				return storage.index( 0 );
		}
		
		ElementType back() const
		in {
			assert( !empty, emptyMsg_("back") );
		} body {
			static if( is( typeof(Storage.init.back) ) )
				return storage.back;
			else
				return storage.index( storage.length - 1 );
		}
	}
	
	string toString() const {
		if( empty )
			return "[]";
		
		auto r = appender!string("[");
		r.put( to!string( this[ 0 ] ) );
		foreach( i ; 1 .. length ) {
			r.put( ", " );
			r.put( to!string( this[i] ) );
		}
		r.put( "]" );
		return r.data();
	}
	
	alias toString pretty;
	
	static if( vectorType == VectorType.Column ) mixin Operand!( Closure.ColumnVector );
	else                                         mixin Operand!( Closure.RowVector    );
	
	template Promote( T ) {
		static if( isVector!T ) {
			alias BasicVector!( Promotion!(Storage,T.Storage) ) Promote;
		} else static if( isMatrix!T ) {
			alias BasicVector!( Promotion!(Storage,T.Storage) ) Promote;
		} else static if( isScalar!T ) {
			alias BasicVector!( Promotion!(Storage,T) ) Promote;
		}
	}
	
	Storage storage;

private:
	mixin ArrayErrorMessages;	
}

unittest {
	// TODO: Write tests for Vector.
}