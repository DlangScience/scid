module scid.vector;

import scid.storage.array;
import scid.storage.arrayview;
import scid.common.traits;
import scid.ops.expression;
import scid.ops.eval;
import scid.matrix;

import std.traits, std.range, std.algorithm, std.conv;

enum VectorType {
	Row, Column
}

template Vector( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	static if( isFortranType!ElementOrStorage )
		alias BasicVector!( ArrayStorage!( ElementOrStorage, vectorType ) ) Vector;
	else
		alias BasicVector!( ElementOrStorage )              Vector;
}

template VectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
			
	alias BasicVector!( ArrayViewStorage!( ElementOrStorage, vectorType ) ) VectorView;
}

template StridedVectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	alias BasicVector!( StridedArrayViewStorage!( ElementOrStorage, vectorType ) ).View StridedVectorView;
}

auto vectorWithStorage( S )( auto ref S storage ) {
	return BasicVector!S( storage );	
}

template isVector( T ) {
	static if( is( typeof( T.init[0]          ) ) &&
			   is( typeof( T.init[0..1]       ) ) &&
			   is( typeof( T.init.storage     ) ) &&
			   is( typeof( T.init.view(0,0)   ) ) &&
			   is( typeof( T.init.view(0,0,1) ) ) &&
			   isInputRange!T )
		enum isVector = true;
	else
		enum isVector = false;
}

template signOfOp( string op, T ) {
	static if( op == "+" )
		enum T signOfOp = 1;
	else static if( op == "-" )
		enum T signOfOp = -1;
}

struct BasicVector( Storage_ ) {
	alias BaseElementType!Storage                          ElementType;
	alias Storage_                                         Storage;
	alias Storage.vectorType                               vectorType;
	alias BasicVector!( typeof(Storage.init.view(0,0)) )   View;
	alias BasicVector!( typeof(Storage.init.view(0,0,0)) ) StridedView;
	alias BasicVector!( Storage.Transposed )               Transposed;
	alias storage                                          this;
	
	static if( is( Storage.Temporary ) )
		alias BasicVector!( Storage.Temporary ) Temporary;
	else
		alias typeof( this ) Temporary;
	
	//static assert( isVectorStorage!Storage );
	
	static if( isReference!Storage )
		alias BasicVector!( ReferencedBy!Storage ) Referenced;

	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isVector!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage(args);
	}
	
	this( Expr )( Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	this( A )( BasicVector!A other ) {
		static if( is( A : Storage ) ) move( other.storage, storage );
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
		move( rhs.storage, storage );
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

	void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
		evalCopy( rhs, this );
	}
	
	void opSliceAssign( Rhs )( Rhs rhs, size_t start, size_t end ) {
		auto v = view( start, end );
		evalCopy( rhs, v );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "-" ) {
		evalScaledAddition( signOfOp!(op,ElementType), rhs, this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1.0 / rhs;
		evalScaling( rhs, this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if( op == "+" || op == "-" ) {
		auto v = view( start, end );
		evalScaledAddition( signOfOp!(op,ElementType), rhs, v );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if(op == "*" || op == "/" ) {
		auto v = view( start, end );
		static if( op == "/" ) rhs = 1.0 / rhs;
		evalScaling( rhs, v );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( storage.view( start, end ) );
	}
	
	StridedView view( size_t start, size_t end, size_t stride ) {
		return typeof( return )( storage.view( start, end, stride ) );	
	}
	
	static if( isInputRange!Storage ) {
		void popFront() { storage.popFront(); }
		void popBack()  { storage.popBack(); }
	}
	
	@property {
		bool          empty()   const { return storage.empty; }
		ElementType   front()   const { return storage.front; }
		ElementType   back()    const { return storage.back; }
		size_t        length()  const { return storage.length; }
		typeof(this)* addr()          { return &this; }
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
		static if( is( T S : BasicVector!S ) ) {
			alias BasicVector!( Promotion!(Storage,S) ) Promote;
		} else static if( is( T S : BasicMatrix!S ) ) {
			alias BasicVector!( Promotion!(Storage,S) ) Promote;
		} else static if( isFortranType!T ) {
			alias BasicVector!( Promotion!(Storage,T) ) Promote;
		}
	}
	
	Storage storage;
}

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


unittest {
	// TODO: Write tests for Vector.
}