module scid.vector;

import scid.storage.array;
import scid.storage.arrayview;
import scid.common.traits;
import scid.internal.expression;
import scid.internal.hlblas;

import std.traits, std.range, std.algorithm, std.conv;

enum VectorType {
	Row, Column
}

enum Transpose : bool {
	no = false,
	yes = true
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
	alias storage                                            this;

	
	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isVector!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage(args);
	}
	
	this( Expr )( auto ref Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	this( A )( auto ref BasicVector!A other ) {
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
		move( rhs.storage, storage );
		return this;
	}
	
	ref typeof(this) opAssign( ElementType[] rhs ) {
		//storage.copyFrom( rhs );
		return this;
	}
	/*
	ref typeof(this) opOpAssign( string op, Rhs )( Rhs rhs ) if( op == "+" || op == "-" ) {
		hlAxpy( signOfOp!(op,ElementType), rhs, this );
		return this;
	}
	
	ref typeof(this) opOpAssign( string op )( ElementType rhs ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1.0 / rhs;
		hlScal( rhs, this );
		return this;
	}
	*/
	typeof( this ) opSlice() {
		return typeof(this)( storage );
	}
	
	typeof( this ) opSlice( size_t start, size_t end ) {
		return typeof(this)( storage.slice( start, end ) );	
	}
	
	bool opEquals( Rhs )( Rhs rhs ) if( isInputRange!Rhs ) {
		size_t i = 0;
		foreach( x ; rhs ) {
			if( i >= length || this[ i ++ ] != x )
				return false;
		}
		
		return true;
	}

	void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
		hlCopy( rhs, this );
	}
	
	void opSliceAssign( Rhs )( Rhs rhs, size_t start, size_t end ) {
		hlCopy( rhs, view( start, end ) );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "-" ) {
		hlAxpy( signOfOp!(op,ElementType), rhs, this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1.0 / rhs;
		hlScal( rhs, this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if( op == "+" || op == "-" ) {
		hlAxpy( signOfOp!(op,ElementType), rhs, view(start,end) );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1.0 / rhs;
		hlScal( rhs, view( start, end ) );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( storage.view( start, end ) );
	}
	
	StridedView view( size_t start, size_t end, size_t stride ) {
		return typeof( return )( storage.view( start, end, stride ) );	
	}
	
	static if( isBidirectionalRange!Storage ) {
		void popFront() { storage.popFront(); }
		void popBack()  { storage.popBack(); }
	}
	
	@property {
		bool                 empty()   const { return storage.empty; }
		ElementType          front()   const { return storage.front; }
		ElementType          back()    const { return storage.back; }
		size_t               length()  const { return storage.length; }
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
	
	static if( vectorType == VectorType.Column ) mixin Literal!( Closure.ColumnVector );
	else                                         mixin Literal!( Closure.RowVector    );
	
	Storage storage;
}

unittest {
	// TODO: Write tests for Vector.
}