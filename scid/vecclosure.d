import vecclosure;
import scalclosure;
import expressions;

import scid.matrix;
import scid.common.traits;

/** Convinience function to create a VectorClosure while
	avoiding typing the template parameter.
*/
VectorClosure!T vectorClosure( T )( T value ) {
	return typeof( return )( value );
}

/** Convinience function to create a BinaryExpression
	representing vector addition.
*/
BinaryExpression!( ExpressionKind.VectorSum, L, R ) vectorSum( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}

/** Convinience function to create a VectorUnary
	expression.
*/
VectorUnary!(T) vectorUnary( T )( T vec ) {
	return typeof( return )( vec, ScalarClosure!(BaseElementType!T)(1) );
}

/// ditto
VectorUnary!(T) vectorUnary( T )( T vec, ScalarClosure!(BaseElementType!T) scal ) {
	return typeof( return )( vec, scal );
}

struct VectorClosure( T ) {
	/// Convinence check for same base element type (all the operators require it).
	template hasSameType( alias U ) {
		enum hasSameType = haveSameElementType!(T,U);
	}
	
	/// Type of the enclosed expression.
	alias T						EnclosedType;
	alias BaseElementType!T		ElementType;
	alias expressionKindOf!T	expressionKind;
	alias T.isRow				isRow;
	
	/// The enclosed expression.
	T expr;
	
	/// Construct a closure from an expression.
	this( T expr ) {
		this.expr = expr;
	}
	
	// -- operators for closure arguments --
	/// Vector negation.
	auto opUnary( string op )() if( op == "-" ) {
		static if( expressionKind == ExpressionKind.VectorUnary )
			return vectorClosure( vectorUnary( expr.lhs, -expr.rhs ) );
		else
			return vectorClosure( vectorUnary( expr, ElementType(-1) ) );
	}
	
	/// Vector addition.
	auto opBinary( string op, Rhs )( VectorClosure!Rhs rhs ) if( op == "+" && hasSameType!Rhs ) {
		return vectorClosure( vectorSum( expr, rhs.expr ) );
	}
	
	/// Vector subtraction. Implemented as addition with -1 * rhs;
	auto opBinary( string op, Rhs )( VectorClosure!Rhs rhs ) if( op == "-" && hasSameType!Rhs ) {
		return this + (-rhs);
	}
	
	/// Vector dot product.
	auto opBinary( string op, Rhs )( VectorClosure!Rhs rhs ) if( op == "*" && hasSameType!Rhs && isRow && !rhs.isRow ) {
		return scalarClosure( dotProduct( expr, rhs.expr ) );	
	}
	
	/// Vector-Scalar multiplication.
	auto opBinary( string op, Rhs )( ScalarClosure!Rhs rhs ) if( op == "*" && hasSameType!Rhs ) {
		static if( expressionKind == ExpressionKind.VectorUnary )
			return vectorClosure( vectorUnary( expr.lhs, rhs * expr.rhs ) );
		else
			return vectorClosure( vectorUnary( expr, rhs ) );
	}
	
	/// Scalar-Vector multiplication.
	auto opBinaryRight( string op, Lhs )( ScalarClosure!Lhs lhs ) if( op == "*" && hasSameType!Lhs ) {
		return this * lhs;
	}
	
	// -- operators for non-closure arguments --
	/// Binary operations for nonclosure arguments. They get forwarded to their closure variants.
	auto opBinary( string op, int stride )( VectorView!(ElementType, false, stride) rhs ) if( op == "*" || op == "+" || op == "-" ) {
		return mixin( "this" ~ op ~ "vectorClosure( rhs )" );
	}
	
	/// ditto
	auto opBinaryRight( string op, int stride )( VectorView!(ElementType, true, stride) rhs ) if( op == "*"  || op == "+" || op == "-" ) {
		return mixin( "vectorClosure( rhs )" ~ op ~ "this" );
	}
	
	/// ditto
	auto opBinary( string op )( ElementType rhs ) if( op == "*" ) {
		return this * scalarClosure( rhs );	
	}
	
	/// ditto
	auto opBinaryRight( string op )( ElementType rhs ) if( op == "*" ) {
		return this * rhs;
	}
}

/// Struct representing vector scaling & transposition.
struct VectorUnary( VectorT_, bool Row_ ) {
	alias	VectorT_					OperandType;
	alias	BaseElementType!VectorT_	ElementType;
	alias	ScalarClosure!ElementType	ScalarT;
	alias	Row_						isRow;
	
	VectorT vector;
	ScalarT scalar;
	
	this( VectorT vector, ScalarT scalar ) {
		this.vector = vector;
		this.scalar = scalar;
	}
}
	

// -- free function operations --
/// Vector transposition.
VectorClosure!T transpose( T )( VectorClosure!T v ) {
	static if( v.expressionKind == ExpressionKind.VectorUnary )
		return vectorClosure( vectorUnary!( !v.isRow )( m.expr.lhs, m.expr.rhs ) );
	else
		return vectorClosure( vectorUnary!( v.isRow )( m.expr, v.ElementType(1) ) );
}

/// ditto
auto transpose( T, bool Row, int Stride )( VectorView!(T, Row, Stride) m ) {
	return transpose( vectorClosure( m ) );
}