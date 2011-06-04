import expressions;

import scid.common.traits;

import std.conv;

/** Convinience function to create a ScalarClosure while
	avoiding typing the template parameter.
*/
ScalarClosure!T scalarClosure( T )( T value ) {
	return typeof( return )( value );
}

/** Convinience function to create a BinaryExpression
	representing scalar addition.
*/
BinaryExpression!( ExpressionKind.ScalarSum, L, R ) scalarSum( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}

/** Convinience function to create a BinaryExpression
	representing scalar multiplication.
*/
BinaryExpression!( ExpressionKind.ScalarProduct, L, R ) scalarProduct( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}

/** Convinience function to create a BinaryExpression
	representing scalar division.
*/
BinaryExpression!( ExpressionKind.ScalarDivision, L, R ) scalarDivision( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}

/** Convinience function to create a BinaryExpression
	representing scalar addition.
*/
BinaryExpression!( ExpressionKind.ScalarSum, L, R ) scalarSum( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}


struct ScalarClosure( T ) {
	/// Convinence check for same base element type (all the operators require it).
	template hasSameType( U ) {
		enum hasSameType = haveSameElementType!(T,U);
	}
	
	/// Type of the enclosed expression.
	alias T						EnclosedType;
	alias BaseElementType!T		ElementType;
	alias expressionKindOf!T	expressionKind;
	
	/// The enclosed expression.
	T expr;
	
	this( T expr ) {
		this.expr = expr;
	}
	
	this( this ) {
		this.expr = expr;
	}
	
	// -- operators for closure arguments --
	/// Scalar negation.
	auto opUnary( string op )() if( op == "-" ) {
		return scalarClosure( scalarNegation( expr ) );
	}
	
	/// Scalar addition.
	auto opBinary( string op, Rhs )( ScalarClosure!Rhs rhs ) if( op == "+" && hasSameType!Rhs ) {
		return scalarClosure( scalarSum( expr, rhs.expr ) );
	}
	
	/// Scalar subtraction. Implemented as addition with -rhs.
	auto opBinary( string op, Rhs )( ScalarClosure!Rhs rhs ) if( op == "-" && hasSameType!Rhs ) {
		return this + (-rhs);
	}
	
	/// Scalar multiplication.
	auto opBinary( string op, Rhs )( ScalarClosure!Rhs rhs ) if( op == "*" && hasSameType!Rhs ) {
		return scalarClosure( scalarProduct( expr, rhs.expr ) );	
	}
	
	/// Scalar division.
	auto opBinary( string op, Rhs )( ScalarClosure!Rhs rhs ) if( op == "/" && hasSameType!Rhs ) {
		return scalarClosure( scalarDivision( expr, rhs.expr ) );	
	}
	
	// -- operators for non-closure arguments --
	/// Binary operations for nonclosure arguments. They get forwarded to their closure variants.
	auto opBinary( string op )( ElementType rhs ) {
		return mixin("this" ~ op ~ "scalarClosure( rhs )");	
	}
	
	/// ditto
	auto opBinaryRight( string op )( ElementType rhs ){
		return mixin("scalarClosure( rhs )" ~ op ~ "this");
	}
	
	string toString() {
		return "ScalarClosure(" ~ to!string(expr) ~ ")";
	}
}