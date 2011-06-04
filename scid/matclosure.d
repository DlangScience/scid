
import vecclosure;
import scalclosure;
import expressions;

import scid.common.traits;
import scid.matrix;

import std.conv;


/** Convinience function to create a MatrixClosure while
	avoiding typing the template parameter.
*/
MatrixClosure!T matrixClosure( T )( T value ) {
	return typeof( return )( value );
}

/** Convinience function to create a BinaryExpression
	representing matrix addition.
*/
BinaryExpression!( ExpressionKind.MatrixSum, L, R ) matrixSum( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );
}

/** Convinience function to create a BinaryExpression
	representing matrix multiplication.
*/
BinaryExpression!( ExpressionKind.MatrixProduct, L, R ) matrixProduct( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );	
}

/** Convinience function to create a BinaryExpression
	representing matrix-vector multiplication.
*/
BinaryExpression!( ExpressionKind.MatrixVectorProduct, L, R ) matrixVectorProduct( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );	
}

/** Convinience function to create a BinaryExpression
	representing vector-matrix multiplication.
*/
BinaryExpression!( ExpressionKind.VectorMatrixProduct, L, R ) vectorMatrixProduct( L, R )( L lhs, R rhs ) if( haveSameElementType!(L,R) ) {
	return typeof( return )( lhs, rhs );	
}

/** Convinience function to create a MatrixUnary
	expression.
*/

MatrixUnary!(T,ScalarClosure!(BaseElementType!T),tr,inv) matrixUnary( bool tr = false, bool inv = false, T )( T mat ) {
	return typeof( return )( mat, ScalarClosure!(BaseElementType!T)(1) );
}

/// ditto
MatrixUnary!(T,ScalarClosure!S,tr,inv) matrixUnary(bool tr = false, bool inv = false, T, S )( T mat, ScalarClosure!S scal ) {
	return typeof( return )( mat, scal );
}

/** This struct wraps all operations (expression trees) whose
	result is a matrix:
		MatrixSum
		MatrixProduct
		MatrixUnary
		MatrixView
*/
struct MatrixClosure( T ) {
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
	
	/// Construct a closure from an expression.
	this( T expr ) {
		this.expr = expr;
	}
	
	// -- operators for closure arguments --
	/// Matrix negation. Implemented as -1 * this
	auto opUnary( string op )() if( op == "-" ) {
		static if( expressionKind == ExpressionKind.MatrixUnary )
			return matrixClosure( matrixUnary!(expr.isTransposed,expr.isInverted)( expr.matrix, -expr.scalar ) );
		else
			return matrixClosure( matrixUnary( expr, ScalarClosure!ElementType(-1) ) );
	}
	
	/// Matrix addition.
	auto opBinary( string op, Rhs )( MatrixClosure!Rhs rhs ) if( op == "+" && hasSameType!Rhs ) {
		return matrixClosure( matrixSum( expr, rhs.expr ) );
	}
	
	/// Matrix subtraction. Implemented as addition with -1 * rhs.
	auto opBinary( string op, Rhs )( MatrixClosure!Rhs rhs ) if( op == "-" && hasSameType!Rhs ) {
		return this + (-rhs);
	}
	
	/// Matrix multiplication.
	auto opBinary( string op, Rhs )( MatrixClosure!Rhs rhs ) if( op == "*" && hasSameType!Rhs ) {
		return matrixClosure( matrixProduct( expr, rhs.expr ) );
	}
	
	/// Matrix-Vector multiplication.
	auto opBinary( string op, Rhs )( VectorClosure!Rhs rhs ) if( op == "*" && hasSameType!Rhs && !rhs.isRow ) {
		return vectorClosure( matrixVectorProduct( expr, rhs.expr ) );
	}
	
	/// Vector-Matrix multiplication.
	auto opBinaryRight( string op, Lhs )( VectorClosure!Lhs lhs ) if( op == "*" && hasSameType!Lhs && lhs.isRow ) {
		return vectorClosure( matrixVectorProduct( expr, transpose(lhs) ) );
	}
	
	/// Matrix-Scalar multiplication.
	auto opBinary( string op )( ScalarClosure!ElementType rhs ) if( op == "*"  ) {
		pragma(msg,BaseElementType!T);
		pragma(msg,typeof(rhs));
		static if( expressionKind == ExpressionKind.MatrixUnary )
			return matrixClosure( matrixUnary( expr.matrix, rhs * expr.scalar ) );
		else
			return matrixClosure( matrixUnary( expr, rhs ) );
	}
	
	/// Scalar-Matrix multiplication.
	auto opBinaryRight( string op, Lhs )( ScalarClosure!Lhs lhs ) if( op == "*" && hasSameType!Lhs ) {
		return this * lhs;
	}
	
	// -- operators for non-closure arguments --
	/// Binary operations for nonclosure arguments. They get forwarded to their closure variants.
	auto opBinary( string op, Storage stor, Triangle tri )( MatrixView!(ElementType,stor,tri) rhs ) if( op == "+" || op == "-" || op == "*" ) {
		return mixin( "this" ~ op ~ "matrixClosure( rhs )" );
	}
	
	/// ditto
	auto opBinary( string op, int stride )( VectorView!(ElementType, false, stride) rhs ) if( op == "*" ) {
		return mixin( "this" ~ op ~ "vectorClosure( rhs )" );
	}
	
	/// ditto
	auto opBinaryRight( string op, int stride )( VectorView!(ElementType, true, stride) rhs ) if( op == "*" ) {
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
	
	string toString() {
		return "MatrixClosure(" ~ to!string(expr) ~ ")";
	}
}

/// Struct represesnting any composition of scalings, transpositions and inversions
struct MatrixUnary( MatrixT_, ScalarT_, bool Tr_, bool Inv_ ) {
	alias	MatrixT_					OperandType;
	alias	BaseElementType!MatrixT_	ElementType;
	alias	ScalarT_					ScalarT;
	alias	Tr_							isTransposed;
	alias	Inv_						isInverted;
	
	enum	expressionKind = ExpressionKind.MatrixUnary;
	
	OperandType matrix;
	ScalarT     scalar;
	
	this( OperandType matrix, ScalarT scalar ) {
		this.matrix = matrix;
		this.scalar = scalar;
	}
	
	string toString() {
		static if( Tr_ && Inv_ )
			static string e = "^(-T)";
		else static if( Tr_ )
			static string e = "^T";
		else static if( Inv_ )
			static string e = "^(-1)";
		else
			static string e = "";
		return "MatrixUnary(" ~ to!string(matrix) ~ e ~ "*" ~ to!string(scalar) ~ ")";
	}
}

// -- free function operations --
/// Matrix transposition.
auto transpose( T )( MatrixClosure!T m ) {
	static if( m.expressionKind == ExpressionKind.MatrixUnary )
		return matrixClosure( matrixUnary!( !m.expr.isTransposed, m.expr.isInverted )( m.expr.matrix, m.expr.scalar ) );
	else
		return matrixClosure( matrixUnary!( true, false )( m.expr ) );
}

/// ditto
auto transpose( T, Storage stor, Triangle tri )( MatrixView!(T, stor, tri) m ) {
	return transpose( matrixClosure( m ) );
}

/// Matrix inversion.
auto inverse( T )( MatrixClosure!T m ) {
	static if( m.expressionKind == ExpressionKind.MatrixUnary )
		return matrixClosure( matrixUnary!( m.expr.isTransposed, !m.expr.isInverted )( m.expr.matrix, inverse(m.expr.scalar) ) );
	else
		return matrixClosure( matrixUnary!( false, true )( m.expr ) );
}

/// ditto
auto inverse( T, Storage stor, Triangle tri )( MatrixView!( T, stor, tri ) m ) {
	return inverse( matrixClosure( m ) );
}

