import std.conv;
import scid.common.traits;

enum ExpressionKind {
	// Matrix result
	MatrixLiteral,
	MatrixUnary,
	MatrixProduct,
	MatrixSum,
	
	// Vector result
	VectorLiteral,
	MatrixVectorProduct,
	VectorMatrixProduct,
	VectorSum,
	VectorUnary,
		
	// Scalar result
	ScalarLiteral,
	VectorDot,
	VectorNorm,
	ScalarSum,
	ScalarProduct,
	ScalarDivision,
	ScalarNegation
}

template isExpression( T ) {
	enum bool isExpression = is( typeof(T.expressionKind) : ExpressionKind );
}

template isMatrixExpression( T ) {
	enum bool isMatrixExpression =
		T.expressionKind == ExpressionKind.MatrixUnary ||
		T.expressionKind == ExpressionKind.MatrixProduct ||
		T.expressionKind == ExpressionKind.MatrixSum;
}

template isVectorExpression( T ) {
	enum bool isVectorExpression =
		T.expressionKind == ExpressionKind.MatrixVectorProduct ||
		T.expressionKind == ExpressionKind.VectorMatrixProduct ||
		T.expressionKind == ExpressionKind.VectorSum ||
		T.expressionKind == ExpressionKind.VectorUnary;
}

template isScalarExpression( T ) {
	enum bool isScalarExpression =
		T.expressionKind == ExpressionKind.VectorDot ||
		T.expressionKind == ExpressionKind.VectorNorm ||
		T.expressionKind == ExpressionKind.ScalarSum ||
		T.expressionKind == ExpressionKind.ScalarProduct;
}

template expressionKindOf( T ) {
	static if( __traits( compiles, T.expressionKind ) )
		enum ExpressionKind expressionKindOf =	T.expressionKind;
	else //static if( is( T : MatrixView!E ) )
		enum ExpressionKind expressionKindOf = ExpressionKind.MatrixLiteral;
	//else static if( is( T : VectorView!E ) )
		//enum ExpressionKind expressionKindOf = ExpressionKind.VectorLiteral;
	//else static if( isNumeric!T || isComplex!T )
		//enum ExpressionKind expressionKindOf = ExpressionKind.ScalarLiteral;
}

mixin template UnaryExpression( ExpressionKind Ek_, OpT_ ) {
	alias OpT_					OperandType;
	alias BaseElementType!OpT_	ElementType;
	alias Ek_					expressionKind;
	
	OperandType	operand;
	
	this( OperandType operand ) {
		this.operand = operand;	
	}
}

struct BinaryExpression( ExpressionKind Ek_, LhsT_, RhsT_ ) {
	alias LhsT_						LhsType;
	alias RhsT_						RhsType;
	alias BaseElementType!LhsType	ElementType;
	alias Ek_						expressionKind;
	
	static assert( is(BaseElementType!RhsType : ElementType) );
	
	LhsType lhs;
	RhsType rhs;
	
	this( LhsType lhs, RhsType rhs ) {
		this.lhs = lhs;
		this.rhs = rhs;
	}
	
	string toString() {
		return to!string(expressionKind)~"(" ~ to!string(lhs) ~ ", " ~ to!string(rhs) ~ ")";
	}
}


mixin template LiteralExpression( ExpressionKind kind, alias closureFun ) {
	enum ExpressionKind expressionKind = kind;
	
	auto opBinary( string op, Rhs )( Rhs rhs ) if(
			expressionKindOf!Rhs == ExpressionKind.MatrixLiteral ||
			expressionKindOf!Rhs == ExpressionKind.VectorLiteral ||
			expressionKindOf!Rhs == ExpressionKind.ScalarLiteral ) {
		return mixin("closureFun(this)" ~ op ~ "rhs");
	}
	
	auto opUnary( string op )() {
		return closureFun(this).opUnary!( op )();
	}
}
