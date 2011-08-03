/** Contains types and other declarations for storing and reperesenting matrix and vector expressions.
    Solving the actual operations is performed in the scid.ops.eval.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ops.expression;

import scid.common.traits, scid.common.meta;
import std.traits, std.range;
import std.conv;
import std.typecons;

debug( Expression ) {
	import std.stdio;
}

/** An enumeration of all the operations supported. */
enum Operation {
	MatMatSum,    // Matrix = Matrix + Matrix
	MatMatProd,   // Matrix = Matrix * Matrix
	MatScalProd,  // Matrix = Matrix * Scalar
	MatInverse,   // Matrix = Matrix ^^ (-1)
	MatTrans,     // Matrix = Matrix.t
	
	RowRowSum,    // RowVector = RowVector + RowVector
	RowScalProd,  // RowVector = RowVector * Scalar
	RowMatProd,   // RowVector = RowVector * Matrix
	ColTrans,     // RowVector = ColVector.t
		
	ColColSum,    // ColumnVector = ColumnVector + ColumnVector
	ColScalProd,  // ColumnVector = ColumnVector * Scalar
	MatColProd,   // ColumnVector = Matrix * ColumnVector
	RowTrans,     // ColumnVector = RowVector.t
	
	DotProd,      // Scalar = RowVector * ColumnVector
	ScalScalSum,  // Scalar = Scalar +  Scalar
	ScalScalSub,  // Scalar = Scalar -  Scalar
	ScalScalProd, // Scalar = Scalar *  Scalar
	ScalScalDiv,  // Scalar = Scalar /	Scalar
	ScalScalPow,  // Scalar = Scalar ^^ Scalar
}

/** An enumeration of the closure types. */
enum Closure {
	Matrix,
	RowVector,
	ColumnVector,
	Scalar
}

/** Convinience function to create an expression object. */
Expression!( op, Unqual!Lhs, Unqual!Rhs ) expression( string op, Lhs, Rhs )( auto ref Lhs lhs, auto ref Rhs rhs ) {
	static assert( isConvertible!(BaseElementType!Lhs, BaseElementType!Rhs) &&
				   isConvertible!(BaseElementType!Lhs, BaseElementType!Rhs),
		"Invalid types in " ~ to!string(closureOf!Lhs) ~ " " ~ op ~ " " ~ to!string(closureOf!Rhs) ~ ": " ~
		BaseElementType!Lhs.stringof ~ " and " ~ BaseElementType!Rhs.stringof
	);
	return typeof( return )( lhs, rhs );
}

/// ditto	
Expression!( op, Unqual!Lhs ) expression( string op, Lhs )( auto ref Lhs lhs ) {
	return typeof( return )( lhs );
}

/** Struct representing a unary or binary expression. The parameters are the operator and the types involved.
    An expression tree provides a copy of all the operands involved thus creating a closure. The closure type
    refers to the type of the result of the expression. An expression provides the same operator overloads as
    a literal of the same closure type.
    TODO: Indexing & slicing are not implemented yet.
*/
struct Expression( string op_, Lhs_, Rhs_ = void ) {
	/** The type of the left hand side of the operation. */
	alias Lhs_ Lhs;
	
	/** The type of the right hand side of the operation. */
	alias Rhs_ Rhs;
	
	/** The operator used in the expression. */
	alias op_ operator;
	
	/** Whether the operation is binary. False if it's unary. */
	enum isBinary   = !is( Rhs : void );
	
	/** The element type of the result. */
	static if( isBinary )
		alias Promotion!(BaseElementType!Lhs_,BaseElementType!Rhs_) ElementType;
	else
		alias BaseElementType!Lhs_ ElementType;
	
	/** The closure type of the left hand side. */
	enum lhsClosure = closureOf!Lhs;
	
	static if( isBinary ) {
		/** The closure type of the right hand side. */
		enum rhsClosure = closureOf!Rhs;
		
		/** The type of the represented operation. */
		enum operation  = operationOf!(op_, lhsClosure, rhsClosure );
	} else {
		/** The type of the represented operation. */
		enum operation  = operationOf!(op_, lhsClosure );
	}
	
	/** The Operand mixin defines the proper operator overloads for a given closure type. */
	mixin Operand!( closureOf!operation );
	
	// Define methods appropriate for the closure type (e.g length, rows etc.)
	static if( isVectorClosure!closure ) {
		/** The type of the resulting vector. */
		// enum vctorType = closure == Closure.RowVector ? VectorType.Row : VectorType.Column;
		
		/** Length of the resulting vector. */
		size_t length() const @property {
			static if( operation == Operation.RowMatProd ) {	
				return rhs_.columns; 
			} else static if( operation == Operation.MatColProd ) {
				return lhs_.rows; 
			} else static if( isVectorClosure!lhsClosure ) {
				return lhs_.length;
			} else {
				return rhs_.length;
			}
		}
	} else static if( closure == Closure.Matrix ) {
		/** The number of dimensions of the resulting matrix. */
		size_t rows() const @property {
			static if( operation == Operation.MatTrans )
				return lhs_.columns;
			else
				return lhs_.rows;
		}
		
		/// ditto
		size_t columns() const @property {
			static if( operation == Operation.MatTrans )
				return lhs_.rows;
			else static if( operation == Operation.MatMatProd )
				return rhs_.columns;
			else
				return lhs_.columns;
		}
	}
	
	string toString() {
		static if( isScalar!Lhs )
			string lhsString = toImpl!(string,Lhs)(lhs);
		else
			string lhsString = lhs.toString;
		
		static if( isBinary ) {
			static if( isScalar!Rhs )
				string rhsString = toImpl!(string,Rhs)(rhs);
			else
				string rhsString = rhs.toString;
			return "( " ~ lhsString ~ " " ~ operator ~ " " ~ rhsString ~ " )";
		} else
			return lhsString ~ "." ~ operator;
	}
	
	/** The left hand side of the operation. */
	alias lhs_ lhs;
	Lhs   lhs_;
	
	static if( isBinary ) {
		/** The right hand side of the operation. */
		alias rhs_ rhs;
		Rhs   rhs_;
	}
}

/** Mixin that adds the appropriate operator overloads for any operand - literal or expression. */
template Operand( Closure closure_ ) {
	import scid.common.meta : Zero, One, MinusOne;
	
	/** The closure type. */
	alias closure_ closure;
	
	static if( closure == Closure.Scalar ) { 
		auto opBinaryRight( string op, E )( E newLhs )
				if( isConvertible!(E,ElementType) && 
				   (op == "^^" || op == "/" || op == "-" || op == "*") ) {
			return expression!op( newLhs, this );
		}
		
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs )
				if( op == "^^" || op == "/" || op == "-" || op == "+" ) {
			return expression!op( this, newRhs );
		}
	
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs )
				if( op == "*" ) {
			return expression!op( newRhs, this );
		}
	} else {
		@property auto ref t() {
			// this avoids an infinite compile time recursion
			static if( isExpression!(typeof(this)) && isTransposition!operation )
				auto r = this.lhs;
			else
				auto r = expression!"t"( this );
			return r;
		}
		
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "/" ) {
			return expression!"*"( this, expression!"/"(One!ElementType, newRhs) );
		}
	
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "*" || op == "+" ) {
			return expression!op( this, newRhs );
		}
		
		auto opBinaryRight( string op, E )( E newLhs ) if( isConvertible!(E,ElementType) && op == "*" ) {
			return expression!op( this, newLhs );
		}
		
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "-" ) {
			return expression!"+"( this, expression!"*"(newRhs, MinusOne!ElementType) );
		}
	}
	
	auto opUnary( string op )() if( op == "-" ) {
		return expression!"*"( this, MinusOne!ElementType );
	}
}

/** Get the type of the result of an expression. */
template ExpressionResult( E ) {
	static if( isLeafExpression!E ) {
		// if we've reached a leaf return it
		alias E ExpressionResult;
	} else static if( E.closure == Closure.Scalar ) {
		// if the node results in a scalar then the result is of the same type as the element type
		alias E.ElementType ExpressionResult;
	} else static if( isTransposition!(E.operation) ) {
		// if the node is a transposition then the result is the Transposed of the child node
		alias ExpressionResult!(E.Lhs).Transposed ExpressionResult;
	} else static if( E.isBinary ) {
		alias Promotion!(ExpressionResult!(E.Lhs),ExpressionResult!(E.Rhs)) ExpressionResult;
	} else {
		alias ExpressionResult!(E.Lhs) ExpressionResult;
	}
}

/** Operation trait. */
template isScaling( Operation op ) {
	enum isScaling =
		op == Operation.MatScalProd ||
		op == Operation.RowScalProd ||
		op == Operation.ColScalProd ||
		op == Operation.ScalScalProd;
}

/// ditto
template isInnerProduct( Operation op ) {
	enum isInnerProduct = 
		op == Operation.MatMatProd ||
		op == Operation.RowMatProd ||
		op == Operation.MatColProd;
}

/// ditto
template isAddition( Operation op ) {
	enum isAddition =
		op == Operation.MatMatSum ||
		op == Operation.RowRowSum ||
		op == Operation.ColColSum;
}

/// ditto
template isTransposition( Operation op ) {
	enum isTransposition =
		op == Operation.RowTrans ||
		op == Operation.ColTrans ||
		op == Operation.MatTrans;	
}

/** Get the closure type of an expression or a literal. */
template closureOf( T ) {
	static if( isScalar!(Unqual!T) ) { 
		enum closureOf = Closure.Scalar;
	} else static if( is( typeof(T.closure) : Closure ) ) {
		enum closureOf = T.closure;
	} else {
		static assert( false, T.stringof ~ " is not a valid expression." );
	}
}

/** Closure trait. */
template isVectorClosure( Closure c ) {
	enum isVectorClosure =
		c == Closure.RowVector ||
		c == Closure.ColumnVector;
}

/** Get the closure type of an operation. */
template closureOf( Operation op ) {
	enum closureOf = operationClosures[ op ];
}

/** Find the more general type of two. */
template Promotion( A, B ) {
	alias PromotionImpl!(A,B).Result Promotion;	
}

template isRefCounted( T ) {
	static if( is( T E : RefCounted!(E, autoInit), uint autoInit ) )
		enum isRefCounted = true;
	else
		enum isRefCounted = false;
}

template refCountedAutoInitializeOf( T ) {
	static if(	is( T E : RefCounted!(E, autoInit), uint autoInit ) )
		enum refCountedAutoInitializeOf = autoInit;
}

private template PromotionImpl( A, B ) {
	static if( isScalar!A && isScalar!B ) {
		alias typeof( A.init + B.init ) Result;
	} else static if( is( A E : RefCounted!(E, autoInit), uint autoInit ) ) {
		alias Promotion!(E,B) Result;
	} else static if( is( B E : RefCounted!(E, autoInit), uint autoInit ) )  {
		alias Promotion!(A,E) Result;
	} else static if( is(A.Promote!B R) ) {
		alias R Result;
	} else static if( is(B.Promote!A R) ) {
		alias R Result;
	} else {
		static assert( false, "Types '" ~ A.stringof ~ "' and '" ~ B.stringof ~ "' do not define a promotion. " );
	}
}


// This gives a small compile time boost when checking the closure of a given operation.
private enum operationClosures = [
	Closure.Matrix,           // MatMatSum,
	Closure.Matrix,           // MatMatProd,
	Closure.Matrix,			  // MatScalProd,
	Closure.Matrix,			  // MatInverse,
	Closure.Matrix,			  // MatTrans,
							  
	Closure.RowVector,		  // RowRowSum,
	Closure.RowVector,		  // RowScalProd,
	Closure.RowVector,		  // RowMatProd,
	Closure.RowVector,		  // ColTrans,
							  	 	
	Closure.ColumnVector,	  // ColColSum,
	Closure.ColumnVector,	  // ColScalProd,
	Closure.ColumnVector,	  // MatColProd,
	Closure.ColumnVector,	  // RowTrans,
							  
	Closure.Scalar,			  // DotProd,
	Closure.Scalar,			  // ScalScalSum,
	Closure.Scalar,			  // ScalScalSub,
	Closure.Scalar,           // ScalScalProd,
	Closure.Scalar,           // ScalScalDiv,
	Closure.Scalar            // ScalScalPow,
];


// Find the operation type given the operator and the closure types of the operands
private template operationOf( string op, Closure l, Closure r ) if( op == "*" ) {
	
	     static if( l == Closure.Matrix       && r == Closure.Matrix )          { enum operationOf = Operation.MatMatProd;   }
	else static if( l == Closure.Matrix       && r == Closure.Scalar )          { enum operationOf = Operation.MatScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.Matrix )          { enum operationOf = Operation.MatScalProd;  }
	else static if( l == Closure.Matrix       && r == Closure.ColumnVector )    { enum operationOf = Operation.MatColProd;   }
	else static if( l == Closure.RowVector    && r == Closure.Matrix )          { enum operationOf = Operation.RowMatProd;   }
	else static if( l == Closure.RowVector    && r == Closure.ColumnVector )    { enum operationOf = Operation.DotProd;      }
	
	else static if( l == Closure.RowVector    && r == Closure.Scalar )          { enum operationOf = Operation.RowScalProd;  }
	else static if( l == Closure.ColumnVector && r == Closure.Scalar )          { enum operationOf = Operation.ColScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.ColumnVector )    { enum operationOf = Operation.RowScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.RowVector )       { enum operationOf = Operation.ColScalProd;  }
	
	else static if( l == Closure.Scalar       && r == Closure.Scalar )          { enum operationOf = Operation.ScalScalProd; }
	
	else static if( l == Closure.ColumnVector && r == Closure.RowVector )       { static assert( false, "ColumnVector * RowVector is not supported yet." );  }
	
	else static assert( false, "Invalid multiplication between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l, Closure r ) if( op == "+" ) {
	     static if( l == Closure.Matrix       && r == Closure.Matrix )          { enum operationOf = Operation.MatMatSum;    }
	else static if( l == Closure.RowVector    && r == Closure.RowVector )       { enum operationOf = Operation.RowRowSum;    }
	else static if( l == Closure.ColumnVector && r == Closure.ColumnVector )    { enum operationOf = Operation.ColColSum;    }
	else static if( l == Closure.Scalar       && r == Closure.Scalar )          { enum operationOf = Operation.ScalScalSum;  }
	
	else static assert( false, "Invalid addition between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l ) if( op == "t" ) {
	     static if( l == Closure.Matrix       )  { enum operationOf = Operation.MatTrans;  }
	else static if( l == Closure.RowVector    )  { enum operationOf = Operation.RowTrans;  }
	else static if( l == Closure.ColumnVector )  { enum operationOf = Operation.ColTrans;  }
	
	else static assert( false, "Invalid transposition of " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l, Closure r ) if( op == "/" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalDiv;  }
	else static assert( false, "Invalid division between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l, Closure r ) if( op == "^^" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalPow;  }
	else static assert( false, "Invalid pow between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l, Closure r ) if( op == "-" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalSub;  }
	else static assert( false, "Invalid subtraction between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

private template operationOf( string op, Closure l ) if( op == "inv" ) {
	static if( l == Closure.Matrix )  { enum operationOf = Operation.MatInverse;  }
	else static assert( false, "Invalid inversion of " ~ to!string(l) );
}