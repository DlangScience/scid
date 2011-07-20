module scid.internal.expression;

import scid.common.traits, scid.common.meta;
import std.traits, std.range;

debug = Expression;
// version = ExpressionsUsePointers;

debug( Expression ) {
	import std.conv;
}

enum Operation {
	MatMatSum,
	MatMatProd,
	MatScalProd,
	MatColProd,
	RowMatProd,
	MatInverse,
	MatTrans,
	
	RowRowSum,
	RowScalProd,
	ColTrans,
		
	ColColSum,
	ColScalProd,
	RowTrans,
	
	DotProd,
	ScalScalSum,
	ScalScalSub,
	ScalScalProd,
	ScalScalDiv,
	ScalScalPow,
}

enum Closure {
	Matrix,
	RowVector,
	ColumnVector,
	Scalar
}

template isExpression( T ) {
	static if( is( typeof(T.operation) : Operation ) )
		enum isExpression = true;
	else
		enum isExpression = false;
}

template isLiteral( T ) {
	static if( is( T E : E* ) ) {
		enum isLiteral = isLiteral!E;
	} else static if( is(typeof(T.operation)) ) {
		enum isLiteral = false;
	} else {
		enum isLiteral = true;
	}
}

template ReferencedBy( T ) {
	static if( is( T.Referenced ) )
		alias T.Referenced ReferencedBy;
	else
		alias T ReferencedBy;
}

version( ExpressionsUsePointers ) {
	Expression!( op, Unqual!Lhs, Unqual!Rhs ) expression( string op, Lhs, Rhs )( auto ref Lhs lhs, auto ref Rhs rhs )  {
		return typeof( return )( toRef(lhs), toRef(rhs) );
	}

	Expression!( op, Unqual!Lhs ) expression( string op, Lhs )( auto ref Lhs lhs )  {
		return typeof( return )( toRef(lhs) );
	}

	template Operand( T ) {
		alias Select!( is(T == struct), T*, T ) Operand;
	}

	auto toRef( T )( auto ref T v ) {
		static if( is( T == struct ) ) {
			return &v;
		} else {
			return v;
		}
	}

	ref T unref( T )( T* v ) { return *v; }
	T unref( T )( T v ) if( !is( T E : E* ) ) { return v; }
} else {
	Expression!( op, Unqual!Lhs, Unqual!Rhs ) expression( string op, Lhs, Rhs )( auto ref Lhs lhs, auto ref Rhs rhs )  {
		return typeof( return )( lhs, rhs );
	}

	Expression!( op, Unqual!Lhs ) expression( string op, Lhs )( auto ref Lhs lhs )  {
		return typeof( return )( lhs );
	}
}

struct Expression( string op_, Lhs_, Rhs_ = void ) {
	alias BaseElementType!Lhs_ ElementType;
	alias Lhs_                 Lhs;
	alias Rhs_                 Rhs;
	alias op_                  operator;
	
	template isPrimitiveClosure( T ) {
		enum isPrimitiveClosure = is( T : ElementType[] ) || is( T : ElementType[][] ) || is( T : ElementType );
	}
	
	enum isBinary   = (!is( Rhs_ : void ));
	enum lhsClosure = closureOf!Lhs;
	
	static if( isBinary ) {
		enum rhsClosure = closureOf!Rhs;
		enum operation  = operationOf!(op_, lhsClosure, rhsClosure );
		this()( auto ref Lhs lhs, auto ref Rhs rhs ) {
			this.lhs = lhs; this.rhs = rhs;
		}
	} else {
		enum operation  = operationOf!(op_, lhsClosure );
		this()( auto ref Lhs lhs ) {
			this.lhs = lhs;
		}
	}
	
	enum closure = operationClosures[ operation ];
	
	auto opUnary( string op )() if( op == "-" ) {
		return this * MinusOne!ElementType;
	}
	
	auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs )  if( op == "-" ) {
		return expression!"+"( this, newRhs * MinusOne!ElementType );
	}
	
	auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs )  if( op == "+" ) {
		return expression!"+"( this, newRhs );
	}
	
	static if( closure == Closure.Scalar ) {
		auto opBinaryRight( string op, S )( S newLhs )
				if( (op == "^^" || op == "/" || op == "-" || op == "*") && isPrimitiveClosure!S ) {
			return expression!op( newLhs, this );
		}
		
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "*" ) {
			return newRhs * this;
		}
		
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "^^" || op == "/" || op == "-" ) {
			return expression!op( this, newRhs );
		}
	} else {
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs )  if( op == "*" ) {
			static if( is( NewRhs : ElementType ) && isScalingOperation!operation ) {
				static if( is( Rhs : ElementType ) ) {
					return typeof(this)( lhs_, rhs_*newRhs );
				} else {
					return expression!"*"( this, newRhs );
				}
			} else
				return expression!"*"( this, newRhs  );
		}
	
		auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "/" ) {
			return this * ( One!ElementType / newRhs );
		}
		
		auto opBinaryRight( string op, S )( S newLhs ) if( op == "*" && isPrimitiveClosure!S ) {
			return mixin( "this" ~ op ~ "newLhs" );
		}
	}
	
	static if( closure != Closure.Scalar ) {
		@property auto t() {
			static if( op_ == "t" )
				return lhs_;
			else static if( isScalingOperation!operation ) {
				return expression!"*"( lhs_.t, rhs_ );
			} else static if( operation == Operation.DotProd ) {
				return Expression!("*", typeof(rhs_.t), typeof(lhs_.t))( rhs_.t, lhs_.t );
			} else
				return expression!"t"( this );
		}	
	}
	
	static if( closure == Closure.RowVector || closure == Closure.ColumnVector ) {
		@property size_t length() const { return lhs_.length; }
	} else static if( closure == Closure.Matrix ) {
		static if( operation == Operation.MatTrans ) {
			@property size_t rows()    const { return lhs_.columns; }
			@property size_t columns() const { return lhs_.rows; }
		} else static if( operation == Operation.MatMatProd ) {
			@property size_t rows()    const { return lhs_.rows; }
			@property size_t columns() const { return rhs_.columns; }
		} else {
			@property size_t rows()    const { return lhs_.rows; }
			@property size_t columns() const { return lhs_.columns; }
		}
	}
	
	version( ExpressionsUsePointers ) {
		@property ref auto lhs() { return unref(lhs_); }
		
		Operand!Lhs lhs_;
		static if( isBinary ) {
			@property ref auto rhs() { return unref(rhs_); }
			Operand!Rhs rhs_;
		}
	} else {
		// @property ref Lhs lhs() { return lhs_; }
		
		//string toString() const {
		//	static if( isBinary )
		//		return to!string(operation) ~ "(" ~ to!string(lhs_) ~ ", " ~ to!string(rhs_) ~ ")";
		//	else						    
		//		return to!string(operation) ~ "(" ~ to!string(lhs_) ~ ")";
		//}
		
		Lhs lhs;
		static if( isBinary ) {
			// @property ref Rhs rhs() { return rhs_; }
			Rhs rhs;
			alias rhs rhs_;
		}
		
		alias lhs lhs_;
		
	}
}

template ExpressionResult( E ) {
	static if( isLiteral!E )
		alias ReferencedBy!E ExpressionResult;
	else static if( E.closure == Closure.Scalar )
		alias E.ElementType ExpressionResult;
	else static if( E.operation == Operation.RowTrans || E.operation == Operation.ColTrans ) {
			alias ExpressionResult!(E.Lhs).Transposed ExpressionResult;
	} else static if( E.operation == Operation.RowMatProd || E.operation == Operation.MatColProd ) {
		static if( E.lhsClosure == Closure.Matrix )
			alias ExpressionResult!(E.Rhs) ExpressionResult;
		else
			alias ExpressionResult!(E.Lhs) ExpressionResult;
	} else static if( E.operation == Operation.RowRowSum || E.operation == Operation.ColColSum || E.operation == Operation.RowScalProd || E.operation == Operation.ColScalProd ) {
		alias ExpressionResult!(E.Lhs) ExpressionResult;
	} else static if( E.operation == Operation.MatMatProd ) {
		alias ExpressionResult!(E.Lhs) ExpressionResult;
	} else static assert( false, E.stringof ~ " hasn't got a result." );
}

@property auto t( T )( T[] array )
		if( closureOf!T == Closure.Scalar ||
			closureOf!T == RowVector ||
		    closureOf!T == ColumnVector
		) {
	return expression!"t"( array );
}

template Literal( Closure closure_ ) {
	alias closure_ closure;
	
	import scid.common.meta;
	
	static if( closure != Closure.Scalar ) {
		@property auto t()  {
			return expression!"t"( this );
		}
	}
	
	public auto opUnary( string op )() if( op == "-" ) {
		return this * MinusOne!ElementType;
	}
	
	public auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "-" ) {
		return expression!"+"( this, newRhs * MinusOne!ElementType );
	}
	
	public auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "+" || op == "*" ) {
		return expression!op( this, newRhs);
	}
	
	auto opBinary( string op, NewRhs )( auto ref NewRhs newRhs ) if( op == "/" ) {
		return this * ( One!ElementType / newRhs );
	}
	
	
	public auto opBinaryRight( string op )( ElementType newLhs ) if( op == "*" || (closure==Closure.Scalar && op == "+") ) {
		return mixin( "this" ~ op ~ "newLhs" );
	}
}

template isScalingOperation( Operation op ) {
	enum isScalingOperation = ( op == Operation.MatScalProd || op == Operation.RowScalProd || op == Operation.ColScalProd || op == Operation.ScalScalProd );	
}

template isInnerProduct( Operation op ) {
	enum isInnerProduct =  op == Operation.MatMatProd || op == Operation.RowMatProd || op == Operation.MatColProd;
}

template isAddition( Operation op ) {
	enum isAddition = op == Operation.MatMatSum || op == Operation.RowRowSum || op == Operation.ColColSum;
}

template isTransposition( Operation op ) {
	enum isTransposition = op == Operation.RowTrans || op == Operation.ColTrans || op == Operation.MatTrans;	
}

enum operationClosures = [
	Closure.Matrix,           // MatMatSum,
	Closure.Matrix,           // MatMatProd,
	Closure.Matrix,			  // MatScalProd,
	Closure.ColumnVector,	  // MatColProd,
	Closure.RowVector,		  // RowMatProd,
	Closure.Matrix,			  // MatInverse,
	Closure.Matrix,			  // MatTrans,
							  
	Closure.RowVector,		  // RowRowSum,
	Closure.RowVector,		  // RowScalProd,
	Closure.RowVector,		  // ColTrans,
							  	 	
	Closure.ColumnVector,	  // ColColSum,
	Closure.ColumnVector,	  // ColScalProd,
	Closure.ColumnVector,	  // RowTrans,
							  
	Closure.Scalar,			  // DotProd,
	Closure.Scalar,			  // ScalScalSum,
	Closure.Scalar,			  // ScalScalSub,
	Closure.Scalar,           // ScalScalProd,
	Closure.Scalar,           // ScalScalDiv,
	Closure.Scalar            // ScalScalPow,
];




template operationOf( string op, Closure l, Closure r ) if( op == "*" ) {
	
	     static if( l == Closure.Matrix       && r == Closure.Matrix )          { enum operationOf = Operation.MatMatProd;   }
	else static if( l == Closure.Matrix       && r == Closure.Scalar )          { enum operationOf = Operation.MatScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.Matrix )          { enum operationOf = Operation.MatScalProd;  }
	else static if( l == Closure.Matrix       && r == Closure.ColumnVector )    { enum operationOf = Operation.MatColProd;   }
	else static if( l == Closure.RowVector    && r == Closure.Matrix )          { enum operationOf = Operation.RowMatProd;   }
	else static if( l == Closure.RowVector    && r == Closure.ColumnVector )    { enum operationOf = Operation.DotProd;      }
	else static if( l == Closure.ColumnVector && r == Closure.RowVector )       { enum operationOf = Operation.DotProd;      }
	
	else static if( l == Closure.RowVector    && r == Closure.Scalar )          { enum operationOf = Operation.RowScalProd;  }
	else static if( l == Closure.ColumnVector && r == Closure.Scalar )          { enum operationOf = Operation.ColScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.ColumnVector )    { enum operationOf = Operation.RowScalProd;  }
	else static if( l == Closure.Scalar       && r == Closure.RowVector )       { enum operationOf = Operation.ColScalProd;  }
	
	else static if( l == Closure.Scalar       && r == Closure.Scalar )          { enum operationOf = Operation.ScalScalProd; }
	
	else static assert( false, "Invalid multiplication between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l, Closure r ) if( op == "+" ) {
	     static if( l == Closure.Matrix       && r == Closure.Matrix )          { enum operationOf = Operation.MatMatSum;    }
	else static if( l == Closure.RowVector    && r == Closure.RowVector )       { enum operationOf = Operation.RowRowSum;    }
	else static if( l == Closure.ColumnVector && r == Closure.ColumnVector )    { enum operationOf = Operation.ColColSum;    }
	else static if( l == Closure.Scalar       && r == Closure.Scalar )          { enum operationOf = Operation.ScalScalSum;  }
	
	else static assert( false, "Invalid addition between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l ) if( op == "t" ) {
	     static if( l == Closure.Matrix       )  { enum operationOf = Operation.MatTrans;  }
	else static if( l == Closure.RowVector    )  { enum operationOf = Operation.RowTrans;  }
	else static if( l == Closure.ColumnVector )  { enum operationOf = Operation.ColTrans;  }
	
	else static assert( false, "Invalid transposition of " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l, Closure r ) if( op == "/" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalDiv;  }
	else static assert( false, "Invalid division between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l, Closure r ) if( op == "^^" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalPow;  }
	else static assert( false, "Invalid pow between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l, Closure r ) if( op == "-" ) {
	static if( l == Closure.Scalar && r == Closure.Scalar ) { enum operationOf = Operation.ScalScalSub;  }
	else static assert( false, "Invalid subtraction between " ~ to!string(l) ~ " and " ~ to!string(r) );
}

template operationOf( string op, Closure l ) if( op == "inv" ) {
	     static if( l == Closure.Matrix       )  { enum operationOf = Operation.MatInverse;  }
	
	else static assert( false, "Invalid inversion of " ~ to!string(l) ~ " and " ~ to!string(r) );
}


template closureOf( T ) {
	static if( isFortranType!(Unqual!T) ) { 
		enum closureOf = Closure.Scalar;
	} else static if( is( typeof(T.closure) : Closure ) ) {
		enum closureOf = T.closure;
	} else static if( is( T E : E[] ) && isFortranType!(Unqual!E) ) {
		enum closureOf = Closure.ColumnVector;
	} else static if( is( T E : E[][] ) && isFortranType!(Unqual!E) ) {
		enum closureOf = Closure.Matrix;
	} else static if( isInputRange!T && hasLength!T ) {
		enum c_ = closureOf!(typeof(T.front));
		static if( c_ == Closure.Scalar )
			enum closureOf = Closure.ColumnVector;
		else static if( c_ == Closure.ColumnVector || c_ == Closure.RowVector )
			enum closureOf = Closure.Matrix;
		else
			static assert( false, T.stringof ~ " invalid range for expression." );
	} else {
		static assert( false, T.stringof ~ " is not a valid expression." );
	}
}

