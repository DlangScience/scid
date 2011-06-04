/** Templates for compile-time introspection.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.common.traits;


import std.complex;
import std.traits;




/** Detect whether T is a complex floating-point type. */
template isComplex(T)
{
    enum bool isComplex = is(T==cfloat) || is(T==cdouble) || is(T==creal)
        || is(T==Complex!float) || is(T==Complex!double) || is(T==Complex!real);
}

version(unittest)
{
    static assert (isComplex!cdouble);
    static assert (isComplex!(Complex!double));
    static assert (!isComplex!double);
    static assert (!isComplex!int);
}


/** Detect whether T is a FORTRAN-compatible floating-point type.
    The FORTRAN-compatible types are: float, double, cfloat, cdouble.
*/
template isFortranType(T) if (!isComplex!T)
{
    enum bool isFortranType = is(T==float) || is(T==double);
}

template isFortranType(T) if (isComplex!T)
{
    enum bool isFortranType = isFortranType!(typeof(T.re));
}

version(unittest)
{
    static assert (isFortranType!double);
    static assert (isFortranType!cdouble);
    static assert (isFortranType!(Complex!double));
    static assert (!isFortranType!real);
    static assert (!isFortranType!creal);
    static assert (!isFortranType!(Complex!real));
}


/** Evaluates to true if T is a callable type, i.e. a function, delegate
    or functor type.
*/
template isCallable(T)
{
    enum bool isCallable = is (T == return)  ||  isFunctor!(T);
}


unittest
{
    struct FunctorT
    {
        void opCall(int i) { return; }
    }

    static assert (isCallable!(void function(int)));
    static assert (isCallable!(void delegate(int)));
    static assert (isCallable!(FunctorT));
}




/** Evaluates to true if T is a functor, i.e. a class or struct with
    an opCall method.
*/
template isFunctor(T)
{
    enum bool isFunctor = is (typeof(T.opCall) == return);
}


unittest
{
    struct FunctorT
    {
        void opCall(int i) { return; }
    }

    static assert (!isFunctor!(void function(int)));
    static assert (!isFunctor!(void delegate(int)));
    static assert (isFunctor!(FunctorT));
}




/** Evaluates to true if the following compiles:
    ---
    F f;
    ArgT x;
    RetT y = f(x);
    ---
*/
template isUnaryFunction(F, RetT, ArgT)
{
    enum isUnaryFunction = __traits(compiles,
    {
        F f;
        ArgT x;
        RetT y = f(x);
    });
}


unittest
{
    real f(real x) { return x*1.0; }
    static assert (isUnaryFunction!(typeof(&f), real, real));
    static assert (isUnaryFunction!(typeof(&f), double, int));
    static assert (!isUnaryFunction!(typeof(&f), int, real));
    static assert (!isUnaryFunction!(typeof(&f), real, string));
}




/** Evaluates to true if the following compiles:
    ---
    F f;
    ArgT x;
    f(x);       // no return type check
    ---
*/
template isUnaryFunction(F, ArgT)
{
    enum isUnaryFunction = __traits(compiles,
    {
        F f;
        ArgT x;
        f(x);
    });
}


unittest
{
    real f(real x) { return x*1.0; }
    static assert (isUnaryFunction!(typeof(&f), real));
    static assert (isUnaryFunction!(typeof(&f), double));
    static assert (isUnaryFunction!(typeof(&f), int));
    static assert (!isUnaryFunction!(typeof(&f), string));
}




/** Evaluates to true if FuncType is a vector field, i.e. a callable type
    that takes an ArgType[] array as input and returns a RetType[] array.
*/
template isVectorField(FuncType, ArgType, RetType = ArgType)
{
    enum bool isVectorField = isBufferVectorField!(FuncType, ArgType, RetType)
        || __traits(compiles, 
    {
        ArgType[] x;
        FuncType f;
        RetType[] y = f(x);
    });
}

version (unittest)
{
    static assert (isVectorField!(real[] delegate(real[]), real));
    static assert (isVectorField!(real[] function(real[]), real));
    static assert (!isVectorField!(real delegate(real[]), real));
    static assert (!isVectorField!(real[] delegate(real), real));
    static assert (!isVectorField!(real[] delegate(real[]), int));

    static assert (isVectorField!(int[] delegate(real[]), real, int));
    static assert (!isVectorField!(int[] delegate(real[]), int, int));
    static assert (!isVectorField!(int[] delegate(real[]), real, real));
}


/** Evaluates to true if FuncType is a vector field which takes an additional
    (but often optional) buffer of type RetType[] as the second argument.
*/
template isBufferVectorField(FuncType, ArgType, RetType = ArgType)
{
    enum bool isBufferVectorField = __traits(compiles,
    {
        ArgType[] x;
        RetType[] buf;
        FuncType f;
        RetType[] y = f(x, buf);
    });
}


version (unittest)
{
    static assert (isBufferVectorField!(real[] function(real[], real[]), real));
    static assert (!isBufferVectorField!(real[] function(real[], real), real));
    static assert (!isBufferVectorField!(real[] function(real[]), real));
}




/** Evaluates to true if T is an array type.


    If the second template parameter is specified it is also checked
    whether the elements of T are of type U.
*/
template isArray(T)
{
    static if (is (T U : U[]))
    {
        enum bool isArray = true;
    }
    else
    {
        enum bool isArray = false;
    }
}


/// ditto
template isArray(T, U)
{
    static if (is (T V : V[]))
    {
        enum bool isArray = is (V == U);
    }
    else
    {
        enum bool isArray = false;
    }
}


unittest
{
    static assert (isArray!(int[]));
    static assert (isArray!(int[3]));
    static assert (isArray!(int[][]));

    static assert (isArray!(int[], int));
    static assert (isArray!(int[3], int));
    static assert (isArray!(int[][], int[]));

    static assert (!isArray!(int));
    static assert (!isArray!(int[], long));
}




/** Evaluates to true if T is a one-dimensional array type.

    If the second template parameter is specified it is also checked
    whether the elements of T are of type U.
*/
template is1DArray(T)
{
    static if (is (T V : V[]))
    {
        enum bool is1DArray = !isArray!(V);
    }
    else
    {
        enum bool is1DArray = false;
    }
}


/// ditto
template is1DArray(T, U)
{
    static if (is (T V : V[]))
    {
        enum bool is1DArray = is (V == U);
    }
    else
    {
        enum bool is1DArray = false;
    }
}


unittest
{
    static assert(is1DArray!(double[]));
    static assert(is1DArray!(double[], double));
    static assert(!is1DArray!(double[], float));
    static assert(!is1DArray!(double));
    static assert(!is1DArray!(double[][]));
}




/** Evaluates to true if T is a two-dimensional array type.

    If the second template parameter is specified it is also checked
    whether the elements of T are of type U.
*/
template is2DArray(T)
{
    static if (is (T V : V[][]))
    {
        enum bool is2DArray = !isArray!(V);
    }
    else
    {
        enum bool is2DArray = false;
    }
}


/// ditto
template is2DArray(T, U)
{
    static if (is (T V : V[][]))
    {
        enum bool is2DArray = is (V == U);
    }
    else
    {
        enum bool is2DArray = false;
    }
}


unittest
{
    static assert(is2DArray!(double[][]));
    static assert(is2DArray!(double[][], double));
    static assert(!is2DArray!(double[][], float));
    static assert(!is2DArray!(double[]));
    static assert(!is2DArray!(double[][][]));
}





/** Evaluates to a type tuple of the argument types of T, where
    T must be either a function, delegate or functor type.
*/
template ArgumentTypeTuple(T)
{
    static if (is (T A == function))
    {
        alias A ArgumentTypeTuple;
    }
    else static if (is (T A == delegate))
    {
        alias ArgumentTypeTuple!(A) ArgumentTypeTuple;
    }
    else static if (is (T A == A*))
    {
        alias ArgumentTypeTuple!(A) ArgumentTypeTuple;
    }
    else static if (isFunctor!(T))
    {
        alias ArgumentTypeTuple!(typeof(T.opCall)) ArgumentTypeTuple;
    }
    else
    {
        static assert (false,
            "ArgumentTypeTuple: Not a callable type: "~T.stringof);
    }
}


private template UTTuple(T...) { alias T UTTuple; }
unittest
{
    struct FunctorT
    {
        real opCall(int i, bool b) { return 0.0; }
    }

    static assert (is (ArgumentTypeTuple!(real function(int, bool)) == UTTuple!(int, bool)));
    static assert (is (ArgumentTypeTuple!(real delegate(int, bool)) == UTTuple!(int, bool)));
    static assert (is (ArgumentTypeTuple!(FunctorT) == UTTuple!(int, bool)));
}




/** Evaluates to the inner element type of the array, vector, matrix
    or expression type T. If T is none of these, BaseElementType
	evaluates to T.
*/
template BaseElementType(T)
{
    static if (is (T E : E[]))
    {
        alias BaseElementType!(E) BaseElementType;
    }
    else static if (is (T E : E*))
    {
        alias BaseElementType!E BaseElementType;
    }
    else static if (is (typeof(&T.opIndex) E == return))
    {
        alias BaseElementType!E BaseElementType;
    }
    else static if ( __traits(compiles,T.ElementType.init) )
	{
		alias T.ElementType BaseElementType;	
	}
	else
    {
        alias T BaseElementType;
    }
}

unittest
{
    struct VectorS
    {
        size_t length;
        int opIndex(size_t i) { return 0; }
    }
    struct MatrixS
    {
        size_t rows,cols;
        int opIndex(size_t i, size_t j) { return 0; }
    }
	struct IntExpr
	{
		alias int ElementType;
	}

    static assert (is(BaseElementType!(int) == int));
    static assert (is(BaseElementType!(int[]) == int));
    static assert (is(BaseElementType!(int[][]) == int));

    static assert (is(BaseElementType!(int*) == int));
    static assert (is(BaseElementType!(int**) == int));
    static assert (is(BaseElementType!(int[]*) == int));

    static assert (is(BaseElementType!(VectorS) == int));
    static assert (is(BaseElementType!(MatrixS) == int));
	static assert (is(BaseElementType!(IntExpr) == int));
}

/** Tests if two types A and B have the same base element type. */
template haveSameElementType( A, B ) {
	enum haveSameElementType = is(BaseElementType!A : BaseElementType!B);
}

unittest {
	struct IntExpr
	{
		alias int ElementType;
	}
	
	struct IntExprExpr {
		alias BaseElementType!IntExpr ElementType;
	}
	
	static assert (haveSameElementType!(int**, int));
	static assert (haveSameElementType!(int[][], int*));
	
	static assert (haveSameElementType!(IntExpr, int));
	static assert (haveSameElementType!(IntExprExpr, IntExpr));
}




/** Checks whether all the types U... are implicitly
    convertible to T.
*/
template allConvertibleTo(T, U...) if (U.length > 0)
{
    static if (U.length == 1)
        enum allConvertibleTo = is(U[0] : T);
    else
        enum allConvertibleTo = is(U[0] : T) && allConvertibleTo!(T, U[1 .. $]);
}


unittest
{
    static assert (allConvertibleTo!(dchar, char, wchar, dchar));
    static assert (!allConvertibleTo!(wchar, char, wchar, dchar));
}
