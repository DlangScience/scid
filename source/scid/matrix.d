/**
This module contains the MatrixView class as well as some
functions related to it.

Authors:    Lars Tandle Kyllingstad
Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
License:    Boost License 1.0
*/
module scid.matrix;

import std.array;
import std.string: format;
import std.traits;

import scid.core.meta;
import scid.core.traits;


/// Various matrix representations.
enum Storage
{
    General,        /// General (dense) matrices
    Triangular,     /// Packed storage of triangular matrices
    Symmetric,      /// Packed storage of symmetric matrices
}


/**
In packed storage (triangular, symmetric, and Hermitian matrices),
one can choose to store either the upper or lower triangle.
*/
enum Triangle : char
{
    Upper = 'U',    /// Store upper triangle
    Lower = 'L'     /// Store lower triangle
}


/**
A convenience function that allocates memory for a matrix (using the
GC), optionally sets the values of the matrix elements, and returns
a MatrixView of the allocated memory.

Examples:
---
// Allocate general dense 3x4 matrix:
auto denseMatrix = matrix!real(3, 4);

// Allocate dense 3x2 zero-filled matrix:
auto denseZeroMatrix = matrix!real(3, 2, 0.0L);

// Allocate lower triangular 3x3 matrix:
auto loMatrix = matrix!(real, Storage.Triangular, Triangle.Lower)(3);

// Allocate upper triangular 2x2 matrix where the upper
// triangular elements are set to 3.14.
auto upMatrix = matrix!(real, Storage.Triangular)(2, 3.14L);
---
*/
MatrixView!(T) matrix (T) (size_t rows, size_t cols) pure
{
    return typeof(return)(new T[rows*cols], rows, cols);
}

/// ditto
MatrixView!(T) matrix(T) (size_t rows, size_t cols, T init) pure
{
    auto array = new T[rows*cols];
    array[] = init;
    return typeof(return)(array, rows, cols);
}

/// ditto
MatrixView!(T, stor, tri) matrix
    (T, Storage stor, Triangle tri = Triangle.Upper)
    (size_t n, T init=T.init)
    pure
    if (stor == Storage.Triangular)
{
    auto array = new T[(n*n+n)/2];
    if (init != T.init) array[] = init; // Because of DMD bug #3576 this can't
                                        // be done with a function overload.
    return typeof(return)(array, n);
}

unittest
{
    import std.math;
    auto dense1 = matrix!real(3, 4);
    assert (dense1.rows == 3  &&  dense1.cols == 4);
    assert (isNaN(dense1[1,2]));

    auto dense2 = matrix(4, 3, 1.0);
    assert (dense2.rows == 4  &&  dense2.cols == 3);
    assert (dense2[1,2] == 1.0);

    auto upTri = matrix!(real, Storage.Triangular)(3);
    assert (upTri.rows == 3  &&  upTri.cols == 3);
    assert (isNaN(upTri[0,2])  &&  upTri[2,0] == 0);

    auto loTri = matrix!(double, Storage.Triangular, Triangle.Lower)(3, 1.0);
    assert (loTri.rows == 3  &&  loTri.cols == 3);
    assert (loTri[0,2] == 0  &&  loTri[2,0] == 1.0);
}


/**
A convenience function that creates a copy of the input matrix. Memory for
the copy is allocated using the GC.
*/
MatrixView!(T, stor, tri) copy(T, Storage stor, Triangle tri)
    (const MatrixView!(T, stor, tri) m)
    pure
{
    MatrixView!(T, stor, tri) mcopy;
    mcopy.rows = m.rows;
    mcopy.cols = m.cols;
    mcopy.array = m.array.dup;
    return mcopy;
}

unittest
{
    auto a = matrix!double(2, 2);
    a[1,0] = 1.0;
    auto b = copy(a);
    b[1,0] = 2.0;
    assert (b[1,0] == 2.0  &&  a[1,0] == 1.0);
}


/**
A matrix-like view of the contents of an array.

In order to be compatible with LAPACK routines, the following matrix
representations (i.e., memory layouts) are supported.

General_matrices:
The elements of dense matrices are stored in column-major order.
This means that if the wrapped array contains the elements
---
a b c d e f g h i j k l
---
then a 3x4 dense matrix view of this array looks like this:
---
a d g j
b e h k
c f i l
---

Triangular_matrices:
Triangular matrices are required to be square. If the wrapped
array contains the six elements
---
a b c d e f
---
then the resulting 3x3 upper and lower triangular matrix views
will look like this, respectively:
---
a b d         a 0 0
0 c e   and   b d 0
0 0 f         c e f
---

Symmetric_matrices:
Symmetric matrices are stored in the same way as triangular
matrices. This means that for the array above, the corresponding
symmetric matrix view will be
---
a b d       a b c
b c e   or  b d e
d e f       c e f
---
depending on whether the upper or lower triangle is stored.

Hermitian_matrices:
Hermitian matrices are not implemented yet.

See_also:
LAPACK User's Guide: Matrix storage schemes,
$(LINK http://www.netlib.org/lapack/lug/node121.html)
*/
struct MatrixView (T, Storage stor = Storage.General,
    Triangle tri = Triangle.Upper)
{
    enum Storage storage = stor;
    enum Triangle triangle = tri;

private:
    // Writing fully-qualified names in static ifs gets tiresome, so
    // we introduce a few flags.
    enum : bool
    {
        isGen   = (storage == Storage.General),
        isTri   = (storage == Storage.Triangular),
        isUpTri = (isTri && triangle == Triangle.Upper),
        isLoTri = (isTri && triangle == Triangle.Lower),
        isSym   = (storage == Storage.Symmetric),
        isUpSym = (isSym && triangle == Triangle.Upper),
        isLoSym = (isSym && triangle == Triangle.Lower),
    }


    static if (!isGen)
    {
        static assert (isNumeric!T || isComplex!T,
            "MatrixView: Non-general matrices can only contain numeric values. "
           ~"Non-appropriate type given: "~T.stringof);
    }

    // The zero element in triangular matrices.
    static if (isTri)
    {
        T zero = Zero!T;
    }


public:
    /** The array that is wrapped by this MatrixView. */
    T[] array;


    /** The number of rows in the matrix. */
    size_t rows;

    /** The number of columns in the matrix. */
    size_t cols;

    /** The leading (row) dimension.  Included to support matrix slicing,
        currently just an alias to rows.
    */
    alias rows leading;


    /**
    Wraps a MatrixView with m rows around the given array.

    For general matrices, the number of columns in the matrix
    is set to a.length/m, whereas for triangular and symmetric
    matrices the number of columns is set equal to the number
    of rows.
    */
    this (T[] a, size_t m)  pure nothrow
    in
    {
        static if (isGen)  assert (a.length % m == 0);
    }
    body
    {
        static if (isGen)  this (a, m, a.length/m);
        else static if (isTri || isSym)  this(a, m, m);
    }



    /**
    Wraps a MatrixView with m rows and n columns around the given array.
    For a given set of a, m, and n, the following must be true for
    a general matrix:
    ---
    a.length >= m*n
    ---
    For triangular and symmetric matrices, there are two constraints:
    ---
    m == n
    a.length >= (n*n + n)/2
    ---
    These conditions are only checked in non-release builds.
    */
    this (T[] a, size_t m, size_t n) pure nothrow
    in
    {
        static if (isGen)
            assert (a.length >= m*n);
        else static if (isTri || isSym)
        {
            assert (m == n);
            assert (a.length >= (n*n + n)/2);
        }
    }
    body
    {
        array = a;
        rows = m;
        cols = n;
    }


    /**
    Returns (a reference to) the element at row i, column j.

    Warning:
    For convenience, this method returns values by reference. This
    means that one can do stuff like this:
    ---
    m[1,2] += 3.14;
    ---
    Unfortunately, it also means that in a triangular matrix one
    can change the zero element (which is common for all zero elements).
    ---
    assert ((m[1,0] == 0.0)  &&  m[2,0] == 0.0);
    m[1,0] += 3.14;
    assert (m[2,0] == 3.14);    // passes
    ---
    You are hereby warned.
    */
    ref T opIndex(size_t i, size_t j) pure nothrow
    in
    {
        assert (i < rows  &&  j < cols);
    }
    body
    {
        static if (isGen)
            return array.ptr[i + rows*j];
        else static if (isUpTri)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2];
            else return zero;
        }
        else static if (isLoTri)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2];
            else return zero;
        }
        else static if (isUpSym)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2];
            else return array.ptr[j + (i*i+i)/2];
        }
        else static if (isLoSym)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2];
            else return array.ptr[j + ((rows+rows-i-1)*i)/2];
        }
        else static assert (false);
    }


    /**
    Assigns a value to the element at row i, column j.

    Unlike opIndex(), this method checks that zero elements in
    a triangular matrix aren't assigned to, but only in non-release
    builds.
    */
    T opIndexAssign(T value, size_t i, size_t j) nothrow
    in
    {
        assert (i < rows  &&  j < cols);
        static if (isUpTri)  assert (i <= j);
        static if (isLoTri)  assert (i >= j);
    }
    body
    {
        static if (isGen)
            return array.ptr[i + rows*j] = value;
        else static if (isUpTri)
            return array.ptr[i + (j*j+j)/2] = value;
        else static if (isLoTri)
            return array.ptr[i + ((rows+rows-j-1)*j)/2] = value;
        else static if (isUpSym)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2] = value;
            else  return array.ptr[j + (i*i+i)/2] = value;
        }
        else static if (isLoSym)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2] = value;
            else return array.ptr[j + ((rows+rows-i-1)*i)/2] = value;
        }
        else static assert (false);
    }

    unittest
    {
        alias MatrixView!real GeneralMatrix;
        real[] g = [1.0L, 2, 3, 4, 5, 6];

        auto gm1 = GeneralMatrix(g, 2);
        auto gm2 = GeneralMatrix(g, 3);

        assert (gm1.cols == 3);
        assert (gm2.cols == 2);

        assert (gm1[1,0] == 2);
        assert (gm2[1,0] == 2);

        assert (gm1[1,1] == 4);
        assert (gm2[1,1] == 5);

        gm2[1,1] += 1; assert (gm2[1,1] == 6);
        gm2[1,1] = 10; assert (gm2[1,1] == 10);


        alias MatrixView!(real, Storage.Triangular) UTMatrix;
        real[] u = [1.0, 2, 3, 4, 5, 6];

        auto um1 = UTMatrix(u, 3);
        assert (um1.cols == 3);
        assert (um1[1,0] == 0.0);
        assert (um1[1,1] == 3.0);
        um1[0,2] += 3; assert (u[3] == 7);
        um1[2,2] = 10; assert (u[5] == 10);


        alias MatrixView!(real, Storage.Triangular, Triangle.Lower) LTMatrix;
        real[] l = [1.0, 2, 3, 4, 5, 6];

        auto lm1 = LTMatrix(l, 3);
        assert (lm1.cols == 3);
        assert (lm1[0,1] == 0.0);
        assert (lm1[1,1] == 4.0);
        lm1[2,0] += 4; assert (l[2] == 7);
        lm1[2,2] = 10; assert (l[5] == 10);


        alias MatrixView!(real, Storage.Symmetric) USMatrix;
        real[] us = [1.0, 2, 3, 4, 5, 6];

        auto usm1 = USMatrix(us, 3);
        assert (usm1.cols == 3);
        assert (usm1[1,2] == 5.0);
        foreach (i; 0 .. usm1.rows)
            foreach (j; 0 .. i)
                assert (usm1[i,j] == usm1[j,i]);
        usm1[0,2] += 3; assert (usm1[2,0] == 7);
        usm1[1,2] = 10; assert (usm1[2,1] == 10);


        alias MatrixView!(real, Storage.Symmetric, Triangle.Lower) LSMatrix;
        real[] ls = [1.0, 2, 3, 4, 5, 6];

        auto lsm1 = LSMatrix(ls, 3);
        assert (lsm1.cols == 3);
        assert (lsm1[1,2] == 5.0);
        foreach (i; 0 .. lsm1.rows)
            foreach (j; 0 .. i)
                assert (lsm1[i,j] == lsm1[j,i]);
        lsm1[0,2] += 3; assert (lsm1[2,0] == 6);
        lsm1[1,2] = 10; assert (lsm1[2,1] == 10);
    }

    /**
    Sum, substract, product, division and exponentiation operations.
    */
    ref MatrixView opOpAssign(string op)(MatrixView rhs)
    if (op == "+" || op == "-" || op == "*" || op == "/" || op == "^^")
    in
    {
        assert(rows == rhs.rows && cols == rhs.cols);
    }
    body
    {
        mixin("array[] " ~ op ~ "= rhs.array[];");
        return this;
    }

    unittest
    {
        void test(NumericType)() {
            auto m1 = MatrixView!NumericType([0, 1, 2, 3], 2, 2);
            auto m2 = MatrixView!NumericType([3, 4, 5, 6], 2, 2);
            m1 += m2;
            assert(m1.array == [3, 5, 7, 9]);
            m1 -= m2;
            assert(m1.array == [0, 1, 2, 3]);
            m1 *= m2;
            assert(m1.array == [0, 4, 10, 18]);
            m1 /= m2;
            assert(m1.array == [0, 1, 2, 3]);
            m1 ^^= m2;
            assert(m1.array == [0, 1, 32, 729]);
        }
        test!float();
        test!double();
        test!int();
    }

    /**
    Apply an scalar a matrix by an scalar factor
    */
    ref MatrixView opOpAssign(string op, RightType)(RightType scalar)
    if (op == "+" || op == "-" || op == "*" || op == "/" || op == "^^")
    in 
    {
        static assert(isNumeric!RightType);
    }
    body
    {
        mixin("array[] " ~ op ~ "= scalar;");
        return this;
    }

    unittest
    {
        void test(NumericType)() {
            auto m = MatrixView!NumericType([0, 1, 2, 3], 2, 2);
            m += 1;
            assert(m.array == [1, 2, 3, 4]);
            m -= 1;
            assert(m.array == [0, 1, 2, 3]);
            m *= 2;
            assert(m.array == [0, 2, 4, 6]);
            m /= 2;
            assert(m.array == [0, 1, 2, 3]);
            m ^^= 3;
            assert(m.array == [0, 1, 8, 27]);
        }
        test!float();
        test!double();
        test!int();
    }

    /**
    Use opOpAssign methods to generate equivalent opBinary operators
    */
    MatrixView opBinary(string op, RightType)(RightType rhs)
    body
    {
        auto matrix = copy(this);
        matrix.opOpAssign!op(rhs);
        return matrix;
    }

    /**
    Use opOpAssign methods to generate equivalent opBinaryRight operators

    TODO: simplify this code (it looks like it can be easier).
    */
    MatrixView opBinaryRight(string op, LeftType)(LeftType lhs)
    body
    {
        static if (op == "/" || op == "-" || op == "^^") {
            auto matrix = copy(this);
            for (int i = 0; i < matrix.array.length; ++i) {
                mixin("matrix.array[i] = lhs " ~ op ~ " matrix.array[i];");
            }
            return matrix;
        }
        else {
            return opBinary!op(lhs);
        }
    }

    unittest
    {
        void test(NumericType)() {
            auto m1 = MatrixView!NumericType([1, 1, 2, 3], 2, 2);
            auto m2 = MatrixView!NumericType([3, 4, 6, 6], 2, 2);
            assert((m1 + m2).array == [ 4,  5,  8,  9]);
            assert((m1 * m2).array == [ 3,  4, 12, 18]);
            assert((m1 - m2).array == [-2, -3, -4, -3]);
            assert((m2 / m1).array == [ 3,  4,  3,  2]);
            assert((m1 +  1).array == [ 2,  2,  3,  4]);
            assert((m1 *  2).array == [ 2,  2,  4,  6]);
            assert((m1 -  1).array == [ 0,  0,  1,  2]);
            assert((m1 ^^ 3).array == [ 1,  1,  8, 27]);
            static if (isIntegral!NumericType) {
                assert((m2 / 2).array == [  1, 2, 3, 3]);
            } else {
                assert((m2 / 2).array == [1.5, 2, 3, 3]);
            }
            assert((1 +  m1).array == [ 2,  2,  3,  4]);
            assert((2 *  m1).array == [ 2,  2,  4,  6]);
            assert((1 -  m1).array == [ 0,  0, -1, -2]);
            assert((6 /  m1).array == [ 6,  6,  3,  2]);
            assert((3 ^^ m1).array == [ 3,  3,  9, 27]);
        }
        test!float();
        test!double();
        test!int();
    }

}


/**
Matrix-matrix inner product.

TODO: Extend to other kinds of storage.
TODO: Optimize to iterate the arrays in sequential order.
TODO: Make sure the calls to the opIndex are inlined
TODO: Use BLAS.
*/

MatrixView!(T, Storage.General) dotProduct(T)(MatrixView!(T, Storage.General) a,
                                              MatrixView!(T, Storage.General) b)
in
{
    assert(a.cols == b.rows);
}
body
{
    auto c = MatrixView!T(new T[a.rows * b.cols], a.rows, b.cols);
    for (int i = 0; i < c.rows; ++i) {
        for (int j = 0; j < c.cols; ++j) {
            c[i,j] = 0;
            for (int k = 0; k < a.cols; ++k) {
                c[i,j] += a[i,k] * b[k,j];
            }
        }
    }
    return c;
}


unittest
{
    //             | 3 6 |
    // | 1 2 4 |   |     |   | 35 16 |
    // |       | * | 4 1 | = |       |
    // | 1 3 5 |   |     |   | 45 19 |
    //             | 6 2 |
    auto m1 = MatrixView!double([1, 1, 2, 3, 4, 5], 2, 3);
    auto m2 = MatrixView!double([3, 4, 6, 6, 1, 2], 3, 2);
    auto m3 = dotProduct(m1, m2);
    assert(m3.array == [35, 45, 16, 19]);
}


/**
Evaluates to true if the given type is an instantiation of
MatrixView. Optionally test the element type and/or storage
scheme.
*/
template isMatrixView(MatrixT)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = true;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(MatrixView!int));
    static assert (!isMatrixView!int);
}


/// ditto
template isMatrixView(MatrixT, ElemT)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = is(E == ElemT);
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(MatrixView!int, int));
    static assert (!isMatrixView!(MatrixView!int, float));
}


/// ditto
template isMatrixView(MatrixT, Storage stor)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = S == stor;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(
        MatrixView!(int, Storage.Triangular),
        Storage.Triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.Triangular),
        Storage.Symmetric));
}


/// ditto
template isMatrixView(MatrixT, ElemT, Storage stor)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = is(E == ElemT)  &&  S == stor;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(
        MatrixView!(int, Storage.Triangular),
        int, Storage.Triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.Triangular),
        float, Storage.Triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.Triangular),
        int, Storage.Symmetric));
}
