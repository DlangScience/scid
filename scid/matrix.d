/** This module contains the MatrixView class as well as some
    functions related to it.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.matrix;


//private import std.conv;
import std.string: format, repeat;
import std.traits;

import scid.core.meta;
import scid.core.traits;

version(unittest) {
    import scid.core.testing; 
    import std.math;
}




/** Various matrix representations. */
enum Storage
{
    General,        /// General (dense) matrices
    Triangular,     /// Packed storage of triangular matrices
    Symmetric,      /// Packed storage of symmetric matrices
}


/** In packed storage (triangular, symmetric, and Hermitian matrices),
    one can choose to store either the upper or lower triangle.
*/
enum Triangle : char
{
    Upper = 'U',    /// Store upper triangle
    Lower = 'L'     /// Store lower triangle
}




/** A convenience function that allocates heap memory for a matrix,
    optionally sets the values of the matrix elements, and returns
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


///ditto
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
    auto dense1 = matrix!real(3, 4);
    check (dense1.rows == 3  &&  dense1.cols == 4);
    check (isNaN(dense1[1,2]));

    auto dense2 = matrix(4, 3, 1.0);
    check (dense2.rows == 4  &&  dense2.cols == 3);
    check (dense2[1,2] == 1.0);

    auto upTri = matrix!(real, Storage.Triangular)(3);
    check (upTri.rows == 3  &&  upTri.cols == 3);
    check (isNaN(upTri[0,2])  &&  upTri[2,0] == 0);

    auto loTri = matrix!(double, Storage.Triangular, Triangle.Lower)(3, 1.0);
    check (loTri.rows == 3  &&  loTri.cols == 3);
    check (loTri[0,2] == 0  &&  loTri[2,0] == 1.0);

}




/** A convenience function that creates a copy of the input matrix. */
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
    check (b[1,0] == 2.0  &&  a[1,0] == 1.0);
}




/** This struct provides a matrix-like view of the contents of an
    array. In order to be compatible with LAPACK routines, it supports
    the following matrix representations (i.e. memory layouts).
    
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
    b c e   or  b d c
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



    /** Wrap a MatrixView with m rows around the given array.
        
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



    /** Wrap a MatrixView with m rows and n columns around the given array.
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

    MatrixView!(T) opBinary(string op)(MatrixView!(T) rhs) pure
    in
    {
        static if (op == "+" || op == "-") {
            assert (rows == rhs.rows && cols == rhs.cols);
        } else static if (op == "*") {
            assert (cols == rhs.rows);
        }
        
    }
    body
    {
        static if (op == "+" || op == "-") {
            T[] a = array.dup;
            foreach( ulong i; 0..a.length ) {
                a[i] = mixin("a[i] "~op~" rhs.array[i]");
            }
            return MatrixView!(T)(a, rows, cols);
        } else static if (op == "*") {
            ulong r = rows;
            ulong c = rhs.cols;
            T[] a;
            a.length = r * c;
            foreach( ulong j; 0..r ) {
                foreach( ulong i; 0..c ) {
                    a[j * c + i] = 0;
                    foreach( size_t k; 0..cols ) {
                        a[j * c + i] += array[j * c + k] * rhs.array[k * c + i];
                    }
                }
            }
            return MatrixView!(T)(a, r, c);
        } else static assert(0, "Operator "~op~" not implemented");
    }
    
    MatrixView!(T) opBinary(string op)(T rhs) pure
    in
    {
    }
    body
    {
        static if (op == "*" || op == "/") {
            T[] a = array.dup;
            foreach( ulong i; 0..a.length ) {
                a[i] = mixin("a[i] "~op~" rhs");
            }
            return MatrixView!(T)(a, rows, cols);
        } else static assert(0, "Operator "~op~" not implemented");
    }
    
    /** Return (a reference to) the element at row i, column j.

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


    /** Assign a value to the element at row i, column j.
        
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
}

unittest
{
    alias MatrixView!real GeneralMatrix;
    real[] g = [1.0L, 2, 3, 4, 5, 6];

    auto gm1 = GeneralMatrix(g, 2);
    auto gm2 = GeneralMatrix(g, 3);

    check (gm1.cols == 3);
    check (gm2.cols == 2);

    check (gm1[1,0] == 2);
    check (gm2[1,0] == 2);

    check (gm1[1,1] == 4);
    check (gm2[1,1] == 5);

    gm2[1,1] += 1; check (gm2[1,1] == 6);
    gm2[1,1] = 10; check (gm2[1,1] == 10);


    alias MatrixView!(real, Storage.Triangular) UTMatrix;
    real[] u = [1.0, 2, 3, 4, 5, 6];
    
    auto um1 = UTMatrix(u, 3);
    check (um1.cols == 3);
    check (um1[1,0] == 0.0);
    check (um1[1,1] == 3.0);
    um1[0,2] += 3; check (u[3] == 7);
    um1[2,2] = 10; check (u[5] == 10);


    alias MatrixView!(real, Storage.Triangular, Triangle.Lower) LTMatrix;
    real[] l = [1.0, 2, 3, 4, 5, 6];
    
    auto lm1 = LTMatrix(l, 3);
    check (lm1.cols == 3);
    check (lm1[0,1] == 0.0);
    check (lm1[1,1] == 4.0);
    lm1[2,0] += 4; check (l[2] == 7);
    lm1[2,2] = 10; check (l[5] == 10);


    alias MatrixView!(real, Storage.Symmetric) USMatrix;
    real[] us = [1.0, 2, 3, 4, 5, 6];
    
    auto usm1 = USMatrix(us, 3);
    check (usm1.cols == 3);
    check (usm1[1,2] == 5.0);
    foreach (i; 0 .. usm1.rows)
        foreach (j; 0 .. i)
            check (usm1[i,j] == usm1[j,i]);
    usm1[0,2] += 3; check (usm1[2,0] == 7);
    usm1[1,2] = 10; check (usm1[2,1] == 10);


    alias MatrixView!(real, Storage.Symmetric, Triangle.Lower) LSMatrix;
    real[] ls = [1.0, 2, 3, 4, 5, 6];
    
    auto lsm1 = LSMatrix(ls, 3);
    check (lsm1.cols == 3);
    check (lsm1[1,2] == 5.0);
    foreach (i; 0 .. lsm1.rows)
        foreach (j; 0 .. i)
            check (lsm1[i,j] == lsm1[j,i]);
    lsm1[0,2] += 3; check (lsm1[2,0] == 6);
    lsm1[1,2] = 10; check (lsm1[2,1] == 10);
}




/** Evaluates to true if the given type is an instantiation of
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


