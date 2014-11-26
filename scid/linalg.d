/** This module contains functions related to linear algebra.

    For the time being these functions are just a user-friendly interface
    to the corresponding LAPACK functions. Hence, one can only use matrices
    and vectors where the elements are of FORTRAN-compatible types, namely:
    ---
    float
    double
    cfloat
    cdouble
    ---
    Specifically, real- and creal-valued matrices/vectors cannot be used.

    Note:
    Some of the functions in this module come in two forms with the same
    function signature, but where the name of one ends with an underscore.
    An example is
    ---
    solve(a, b)
    solve_(a, b)
    ---
    The difference between these is that, for performance reasons, the
    underscore-suffixed functions use some or all of the input
    matrices/vectors as a workspace, and one can't expect them to
    contain the same values after the function returns. The functions
    without an underscore suffix simply copy the input data and calls the
    high-performance function using the copied data as input.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.linalg;


import std.algorithm;
import std.complex;
import std.conv;
import std.exception;
import std.math;
import std.traits;

import scid.bindings.lapack.dlapack;
import scid.core.fortran;
import scid.core.memory;
import scid.core.meta;
import scid.core.testing;
import scid.core.traits;
import scid.matrix;
import scid.util;

pragma(lib, "blas");
pragma(lib, "lapack");




/** Solve one or more systems of n linear equations in n variables.

    The set of equations is given on the form AX=B, where A
    is an n-by-n matrix and X and B are either vectors of
    length n (one system of n equations in n variables) or n-by-m
    matrices (m systems of n equations in n variables). Given A
    and B as input, this function returns X.

    Examples:
    Solving a system of equations:
    ---
    MatrixView!double a = ...
    double[] b = ...
    double[] x = solve(a, b);
    ---
    Solving several systems of equations:
    ---
    MatrixView!double a = ...
    MatrixView!double b = ...
    MatrixView!double x = solve(a, b);
    ---
*/
Real[] solve (MatrixViewA, Real)
    (const MatrixViewA a, const Real[] b, Real[] buffer=null)
    if (isMatrixView!MatrixViewA  &&  isFortranType!Real)
{
    mixin (newFrame);

    static assert (is (BaseElementType!MatrixViewA == Real),
        "solve: a and b have different element types");

    // Make a copy of a in TempAlloc space.
    auto adup = tempdup(tailConst(a.array));

    // b is copied into the buffer.
    buffer.length = b.length;
    buffer[] = b[];

    solveImpl!(Real, a.storage, a.triangle)
        (adup, a.rows, buffer, buffer.length, 1);
    return buffer;
}


/// ditto
Real[] solve_ (MatrixViewA, Real)
    (MatrixViewA a, Real[] b)
    if (isMatrixView!MatrixViewA  &&  isFortranType!Real)
{
    static assert (is (BaseElementType!MatrixViewA == Real),
        "solve: a and b have different element types");

    solveImpl!(Real, a.storage, a.triangle)
        (a.array, a.rows, b, b.length, 1);
    return b;
}


/// ditto
MatrixViewB solve (MatrixViewA, MatrixViewB)
    (const MatrixViewA a, const MatrixViewB b)
    if (isMatrixView!(MatrixViewA)
     && isMatrixView!(MatrixViewB, Storage.General))
{
    mixin (newFrame);

    alias BaseElementType!MatrixViewB Real;
    static assert (is (BaseElementType!MatrixViewA == Real),
        "solve: a and b have different element types");
    static assert (isFortranType!Real,
        "solve: Not a FORTRAN-compatible type: "~Real.stringof);

    // Make a copy of a in TempAlloc space.
    auto adup = tempdup(tailConst(a.array));

    // Make an ordinary copy of b.
    auto bdup = copy(b);

    solveImpl!(Real, a.storage, a.triangle)
        (adup, a.rows, bdup.array, bdup.rows, bdup.cols);
    return bdup;
}


/// ditto
MatrixViewB solve_ (MatrixViewA, MatrixViewB)
    (MatrixViewA a, MatrixViewB b)
    if (isMatrixView!(MatrixViewA)
     && isMatrixView!(MatrixViewB, Storage.General))
in
{
    assert (a.rows == a.cols);
    assert (a.cols == b.rows);
}
body
{
    alias BaseElementType!MatrixViewB Real;

    static assert (is (BaseElementType!MatrixViewA == Real),
        "solve: a and b have different element types");
    static assert (isFortranType!Real,
        "solve: Not a FORTRAN-compatible type: "~Real.stringof);

    solveImpl!(Real, a.storage, a.triangle)
        (a.array, a.rows, b.array, b.rows, b.cols);
    return b;
}


// This is where the solve() magic happens.
private void solveImpl(Real, Storage aStorage, Triangle aTriangle)
    (Real[] aMatrix, size_t aRows, Real[] bxMatrix, size_t bRows, size_t bCols)
{
    mixin (newFrame);

    int info;
    static if (aStorage == Storage.General)
    {
        int* ipiv = cast(int*) TempAlloc.malloc(aRows*int.sizeof);
        gesv(
            toInt(aRows),   // N
            toInt(bCols),   // NRHS
            aMatrix.ptr,    // A
            toInt(aRows),   // LDA
            ipiv,           // IPIV
            bxMatrix.ptr,   // B
            toInt(bRows),   // LDB
            info);          // INFO
    }
    else static if (aStorage == Storage.Symmetric)
    {
        int* ipiv = cast(int*) TempAlloc.malloc(aRows*int.sizeof);
        spsv(
            aTriangle,      // UPLO
            toInt(aRows),   // N
            toInt(bCols),   // NRHS
            aMatrix.ptr,    // AP
            ipiv,           // IPIV
            bxMatrix.ptr,   // B
            toInt(bRows),   // LDB
            info);          // INFO
    }
    else static if (aStorage == Storage.Triangular)
    {
        enum char trans = 'N';  // Matrix isn't transposed.
        enum char diag = 'N';   // We don't know that the matrix
                                // is unit triangular.

        tptrs(
            aTriangle,      // UPLO
            trans,          // TRANS
            diag,           // DIAG
            toInt(aRows),   // N
            toInt(bCols),   // NRHS
            aMatrix.ptr,    // AP
            bxMatrix.ptr,   // B
            toInt(bRows),   // LDB
            info);          // INFO

    }
    else static assert (false, "solve: Unsupported matrix storage.");

    assert (info >= 0);
    enforce(info == 0, "solve: matrix is singular");
    return;
}

unittest
{
    alias solve!(MatrixView!float, MatrixView!float) solve_f;
    alias solve!(MatrixView!double, MatrixView!double) solve_d;
    alias solve!(MatrixView!cfloat, MatrixView!cfloat) solve_cf;
    alias solve!(MatrixView!cdouble, MatrixView!cdouble) solve_cd;

    alias solve!(MatrixView!(float, Storage.Symmetric), MatrixView!float) solve_f_s;
    alias solve!(MatrixView!(double, Storage.Symmetric), MatrixView!double) solve_d_s;
    alias solve!(MatrixView!(cfloat, Storage.Symmetric), MatrixView!cfloat) solve_cf_s;
    alias solve!(MatrixView!(cdouble, Storage.Symmetric), MatrixView!cdouble) solve_cd_s;

    alias solve!(MatrixView!(float, Storage.Triangular), MatrixView!float) solve_f_t;
    alias solve!(MatrixView!(double, Storage.Triangular), MatrixView!double) solve_d_t;
    alias solve!(MatrixView!(cfloat, Storage.Triangular), MatrixView!cfloat) solve_cf_t;
    alias solve!(MatrixView!(cdouble, Storage.Triangular), MatrixView!cdouble) solve_cd_t;

    alias solve!(MatrixView!double, double) solve_vec;
}

unittest
{
    alias solve_!(MatrixView!float, MatrixView!float) solve_f;
    alias solve_!(MatrixView!double, MatrixView!double) solve_d;
    alias solve_!(MatrixView!cfloat, MatrixView!cfloat) solve_cf;
    alias solve_!(MatrixView!cdouble, MatrixView!cdouble) solve_cd;

    alias solve_!(MatrixView!(float, Storage.Symmetric), MatrixView!float) solve_f_s;
    alias solve_!(MatrixView!(double, Storage.Symmetric), MatrixView!double) solve_d_s;
    alias solve_!(MatrixView!(cfloat, Storage.Symmetric), MatrixView!cfloat) solve_cf_s;
    alias solve_!(MatrixView!(cdouble, Storage.Symmetric), MatrixView!cdouble) solve_cd_s;

    alias solve_!(MatrixView!(float, Storage.Triangular), MatrixView!float) solve_f_t;
    alias solve_!(MatrixView!(double, Storage.Triangular), MatrixView!double) solve_d_t;
    alias solve_!(MatrixView!(cfloat, Storage.Triangular), MatrixView!cfloat) solve_cf_t;
    alias solve_!(MatrixView!(cdouble, Storage.Triangular), MatrixView!cdouble) solve_cd_t;

    alias solve_!(MatrixView!double, double) solve_vec;
}

unittest
{
    // A is dense.
    double[] daa = [ 1.0, 2, -1, -1, 2, 5, 1, -4, 0 ];
    double[] dba = [ 2.0, -6, 9, 0, -6, 1 ];
    double[] dxa = [ 1.0, 2, 3, -1, 0, 1 ];
    auto da = MatrixView!double(daa, 3);
    auto db = MatrixView!double(dba, 3);
    auto dx = solve(da, db);
    check (approxEqual(dx.array, dxa, double.epsilon));

    // Vector version.
    double[] dbva = [ 2.0, -6, 9 ];
    double[] dxva = [ 1.0,  2, 3 ];
    auto dxv = solve(da, dbva);
    check (approxEqual(dxv, dxva, double.epsilon));

    // A is symmetric, packed storage.
    float[] saa = [1.0, 2, 4, -3, 5, -6];
    float[] sba = [-4.0, 25, -11, -4, 58, -23];
    float[] sxa = [1.0, 2, 3, 4, 5, 6];
    auto sa = MatrixView!(float, Storage.Symmetric)(saa, 3);
    auto sb = MatrixView!(float)(sba, 3);
    auto sx = solve(sa, sb);
    check (approxEqual(sx.array, sxa, float.epsilon));

    // A is triangular, packed storage.
    double[] taa = [1.0, 2, 4, -3, 5, -6];
    double[] tba = [-4.0, 23, -18, -4, 50, -36];
    double[] txa = [1.0, 2, 3, 4, 5, 6];
    auto ta = MatrixView!(double, Storage.Triangular)(taa, 3);
    auto tb = MatrixView!(double)(tba, 3);
    auto tx = solve(ta, tb);
    check (approxEqual(tx.array, txa, double.epsilon));
}




/** Calculate the determinant of a square matrix.

    The type of the return value depends on the element type of
    the matrix. For float or double matrices the determinant is of
    type real, and for cfloat or cdouble matrices the determinant
    is of type creal. The reason for choosing the widest type is that
    determinants are often very big numbers, and therefore tend
    to overflow.

    Examples:
    ---
    import scid.matrix;
    ...
    auto m = matrix!double(2, 2);
    m[0,0] = 1.0;  m[0,1] = 2.0;
    m[1,0] = 3.0;  m[1,1] = 4.0;

    auto d = det(m);
    writeln(d);   // Prints "-2"
    ---
*/
DetType!(MatrixView!(T, stor)) det
    (T, Storage stor)
    (const MatrixView!(T, stor) m)
{
    return det_(copy(m));
}


/// ditto
DetType!(MatrixView!(T, stor)) det_
    (T, Storage stor)
    (MatrixView!(T, stor) m)
in
{
    assert (m.rows == m.cols);
}
body
{
    mixin (newFrame);

    static assert (isFortranType!(T),
        "det: Not a FORTRAN-compatible type: "~T.stringof);

    alias typeof(return) DT;

    // LU factorise matrix.
    int info;
    static if (m.storage == Storage.General)
    {
        auto ipiv = newStack!int(m.rows);
        getrf(toInt(m.rows), toInt(m.cols), m.array.ptr, toInt(m.rows), ipiv.ptr, info);
        assert (info >= 0, "invalid input to getrf");
        
        // If matrix is singular, determinant is zero.
        if (info > 0)  return Zero!DT;

        // The determinant is the product of the diagonal entries
        // of the upper triangular matrix. The array ipiv contains
        // the pivots.
        DT d = One!DT;
        for (int i=0; i<m.rows; i++)
        {
            auto p = ipiv[i];
            if (p == i+1) d *= m[i,i];  // i.e. row interchanged with itself
            else          d *= -m[i,i]; // i.e. row interchanged with another
        }

        return d;
    }

    else static if (m.storage == Storage.Symmetric)
    {
        auto ipiv = newStack!int(m.rows);
        sptrf(m.triangle, toInt(m.rows), m.array.ptr, ipiv.ptr, info);
        assert (info >= 0, "invalid input to sptrf");
        
        // If matrix is singular, determinant is zero.
        if (info > 0)  return Zero!DT;

        // Calculate determinant.
        DT d = One!DT;
        for (int k=0; k<m.rows; k++)
        {
            auto p = ipiv[k];
            if (p > 0)  // 1x1 block at m[k,k]
            {
                if (p == k+1) d *= m[k,k];  // row interchanged with self
                else          d *= -m[k,k]; // row interchanged with other
            }
            else // p < 0, 2x2 block at m[k..k+1, k..k+1]
            {
                auto kp1 = k+1;
                auto offDiag = m[k, kp1];
                auto blockDet = m[k,k]*m[kp1,kp1] - offDiag*offDiag;
                if (-p == kp1 || -p == k+2) // row interchanged with self
                    d *= blockDet;
                else                        // row interchanged with other
                    d *= -blockDet;
                k++;
            }
        }

        return d;
    }
    else static if (m.storage == Storage.Triangular)
    {
        DT d = m[0,0];
        foreach (i; 1 .. m.rows)  d *= m[i,i];
        return d;
    }
    else static assert (false, "det: Unsupported matrix storage.");
}


template DetType(MatrixT)
{
    static if (scid.core.traits.isComplex!(BaseElementType!MatrixT))
        alias creal DetType;
    else
        alias real DetType;
}


unittest
{
    // Check for all FORTRAN compatible types.
    alias det!(float, Storage.General) det_fd;
    alias det!(double, Storage.General) det_dd;
    alias det!(cfloat, Storage.General) cfdet_cfd;
    alias det!(cdouble, Storage.General) cddet_cdd;

    alias det!(float, Storage.Symmetric) det_fs;
    alias det!(double, Storage.Symmetric) det_ds;
    alias det!(cfloat, Storage.Symmetric) cfdet_cfs;
    alias det!(cdouble, Storage.Symmetric) cddet_cds;

    // Check for zero-determinant shortcut.
    double[] dsinga = [4.0, 2, 2, 1];   // dense singular
    auto dsing = MatrixView!double(dsinga, 2);
    auto dsingd = det(dsing);
    check (dsingd == 0.0);

    double[] ssinga = [4.0, 2, 1];      // symmetric packed singular
    auto ssing = MatrixView!(double, Storage.Symmetric)(ssinga, 2);
    auto ssingd = det(ssing);
    check (ssingd == 0.0);


    // General dense matrix.
    int dn = 101;
    auto d = matrix!double(dn, dn);
    for (int k=0; k<dn; k++)
    {
        for (int l=0; l<dn; l++)
        {
            if (k == l)
                d[k,l] = (k+1)*(k+1) + 1.0;
            else
                d[k,l] = 2.0*(k+1)*(l+1);
        }
    }
    auto dd = det(d);
    check (approxEqual(dd, 8.972817920259982e319L, sqrt(real.epsilon)));


    // Symmetric packed matrix
    double[] spa = [ 1.0, -2, 3, 4, 5, -6, -7, -8, -9, 10];
    auto sp = MatrixView!(double, Storage.Symmetric)(spa, 4);
    auto spd = det(sp);
    check (approxEqual(spd, 5874.0, sqrt(double.epsilon)));


    // Triangular packed matrix
    double[] tpa = [ 1.0, -2, 3, 4, 5, -6, -7, -8, -9, 10];
    auto tp = MatrixView!(double, Storage.Triangular)(tpa, 4);
    auto tpd = det(tp);
    check (approxEqual(tpd, -180.0, sqrt(double.epsilon)));
}




/** Calculate the eigenvalues of a general dense square matrix.

    If some eigenvalues cannot be calculated, the algorithm
    throws an EigenvalueException containing an array of the ones
    that have been calculated.

    Params:
        m = An n-by-n symmetric matrix.
        buffer = (optional) A buffer for the returned values, must
            have length >= n and type Complex!T[].

    Examples:
    ---
    auto m = matrix!double(3, 3);
    ...
    auto e = eigenvalues(m);
    ---
*/
EigenvalueType!ElementT[] eigenvalues(ElementT, Storage stor, Triangle tri)
    (MatrixView!(ElementT, stor, tri) m)
    if (stor == Storage.General)
{
    return eigenvalues_(copy(m));
}

/// ditto
ComplexT[] eigenvalues(ElementT, ComplexT, Storage stor, Triangle tri)
    (MatrixView!(ElementT, stor, tri) m, ComplexT[] buffer)
    if (stor == Storage.General)
{
    return eigenvalues_(copy(m), buffer);
}

/// ditto
EigenvalueType!ElementT[] eigenvalues_(ElementT, Storage stor, Triangle tri)
    (MatrixView!(ElementT, stor, tri) m)
    if (stor == Storage.General)
{
    static assert (isFortranType!ElementT,
        "eigenvalues: Not a FORTRAN-compatible type: "~T.stringof);

    static if (scid.core.traits.isComplex!ElementT)
        return eigenvaluesComplex_!ElementT(m, null);
    else
        return eigenvaluesReal_!ElementT(m, null);
}

/// ditto
ComplexT[] eigenvalues_(ElementT, ComplexT, Storage stor, Triangle tri)
    (MatrixView!(ElementT, stor, tri) m, ComplexT[] buffer)
    if (stor == Storage.General)
{
    static assert (scid.core.traits.isComplex!ComplexT,
        "eigenvalues: Not a complex type: "~ComplexT.stringof);
    static assert (is(ElementT == ComplexT)
        || is(ElementT == typeof(ComplexT.re)),
        "eigenvalues: The elements of the matrix must be of type "
        ~ComplexT.stringof~" or "~typeof(ComplexT.re).stringof
        ~", and not "~ElementT.stringof);
    static assert (isFortranType!ComplexT,
        "eigenvalues: Not a FORTRAN-compatible type: "~T.stringof);

    static if (scid.core.traits.isComplex!ElementT)
    {
        return eigenvaluesComplex_(m, buffer);
    }
    else
    {
        return eigenvaluesReal_(m, buffer);
    }
}


template EigenvalueType(T)
{
    static if (isFloatingPoint!T) alias Complex!T EigenvalueType;
    else alias T EigenvalueType;
}


private Complex!T[] eigenvaluesReal_ (T)
    (MatrixView!(T) m, Complex!T[] buffer)
    if (isFloatingPoint!T)
in
{
    assert (m.rows == m.cols, "eigenvalues: Matrix must be square");
}
body
{
    mixin (newFrame);

    immutable int n = toInt(m.rows);
    if (n == 0) return null;    // Empty matrix.
    buffer.length = n;

    // Calculate optimal workspace size.
    int info;
    T optimal;
    geev('N', 'N',      // Not going to compute eigenvectors.
        n, null, n,     // Need matrix info.
        null, null,     // Eigenvalues, not calculated.
        null, 1,        // Left eigenvectors, not calculated.
        null, 1,        // Right eigenvectors, not calculated.
        &optimal, -1,   // Query for optimal workspace size.
        info);

    // Allocate workspace and result arrays.
    immutable npn = n+n;
    auto workspace = newStack!T(npn + to!int(optimal));
    auto wr = workspace[0 .. n];
    auto wi = workspace[n .. npn];
    auto work = workspace[npn .. $];

    // Call LAPACK routine GEEV to calculate eigenvalues.
    geev('N', 'N',              // Don't compute eigenvectors.
        n, m.array.ptr, n,      // Input matrix.
        wr.ptr, wi.ptr,         // Eigenvalues.
        null, 1,                // Left eigenvectors, not calculated.
        null, 1,                // Right eigenvectors, not calculated.
        work.ptr, toInt(work.length),  // Workspace.
        info);

    if (info == 0)          // Success!
    {
        foreach (i, ref e; buffer)
        {
            e.re = wr[i];
            e.im = wi[i];
        }
        return buffer;
    }
    else if (info > 0)      // Only some eigenvalues were computed.
    {
        buffer = buffer[0 .. n-info];
        int i = info;
        foreach (ref e; buffer)
        {
            e.re = wr[i];
            e.im = wi[i];
            i++;
        }
        throw new EigenvalueException!(typeof(buffer))(buffer.dup);
    }

    assert (0);
}


private T[] eigenvaluesComplex_ (T)
    (MatrixView!(T) m, T[] buffer)
    if (scid.core.traits.isComplex!T)
in
{
    assert (m.rows == m.cols, "eigenvalues: Matrix must be square");
}
body
{
    mixin (newFrame);

    // Until std.complex finally replaces cfloat and cdouble, we use this
    // little hack:
    static if (is(typeof(T.re) == float)) alias cfloat cT;
    else static if (is(typeof(T.re) == double)) alias cdouble cT;
    else static assert(0);

    immutable int n = toInt(m.rows);
    if (n == 0) return null;    // Empty matrix.
    buffer.length = n;

    // Calculate optimal workspace size.
    int info;
    cT optimal;
    geev('N', 'N',
        n, null, n,     // Need matrix info.
        null,           // Eigenvalues, not calculated.
        null, 1,        // Left eigenvectors, not calculated.
        null, 1,        // Right eigenvectors, not calculated.
        &optimal, -1,   // Query for optimal workspace size.
        null,           // Second workspace, not needed.
        info);

    // Allocate workspace arrays.
    auto rwork = newStack!(typeof(T.re))(2*n);
    auto work = newStack!T(to!int(optimal.re));

    // Call LAPACK routine GEEV to calculate eigenvalues.
    geev('N', 'N',                              // Don't compute eigenvectors.
        n, cast(cT*) m.array.ptr, n,            // Input matrix.
        cast(cT*) buffer.ptr,                   // Eigenvalues.
        null, 1,                    // Left eigenvectors, not calculated.
        null, 1,                    // Right eigenvectors, not calculated.
        cast(cT*) work.ptr, toInt(work.length), // Workspace 1.
        rwork.ptr,                              // Workspace 2.
        info);

    if (info == 0)          // Success!
    {
        return buffer;
    }
    else if (info > 0)      // Only some eigenvalues were computed.
    {
        throw new EigenvalueException!(typeof(buffer))(buffer[info .. n].dup);
    }

    assert (0);
}

unittest
{
    // General double matrix.
    auto m = matrix!double(3, 3);
    foreach (i; 0 .. m.rows)
        foreach (j; 0 .. m.cols)  m[i,j] = i + 2.0*j;
    auto v = eigenvalues(m);
    check (approxEqual(v[0].re, (9+sqrt(129.0))/2, 1e-10) && v[0].im == 0);
    check (approxEqual(v[1].re, (9-sqrt(129.0))/2, 1e-10) && v[1].im == 0);
    check (approxEqual(v[2].re,               0.0, 1e-10) && v[1].im == 0);
}

unittest
{
    // General complex matrix.
    alias Complex!double C;

    auto m = matrix!C(3, 3);
    foreach (i; 0 .. m.rows)
        foreach (j; 0 .. m.cols)  m[i,j] = C(i, j);

    auto ev = eigenvalues(m);

    double a = 1.5;
    double b = cos(PI_4)*sqrt(42.0)*0.5;
    check (approxEqual(ev[0].re, a+b, 1e-10));
    check (approxEqual(ev[1].re, a-b, 1e-10));
    check (approxEqual(ev[2].re, 0.0, 1e-10));
    foreach (e; ev) check (approxEqual(e.re, e.im, 1e-10));
}




/** Calculate the eigenvalues of a triangular matrix. Note that this
    is a trivial operation - the function just returns the diagonal
    entries of the matrix.

    Params:
        m = An n-by-n triangular matrix.
        buffer = (optional) A buffer for the returned values, must
            have length >= n.
*/
T[] eigenvalues (T, Storage stor, Triangle tri)
    (MatrixView!(T, stor, tri) m, T[] buffer=null)
    if (stor == Storage.Triangular)
{
    immutable int n = toInt(m.rows);
    if (n == 0) return null;    // Empty matrix.
    buffer.length = n;

    foreach (i, ref ev; buffer)  ev = m[i,i];
    return buffer;
}

unittest
{
    double[] a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    auto m = MatrixView!(double, Storage.Triangular)(a, 3);
    auto v = eigenvalues(m);
    assert (v.length == 3);
    check (v[0] == 1.0 && v[1] == 3.0 && v[2] == 6.0);
    foreach (e; v)  check (e.im == 0.0);
}




/** Calculate the eigenvalues of a symmetric matrix.

    Params:
        m = An n-by-n symmetric matrix.
        buffer = (optional) A buffer for the returned values, must
            have length >= n.
*/
T[] eigenvalues (T, Storage stor, Triangle tri)
    (MatrixView!(T, stor, tri) m, T[] buffer=null)
    if (isFloatingPoint!T  &&  stor == Storage.Symmetric)
{
    return eigenvalues_(copy(m), buffer);
}

/// ditto
T[] eigenvalues_ (T, Storage stor, Triangle tri)
    (MatrixView!(T, stor, tri) m, T[] buffer=null)
    if (isFloatingPoint!T  &&  stor == Storage.Symmetric)
{
    mixin (newFrame);

    static assert (isFortranType!T,
        "eigenvalues: Not a FORTRAN-compatible type: "~T.stringof);

    immutable int n = toInt(m.rows);
    if (n == 0) return null;    // Empty matrix.
    buffer.length = n;

    // Allocate workspace.
    auto workspace = cast(T*) TempAlloc.malloc(3*n*T.sizeof);

    // Call LAPACK routine SPEV.
    int info;
    spev('N',                               // Don't compute eigenvectors
        m.triangle, toInt(m.rows), m.array.ptr,    // Input matrix.
        buffer.ptr,                         // Eigenvalue array.
        null, 1,                            // Eigenvectors, not calculated.
        workspace,                          // Workspace.
        info);

    if (info == 0)
    {
        return buffer;
    }
    else if (info > 0)
    {
        throw new Exception(text("The algorithm failed to converge. ",
            info, " off-diagonal elements of an intermediate tridiagonal form "
            ~ "did not converge to zero."));
    }

    assert(0);
}

unittest
{
    // Symmetric double matrix.
    double[] a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    auto m = MatrixView!(double, Storage.Symmetric, Triangle.Lower)(a, 3);
    auto v = eigenvalues(m);
    assert (v.length == 3);
    check (approxEqual(v[0], - 0.5157294716, 1e-10));
    check (approxEqual(v[1],   0.1709151888, 1e-10));
    check (approxEqual(v[2],  11.34481428,   1e-10));
}




/** This exception is thrown when the eigenvalue() function fails to
    compute all eigenvalues. The ones that have been computed are stored
    in the member variable eigenvalues.
*/
class EigenvalueException(T) : Exception
{
    /// The computed eigenvalues.
    T eigenvalues;

    this (T v)
    {
        super ("Failed to compute all eigenvalues. "
            ~to!string(v.length)~" eigenvalue(s) have been computed");
        eigenvalues = v;
    }
}




/** Calculate the inverse of a matrix.

    Currently only defined for general real matrices.
*/
void invert(T, Storage stor)(MatrixView!(T, stor) m)
    if (isFortranType!T  &&  !scid.core.traits.isComplex!T
        &&  stor == Storage.General)
in
{
    assert (m.rows == m.cols, "invert: can only invert square matrices");
}
body
{
    mixin (newFrame);

    // Calculate optimal workspace size.
    int info;
    T optimal;
    getri(
        toInt(m.rows), null, toInt(m.leading),  // Info about M
        null, &optimal, -1,                     // Do workspace query
        info);

    // Allocate workspace memory.
    int* ipiv = cast(int*) TempAlloc.malloc(m.rows*int.sizeof);
    T[] work = newStack!T(to!int(optimal));

    // Calculate LU factorisation.
    getrf(
        toInt(m.rows), toInt(m.cols), m.array.ptr, toInt(m.leading),
        ipiv, info);

    // Invert matrix.
    getri(
        toInt(m.rows), m.array.ptr, toInt(m.leading),   // Matrix
        ipiv, work.ptr, toInt(work.length),             // Workspace
        info);

    assert (info >= 0);
    enforce(info == 0, "invert: matrix is singular");
    return;
}

unittest
{
    double[] aa = [1.0, 0, 2, 2, 2, 0, 0, 1, 1];
    auto a = MatrixView!double(aa, 3, 3);
    invert(a);
    enum : real { _13 = 1.0L/3.0, _16 = 1.0L/6.0, _23 = 2.0L/3.0 }
    real[] ans = [_13, _13, -_23, -_13, _16, _23, _13, -_16, _13];
    foreach (i, e; aa)  check (approxEqual(cast(real) e, ans[i], 1e-10L));
}

unittest
{
    double[] aa = new double[9];
    foreach (i, ref e; aa)  e = i;
    auto a = MatrixView!double(aa, 3, 3);
    
    try { invert(a); check(false, "Matrix should be detected as singular"); }
    catch (Exception e) { check (true); }
}
