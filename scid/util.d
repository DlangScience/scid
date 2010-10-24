/** Various useful utilities that don't really fit anywhere else.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.util;


import std.complex;
import std.math;
import std.range;
import std.traits;
import std.typecons;

import scid.core.testing;
import scid.core.traits;
import scid.exception;




/** Check whether lhs and rhs are equal to within the specified number
    of significant digits.  Both lhs and rhs can be floating-point
    numbers, complex numbers, or input ranges of floating-point or
    complex numbers.
    ---
    assert (matchDigits(0.1234567, 0.1234568));
    ---
*/
bool matchDigits(L, R)(L lhs, R rhs, uint significantDigits = 6)
    if ((isFloatingPoint!L || is(L T == Complex!T)) &&
        (isFloatingPoint!R || is(R U == Complex!U)))
in
{
    static if (isFloatingPoint!L && isFloatingPoint!R)
    {
        assert (significantDigits*LOG2T < CommonType!(L, R).mant_dig,
            "The requested precision is too high for the given type(s)");
    }
}
body
{
    static if (is(L T == Complex!T))
    {
        static if (is(R U == Complex!U))
        {
            // lhs and rhs are complex
            return matchDigits(lhs.re, rhs.re, significantDigits)
                && matchDigits(lhs.im, rhs.im, significantDigits);
        }
        else
        {
            // lhs is complex, rhs is real
            return matchDigits(lhs.re, rhs, significantDigits)
                && matchDigits(lhs.im, 0.0, significantDigits);
        }
    }
    else
    {
        static if (is(R U == Complex!U))
        {
            // lhs is real, rhs is complex
            return matchDigits(lhs, rhs.re, significantDigits)
                && matchDigits(0.0, rhs.im, significantDigits);
        }
        else
        {
            // lhs and rhs are real
            return feqrel!(CommonType!(L, R))(lhs, rhs)
                > significantDigits*LOG2T;
        }
    }
}


unittest
{
    check(matchDigits(1.0, 1.0, 15)); 
    check(matchDigits(0.1234567, 0.1234568, 6));
    check(!matchDigits(0.1234567, 0.1234568, 7));

    auto z = Complex!real(0.1234567, 0.0);
    check(matchDigits(z, 0.1234568, 6));
    check(matchDigits(0.1234568, z, 6));

    auto u = Complex!real(0.1234567, 0.1234567);
    auto v = Complex!real(0.1234568, 0.1234567);
    auto w = Complex!real(0.1234567, 0.1234568);
    check(matchDigits(u, v, 6));
    check(matchDigits(u, w, 6));
    check(matchDigits(v, w, 6));
}




/// ditto
bool matchDigits(L, R)(L lhs, R rhs, uint significantDigits = 6)
    if (isInputRange!L || isInputRange!R)
{
    static if (isInputRange!L)
    {
        static if (isInputRange!R)
        {
            // Both lhs and rhs are ranges.
            for (;; lhs.popFront(), rhs.popFront())
            {
                if (lhs.empty) return rhs.empty;
                if (rhs.empty) return lhs.empty;
                if (!matchDigits(lhs.front, rhs.front, significantDigits))
                    return false;
            }
        }
        else
        {
            // lhs is a range, rhs is a number.
            for (; !lhs.empty; lhs.popFront)
            {
                if (!matchDigits(lhs.front, rhs, significantDigits))
                    return false;
            }
            return true;
        }
    }
    else
    {
        // lhs is a number, rhs is a range.
        return matchDigits(rhs, lhs, significantDigits);
    }
}


unittest
{
    check(matchDigits([0.1234566, 0.1234567, 0.1234568], 0.1234567, 6));
    check(matchDigits(0.1234567, [0.1234566, 0.1234567, 0.1234568], 6));
    check(matchDigits([0.1234566, 0.1234567, 0.1234568], 
        [0.1234567, 0.1234567, 0.1234567], 6));
}




/** Replaces real numbers that are close to zero by exactly zero.
    ---
    assert (chop(1.0) == 1.0);
    assert (chop(1e-20) == 0.0);
    ---
*/
Real chop(Real)(Real x, real threshold = 1e-10L) pure nothrow
    if (isFloatingPoint!(Real))
{
    if (fabs(x) < threshold) return 0.0;
    return x;
}

unittest
{
    check (chop(1.0) == 1.0);
    check (chop(1e-20) == 0.0);
}


/** Replaces all numbers in the given array that are close to zero
    by exactly zero. To chop the array in-place, pass the same array
    as both x and buffer.
    ---
    double[] a = [1.0, 1e-20, 2.0];
    double[] x = [1.0, 0.0, 2.0];
    auto b = chop(a);
    assert (b == x);
    chop(a, 1e-12L, a);
    assert (a == x);
    ---
*/
Real[] chop(Real) (Real[] x, real threshold = 1e-10L, Real[] buffer=null)
    nothrow
    if (isFloatingPoint!Real)
{
    if (buffer.length < x.length) buffer.length = x.length;

    foreach (i; 0 .. x.length)
        if (fabs(x[i]) < threshold)  buffer[i] = 0.0;
        else  buffer[i] = x[i];

    return buffer;
}

unittest
{
    real[] a = [1.0L, real.min_normal, 2.0];
    real[] x = [1.0L, 0.0, 2.0];

    auto b = chop(a);
    check (b == x);

    chop(a, 1e-12L, a);
    check (a == x);
}
