/** Various useful utilities that don't really fit anywhere else.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.util;


import std.algorithm;
import std.complex;
import std.math;
import std.range;
import std.traits;
import std.typetuple: allSatisfy;

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

    Note that numbers which are very close to zero should normally be
    compared using an absolute difference criterion.  This can be
    specified using an optional parameter, and is by default set to 1e-20.
    Set it to zero if you do not want to use the absolute difference at all.
*/
bool matchDigits(L, R)
    (L lhs, R rhs, uint significantDigits = 6, real maxAbsDiff = 1e-20)
    if ((isFloatingPoint!L || is(L T == Complex!T)) &&
        (isFloatingPoint!R || is(R U == Complex!U)))
in
{
    assert (maxAbsDiff >= 0);
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
            return matchDigits(lhs.re, rhs.re, significantDigits, maxAbsDiff)
                && matchDigits(lhs.im, rhs.im, significantDigits, maxAbsDiff);
        }
        else
        {
            // lhs is complex, rhs is real
            return matchDigits(lhs.re, rhs, significantDigits, maxAbsDiff)
                && matchDigits(lhs.im, 0.0, significantDigits, maxAbsDiff);
        }
    }
    else
    {
        static if (is(R U == Complex!U))
        {
            // lhs is real, rhs is complex
            return matchDigits(lhs, rhs.re, significantDigits, maxAbsDiff)
                && matchDigits(0.0, rhs.im, significantDigits, maxAbsDiff);
        }
        else
        {
            // lhs and rhs are real
            return
                feqrel!(CommonType!(L, R))(lhs, rhs) > significantDigits*LOG2T
             || abs(lhs - rhs) <= maxAbsDiff;
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
bool matchDigits(L, R)
    (L lhs, R rhs, uint significantDigits = 6, real maxAbsDiff = 1e-20)
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
                if (!matchDigits(lhs.front, rhs.front, significantDigits, maxAbsDiff))
                    return false;
            }
        }
        else
        {
            // lhs is a range, rhs is a number.
            for (; !lhs.empty; lhs.popFront)
            {
                if (!matchDigits(lhs.front, rhs, significantDigits, maxAbsDiff))
                    return false;
            }
            return true;
        }
    }
    else
    {
        // lhs is a number, rhs is a range.
        return matchDigits(rhs, lhs, significantDigits, maxAbsDiff);
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




/** Create a static array literal without any heap allocation.

    staticArray() automatically deduces its type from the arguments,
    while staticArrayOf() lets you specify the type explicitly.
*/
CommonType!(T)[T.length] staticArray(T...)(T elements)
    @safe pure nothrow
    if (!is(CommonType!T == void))
{
    // Inspired by code posted by David Simcha on the Phobos
    // developers' mailing list.
    typeof(return) a = void;
    foreach (i, e; elements)  a[i] = e;
    return a;
}


/// ditto
T[U.length] staticArrayOf(T, U...)(U elements)
    @safe pure nothrow
    if (allConvertibleTo!(T, U))
{
    typeof(return) a = void;
    foreach (i, e; elements)  a[i] = e;
    return a;
}


unittest
{
    auto a = staticArray(0.0, 1.0, 2.0);
    auto b = staticArrayOf!double(0.0F, 1.0, 2.0L);
    double[3] c = [0.0, 1.0, 2.0];
    check (a == c);
    check (b == c);
}




/** Returns a range that iterates over n equally-spaced floating-point
    numbers in the inclusive interval [a,b].

    This is similar to std.range.iota, except that it allows
    you to specify the number of steps it takes rather than the step size,
    and that the last point is exactly equal to b (unless n = 1, in
    which case a is the first and last point).  This makes it
    more useful than iota for iterating over floating-point numbers.

    Example:
    ---
    int i = 0;
    foreach (x; steps(0.0, 9.0, 10))
    {
        assert (x == i);
        ++i;
    }
    ---
*/
Steps!(CommonType!(T, U)) steps(T, U)(T a, U b, int n)
    if (isFloatingPoint!T && isFloatingPoint!U)
{
    return typeof(return)(a, b, n);
}


struct Steps(T) if (isFloatingPoint!T)
{
private:
    int _i, _resolution;
    T _delta, _start, _stop;

public:
    this(T start, T stop, int resolution)
    in { assert (resolution >= 0); }
    body {
        _start = start;
        _stop = stop;
        if (resolution > 1) _delta = (stop - start) / (resolution - 1);
        _i = resolution - 1;
        _resolution = resolution;
    }


    bool empty() { return _i < 0; }


    T front()
    {
        assert (!empty);
        if (_i == _resolution - 1) return _start;
        return _stop - _i*_delta;
    }


    void popFront() {
        assert (!empty);
        --_i;
    }
}


unittest
{
    auto s1 = steps(-5.0, 5.0, 11);
    check (equal(s1, [-5.0, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]));

    auto s2 = steps(0.0, 1.0, 1);
    check (equal(s2, [0.0]));

    // From the example
    int i = 0;
    foreach (x; steps(0.0, 9.0, 10))
    {
        check (x == i);
        ++i;
    }
}
