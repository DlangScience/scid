/** Various useful utilities.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.util;


import std.math;
import std.traits;

import scid.core.testing;
import scid.core.traits;




/** Replaces real numbers that are close to zero by exactly zero.
    ---
    assert (chop(1.0) == 1.0);
    assert (chop(1e-20) == 0.0);
    ---
*/
pure Real chop(Real)(Real x, real threshold = 1e-10L)
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
