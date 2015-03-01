/** Unit testing tools.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009–2015, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.core.testing;

import scid.core.meta: Zero;
import scid.types;


/** This function is used to check a computed value along with
    its error estimate.  It evaluates to the following:
    ---
        abs(result-expected) <= absError
    &&  absError <= max(abs(result*relAccuracy), absAccuracy)
    ---
*/
bool isAccurate(T)(Result!T result, T expected, T relAccuracy,
    T absAccuracy=Zero!T)
{
    return isAccurate(result.value, result.error, expected, relAccuracy,
        absAccuracy);
}

/// ditto
bool isAccurate(T)(T result, T absError, T expected, T relAccuracy,
    T absAccuracy=Zero!T)
in
{
    assert (absError >= 0.0);
    assert (relAccuracy >= 0.0);
    assert (absAccuracy >= 0.0);
}
body
{
    import std.algorithm: max;
    import std.math: abs;
    return abs(result-expected) <= absError
        &&  absError <= max(abs(result*relAccuracy), absAccuracy);
}

unittest
{
    assert (isAccurate(2.0000001, 0.000001, 2.0, 0.000001));
}
