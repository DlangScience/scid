/** Unit testing tools.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009–2011, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.core.testing;

import std.algorithm: max;
import std.math;
import std.stdio;
import std.string;

import scid.core.meta: Zero;
import scid.types;




/** This function is supposed to be a drop-in replacement for assert()
    in unittests. Instead of halting the entire process when a test fails,
    it prints an error message to STDERR with the file and line number of
    the test that failed.
*/
void check
    (string file=__FILE__, int line=__LINE__)
    (bool test, lazy string msg = null)
{
    if (test)  { _checkReport.success++; return; }

    _checkReport.fail++;
    if (msg == null)
        stderr.writefln("check() failed: %s(%s)", file, line);
    else
        stderr.writefln("check() failed: %s(%s): %s", file, line, msg);
}




/** Returns a report on the tests performed by check(). */
CheckReport checkReport()  { return _checkReport; }

/** Data about the tests performed by check. */
struct CheckReport
{
    /** The number of succeeding check() calls. */
    int success = 0;

    /** The number of failing check() calls. */
    int fail = 0;

    /** The total number of check() calls. */
    @property int total() { return success + fail; }
}
// Global instance
private CheckReport _checkReport;




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
    return abs(result-expected) <= absError
        &&  absError <= max(abs(result*relAccuracy), absAccuracy);
}

unittest
{
    check (isAccurate(2.0000001, 0.000001, 2.0, 0.000001));
}
