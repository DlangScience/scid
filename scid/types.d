/** Various useful types.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.types;


import std.conv;
import std.format;
import std.math;
import std.string: format;
import std.traits;

import scid.core.testing;




/** Struct containing the result of a calculation, along with
    the absolute error in that calculation (x &plusmn; &delta;x).

    It is assumed, but not checked, that the error is nonnegative.
*/
struct Result(V, E=V)
{
    /// Result.
    V value;
    alias value this;

    /// Error.
    E error;

    invariant() { assert (error >= 0); }



    Result opUnary(string op)() const pure nothrow
        if (op == "-" || op == "+")
    {
        mixin ("return Result("~op~"value, error);");
    }


    /** Operators for Result-Result arithmetic.
        
        Note that these operations also perform error calculations, and
        are thus a bit slower than working directly with the value itself.
        It is assumed that the error is much smaller than the value, so
        terms of order O(&delta;x &delta;y) or O(&delta;x&sup2;) are ignored.

        Also note that multiplication and division are only possible when
        the value type V and the error type E are the same (as is the default).

        Through the magic of "alias this", Result!(V, E) is implicitly
        castable to the type V, and thus supports the same operations as V
        in addition to these.

    */
    Result opBinary(string op)(Result rhs) pure nothrow
        if (op == "+" || op == "-")
    {
        auto lhs = this;
        return lhs.opOpAssign!op(rhs);
    }

    /// ditto
    Result opOpAssign(string op)(Result rhs) pure nothrow
        if (op == "+" || op == "-")
    {
        mixin ("value "~op~"= rhs.value;");
        error += rhs.error;
        return this;
    }


    static if (is (V==E))
    {
        /// ditto
        Result opBinary(string op)(Result rhs) pure nothrow
            if (op == "*" || op == "/")
        {
            auto lhs = this;
            return lhs.opOpAssign!op(rhs);
        }

        /// ditto
        Result opOpAssign(string op)(Result rhs) pure nothrow
            if (op == "*")
        {
            value *= rhs.value;
            error = abs(value*rhs.error) + abs(rhs.value*error);
            return this;
        }

        /// ditto
        Result opOpAssign(string op)(Result rhs) pure nothrow
            if (op == "/")
        {
            V inv = 1/rhs.value;
            value *= inv;
            error = (error + abs(rhs.error*value*inv))*abs(inv);
            return this;
        }
    }


    /** Get a string representation of the result. */
    string toString
        (void delegate(const(char)[]) sink = null, string formatSpec = "%s")
        const
    {
        if (sink == null)
        {
            char[] buf;
            buf.reserve(100);
            toString((const(char)[] s) { buf ~= s; }, formatSpec);
            return cast(string) buf;
        }

        formattedWrite(sink, formatSpec, value);
        sink("\u00B1");
        formattedWrite(sink, formatSpec, error);
        return null;
    }
}


unittest
{
    alias Result!(real) rr;
    alias Result!(real, int) rri;
}


unittest
{
    auto r1 = Result!double(1.0, 0.1);
    auto r2 = Result!double(2.0, 0.2);

    check (+r1 == r1);
    
    auto ra = -r1;
    check (ra.value == -1.0);
    check (ra.error == 0.1);

    auto rb = r1 + r2;
    check (abs(rb.value - 3.0) < double.epsilon);
    check (abs(rb.error - 0.3) < double.epsilon);

    auto rc = r1 - r2;
    check (abs(rc.value + 1.0) < double.epsilon);
    check (abs(rc.error - 0.3) < double.epsilon);

    auto rd = r1 * r2;

    auto re = r1 / r2;
}


unittest
{
    auto r1 = Result!double(1.0, 0.1);
    check (r1.toString() == "1±0.1");

    auto r2 = Result!double(0.123456789, 0.00123456789);
    check (r2.toString(null, "%.8e") == "1.23456789e-01±1.23456789e-03");
}




/** An interval [a,b] along the real line, where either endpoint may
    be infinite.

    Examples:
    ---
    auto i1 = interval(1, 5);
    assert (i1.length == 4);
    
    auto i2 = interval(0, real.infinity);
    assert (i2.isInfinite);
    ---
*/
struct Interval(T) if (isSigned!T)
{
    /** The interval endpoints. */
    T a;
    T b;        ///ditto


pure: // TODO: Mark as nothrow as soon as DMD bug 5191 is fixed (DMD 2.051)


    /** The length of the interval, defined as b-a. */
    @property T length() @safe const { return b - a; }


    /** Determine whether the interval is infinite.  This is true if:
        $(UL
            $(LI a is infinite, b is finite)
            $(LI a is finite, b is infinite)
            $(LI a and b are infinite, but have opposite sign.
        )
        If T is an integer type, this is always false.
    */
    @property bool isInfinite() @trusted const
    {
        static if (isFloatingPoint!T) return isInfinity(length);
        else return false;
    }


    /** Determine whether this is an ordered interval, i.e.
        whether a <= b.
    */
    @property bool isOrdered() @safe const { return a <= b; }


    /** If a > b, swap the endpoints. */
    void order() @safe
    {
        if (a > b)
        {
            immutable tmp = a;
            a = b;
            b = tmp;
        }
    }
}


///ditto
Interval!(CommonType!(T, U)) interval(T, U)(T a, U b)
    @safe pure //nothrow
{
    return typeof(return)(a, b);
}


unittest
{
    auto fin = interval(1, 4);
    check(fin.a == 1);
    check(fin.b == 4);
    check(fin.length == 3);
    check(!fin.isInfinite);
    check(fin.isOrdered);
    fin.a = 9;
    check(!fin.isOrdered);
    fin.order();
    check(fin.a == 4 && fin.b == 9);
    check(fin.isOrdered);

    auto uin = interval(0.0, real.infinity);
    check(uin.isInfinite);
    check(uin.isOrdered);
}
