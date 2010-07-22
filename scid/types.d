/** Various useful types.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.types;


import std.conv;
import std.string: format;
import std.math;

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


    Result opUnary(string op)()  if (op == "-")
    {
        mixin ("return Result("~op~"value, error);");
    }

    Result opUnary(string op)()  if (op == "+")
    {
        return this;
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
    Result opBinary(string op)(Result rhs)  if (op == "+" || op == "-")
    {
        auto lhs = this;
        return lhs.opOpAssign!op(rhs);;
    }

    /// ditto
    Result opOpAssign(string op)(Result rhs)  if (op == "+" || op == "-")
    {
        mixin ("return Result(value"~op~"=rhs.value, error+=rhs.error);");
    }


    static if (is (V==E))
    {
        /// ditto
        Result opBinary(string op)(Result rhs)  if (op == "*" || op == "/")
        {
            auto lhs = this;
            return lhs.opOpAssign!op(rhs);
        }

        /// ditto
        Result opOpAssign(string op)(Result rhs)  if (op == "*")
        {
            return Result(
                value*=rhs.value,
                error = abs(value*rhs.error) + abs(rhs.value*error)
            );
        }

        /// ditto
        Result opOpAssign(string op)(Result rhs)  if (op == "/")
        {
            V inv = 1/rhs.value;
            return Result(
                value*=inv,
                error = (error + abs(rhs.error*value*inv))*abs(inv)
            );
        }
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
