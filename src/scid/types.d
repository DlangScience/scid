/** Various useful types.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.types;


import std.conv;
import std.string: format;
import std.math;




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


    /** Negate the value. */
    Result opNeg()
    {
        return Result(-value, error);
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
    Result opAdd(Result r)
    {
        return Result(value+r.value, error+r.error);
    }

    /// ditto
    Result opAddAssign(Result r)
    {
        return Result(value+=r.value, error+=r.error);
    }

    /// ditto
    Result opSub(Result r)
    {
        return Result(value-r.value, error+r.error);
    }

    /// ditto
    Result opSubAssign(Result r)
    {
        return Result(value-=r.value, error+=r.error);
    }

    static if (is (V==E))
    {
        /// ditto
        Result opMul(Result r)
        {
            return Result(
                value*r.value,
                abs(value*r.error) + abs(r.value*error)
            );
        }

        /// ditto
        Result opMulAssign(Result r)
        {
            return Result(
                value*=r.value,
                error = abs(value*r.error) + abs(r.value*error)
            );
        }

        /// ditto
        Result opDiv(Result r)
        {
            V inv = 1/r.value;
            return Result(
                value*inv,
                (error + abs(r.error*value*inv))*abs(inv)
            );
        }

        /// ditto
        Result opDivAssign(Result r)
        {
            V inv = 1/r.value;
            return Result(
                value*=inv,
                error = (error + abs(r.error*value*inv))*abs(inv)
            );
        }
    }
}


unittest
{
    alias Result!(real) rr;
    alias Result!(real, int) rri;
}
