/** Facilities for template metaprogramming.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009–2011, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.core.meta;

import std.traits;
import scid.core.traits;
import std.complex: Complex;


/** Evaluates to the zero value for a given type.
    ---
    assert (Zero!creal == 0.0+0.0i);
    ---
*/
template Zero(T)
{
    static if (isFloatingPoint!T)
        enum T Zero = 0.0;
    else static if (is(T : Complex!E, E))
        enum T Zero = T(0, 0);
    else static if (isComplex!T)
        enum T Zero = 0.0 + 0.0i;
    else static if (isIntegral!T)
        enum T Zero = 0;
    else static assert (false, "Zero: Type has no obvious zero: "~T.stringof);
}

version(unittest)
{
    static assert (Zero!creal == 0.0+0.0i);
}




/** Evaluates to the unit value for a given type.
    ---
    assert (One!creal = 1.0+0.0i);
    ---
*/
template One(T)
{
    static if (isFloatingPoint!T)
        enum T One = 1.0;
    else static if (is(T : Complex!E, E))
        enum T One = T(1, 0);
    else static if (isComplex!T)
        enum T One = 1.0 + 0.0i;
    else static if (isIntegral!T)
        enum T One = 1;
    else static assert (false, "One: Type has no obvious unit element: "
        ~T.stringof);
}

version(unittest)
{
    static assert (One!creal == 1.0+0.0i);
}
