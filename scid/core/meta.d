/** Facilities for template metaprogramming.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.core.meta;


private import std.traits;
private import scid.core.traits;




/** Evaluates to the zero value for a given type.
    ---
    assert (Zero!creal == 0.0+0.0i);
    ---
*/
template Zero(T)
{
    static if (isFloatingPoint!T)
        enum T Zero = 0.0;
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




template RepeatString(string str, int n)
{
    static assert (n>0);
    static if (n==1)
        enum string RepeatString = str;
    else
        enum string RepeatString = str ~ RepeatString!(str, n-1);
}


version(unittest)
{
    static assert (RepeatString!("foo", 3) == "foofoofoo");
}



/** Replace all occurences of from in text with to. */
string replace(string text, string from, string to)
{
    char[] s = text.dup;
    //s[] = (cast(char[])text)[];

    int i=0;
    while (i<s.length)
    {
        int ipl = i + from.length;
        if (ipl > s.length)  break;

        if (s[i .. ipl] == from)
        {
            if (from.length == to.length)
                foreach (j; 0 .. from.length)  s[i+j] = to[j];
            else
                s = s[0 .. i]~to~s[ipl .. $];

            i = ipl;
        }
        else  i++;
    }
    return cast(immutable) s;
}

unittest
{
    enum string test = replace("foo bar baz", "bar", "bob");
    static assert (test == "foo bob baz");
}
