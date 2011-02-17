/** Stuff that is useful when porting FORTRAN code to D.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.core.fortran;


import std.conv;
import std.traits: Unqual;

import scid.core.testing;




/** Wrap a one- or two-dimensional array around the given pointer.
    Meant as a substitute for FORTRAN's dimension statement.
*/
FortranArray!T dimension(T)(T* ptr, size_t len)
{
    FortranArray!T a;
    a._length = len;
    a._ptr = ptr;
    a._ptrm1 = ptr - 1;
    return a;
}


/// ditto
FortranArray2D!T dimension(T)(T* ptr, size_t rows, size_t cols)
{
    FortranArray2D!T a;
    a._length = rows*cols;
    a._ptr = ptr;
    a._ptrm1 = ptr - 1;
    a._rows = rows;
    a._cols = cols;
    return a;
}




/** A simple, lightweight, one-dimensional base-1 array. */
struct FortranArray(T)
{
private:
    // By putting _length and _ptr first, this array is binary compatible
    // with built-in D arrays.
    size_t _length;
    T* _ptr;
    T* _ptrm1;

    void boundsCheck(size_t i, string file, int line) const
    {
        assert (i>0  &&  i<=_length,
            file~":"~to!string(line)~":index out of bounds");
    }


public:
    /// The length of the array.
    @property size_t length() const { return _length; }


    /// A pointer to the first element of the array.
    @property T* ptr() { return _ptr; }


    static if (is(T == const) || is(T == immutable))
    {
        Unqual!T opIndex
            (string file = __FILE__, int line = __LINE__)
            (size_t i)
            const
        {
            boundsCheck(i, file, line);
            return _ptrm1[i];
        }
    }
    else
    {
        // Workaround for DMD bug 2460
        template opIndex(string file = __FILE__, int line = __LINE__) {
        ref T opIndex
            (size_t i)
        {
            boundsCheck(i, file, line);
            return _ptrm1[i];
        }}

        T opIndexAssign
            (string file = __FILE__, int line = __LINE__)
            (T value, size_t i)
        {
            boundsCheck(i, file, line);
            return (_ptrm1[i] = value);
        }

        T opIndexOpAssign
            (string op, string file = __FILE__, int line = __LINE__)
            (T value, size_t i)
        {
            boundsCheck(i, file, line);
            mixin("return _ptrm1[i] "~op~"= value;");
        }
    }
}


unittest
{
    int[] a = [1, 2, 3];
    auto b = dimension(a.ptr, 3);

    // Basic functionality
    check (b[3] == a[2]);
    b[3] = 5;
    check (b[3] == 5);
    b[3] -= 2;
    check (b[3] == 3);

    check (b.ptr[2] == 3);

    // Binary compatible with D arrays.
    int[] d = *(cast(int[]*) &b);
    check (d == a);

    // Works with const and immutable
    immutable int[] ia = [1,2,3];
    auto ib = dimension(ia.ptr, 3);
    check (ib[3] == 3);
    check (!__traits(compiles, { ib[3] = 5; }));

}




/** A simple, lightweight, two-dimensional base-1 array. */
struct FortranArray2D(T)
{
private:
    // By putting length and ptr first, this array is binary compatible
    // with built-in D arrays.
    size_t _length;
    T* _ptr;
    T* _ptrm1;
    size_t _rows;
    size_t _cols;

    void boundsCheck(size_t i, size_t j, string file, int line) const
    {
        assert (i>0  &&  i<=_rows,
            file~":"~to!string(line)~":first index out of bounds");
        assert (j>0  &&  j<=_cols,
            file~":"~to!string(line)~":second index out of bounds");
    }


public:
    /// The number of rows and columns in the array.
    @property size_t rows() const { return _rows; }
    @property size_t cols() const { return _cols; }     /// ditto


    /// The number of elements in the array.
    @property size_t length() const { return _length; }


    /// A pointer to the first element of the array.
    @property T* ptr() { return _ptr; }


    static if (is(T == const) || is(T == immutable))
    {
        Unqual!T opIndex
            (string file = __FILE__, int line = __LINE__)
            (size_t i, size_t j)
            const
        {
            boundsCheck(i, j, file, line);
            return _ptrm1[i + (j-1)*_rows];
        }
    }
    else
    {
        // Workaround for DMD bug 2460
        template opIndex(string file = __FILE__, int line = __LINE__) {
        ref T opIndex
            (size_t i, size_t j)
        {
            boundsCheck(i, j, file, line);
            return _ptrm1[i + (j-1)*_rows];
        }}

        T opIndexAssign
            (string file = __FILE__, int line = __LINE__)
            (T value, size_t i, size_t j)
        {
            boundsCheck(i, j, file, line);
            return (_ptrm1[i + (j-1)*_rows] = value);
        }

        T opIndexOpAssign
            (string op, string file = __FILE__, int line = __LINE__)
            (T value, size_t i, size_t j)
        {
            boundsCheck(i, j, file, line);
            mixin("return _ptrm1[i + (j-1)*_rows] "~op~"= value;");
        }
    }
}

unittest
{
    int[] a = [1, 2, 3, 4, 5, 6];
    auto b = dimension(a.ptr, 3, 2);

    // Basic functionality
    check (b[2,2] == a[4]);
    b[2,2] = 5;
    check (b[2,2] == 5);
    b[2,2] -= 2;
    check (b[2,2] == 3);

    check (b.ptr[4] == 3);

    // Binary compatible with D arrays.
    int[] d = *(cast(int[]*) &b);
    check (d == a);

    // Works with const and immutable
    const int[] ia = [1, 2, 3, 4, 5, 6];
    auto ib = dimension(ia.ptr, 3, 2);
    check (ib[2,2] == ia[4]);
    check (!__traits(compiles, { ib[2,2] = 5; }));
}
