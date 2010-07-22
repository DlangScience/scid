/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0 
*/
module scid.ports.minpack.enorm;


import std.math;



/** Given an n-vector x, this function calculates the
    Euclidean norm of x.

    The Euclidean norm is computed by accumulating the sum of
    squares in three different sums. The sums of squares for the
    small and large components are scaled so that no overflows
    occur. Non-destructive underflows are permitted. Underflows
    and overflows do not occur in the computation of the unscaled
    sum of squarees for the intermediate components.
    The definitions of small, intermediate and large components
    depend on two constants, rdwarf and rgiant. The main
    restrictions on these constants are that rdwarf^2 not
    underflow and rgiant^2 not overflow. The constants
    given here are suitable for every known computer.

    Parameters:
        x = an input array.
*/
Real enorm(Real)(size_t n, Real* x)
{
    enum : Real
    {
        one = 1.0,
        zero = 0.0,
        rdwarf = 3.834e-20,
        rgiant = 1.304e19
    }

    Real s1 = zero;
    Real s2 = zero;
    Real s3 = zero;
    Real x1max = zero;
    Real x3max = zero;
    immutable Real floatn = n;
    immutable Real agiant = rgiant/floatn;

    Real xabs;
    for (size_t i=0; i<n; i++)
    {
        xabs = abs(x[i]);
        if (xabs <= rdwarf  ||  xabs >= agiant)
        {
            if (xabs > rdwarf)
            {
                // Sum for large components.
                if (xabs > x1max)
                {
                    s1 = one + s1*((x1max/xabs)^^2);
                    x1max = xabs;
                }
                else
                {
                    s1 += (xabs/x1max)^^2;
                }
            }
            else
            {
                // Sum for small components.
                if (xabs > x3max)
                {
                    s3 = one + s3*((x3max/xabs)^^2);
                    x3max = xabs;
                }
                else
                {
                    if (xabs != zero) s3 += (xabs/x3max)^^2;
                }
            }
        }
        else
        {
            // Sum for intermediate components.
            s2 += xabs*xabs;
        }
    }

    
    // Calculation of norm.
    Real result;
    if (s1 != zero)
    {
        result = x1max*sqrt(s1+(s2/x1max)/x1max);
    }
    else
    {
        if (s2 != zero)
        {
            if (s2 >= x3max)
                result = sqrt(s2*(one+(x3max/s2)*(x3max*s3)));
            else
                result = sqrt(x3max*((s2/x3max)+(x3max*s3)));
        }
        else
        {
            result = x3max*sqrt(s3);
        }
    }
    
    return result;
}


unittest
{
    double[] v = [ 1.0, 2.0, 3.0 ];
    double norm = enorm(3, v.ptr);
    assert (abs(norm-sqrt(14.0)) < 1e-10);
}

