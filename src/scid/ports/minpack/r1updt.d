/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0 
*/
module scid.ports.minpack.r1updt;


private import std.math: abs, sqrt;




/** Given an m by n lower trapezoidal matrix S, an m-vector u,
    and an n-vector v, the problem is to determine an
    orthogonal matrix Q such that
    ---
                  t
          (S + u v ) Q
    ---
    is again lower trapezoidal.

    This subroutine determines Q as the product of 2*(n - 1)
    transformations
    ---
          gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    ---
    where gv(i), gw(i) are Givens rotations in the (i,n) plane
    which eliminate elements in the i-th and n-th planes,
    respectively. Q itself is not accumulated, rather the
    information to recover the gv, gw rotations is returned.

    Params:

      m = a positive integer input variable set to the number
        of rows of S.

      n = a positive integer input variable set to the number
        of columns of S. n must not exceed m.

      s = an array of length ls. on input s must contain the lower
        trapezoidal matrix S stored by columns. on output s contains
        the lower trapezoidal matrix produced as described above.

      ls = a positive integer input variable not less than
        (n*(2*m-n+1))/2.

      u = an input array of length m which must contain the
        vector u.

      v = an array of length n. on input v must contain the vector
        v. on output v(i) contains the information necessary to
        recover the Givens rotation gv(i) described above.

      w = an output array of length m. w(i) contains information
        necessary to recover the Givens rotation gw(i) described
        above.

      sing = a logical output variable. sing is set true if any
        of the diagonal elements of the output s are zero. otherwise
        sing is set false.
*/
void r1updt(Real)(size_t m, size_t n, Real* s, size_t ls,
    Real* u, Real* v, Real* w, out bool sing)
{
    size_t i, j, jj, l, nmj, nm1;
    Real cos, cotan, sin, tan, tau, temp;
    
    enum : Real
    {
        one = 1.0,
        p5 = 0.5,
        p25 = 0.25,
        zero = 0.0,

        // giant is the largest magnitude
        giant = Real.max
    }


    // Initialize the diagonal element pointer.
    jj = (n*(2*m - n + 1))/2 - (m - n) - 1;


    // Move the nontrivial part of the last column of s into w.
    nm1 = n - 1;
    for (i=nm1, l=jj; i<m; i++, l++)
        w[i] = s[l];


    // Rotate the vector v into a multiple of the n-th unit vector
    // in such a way that a spike is introduced into w.
    if (n > 1)
    {
        for (nmj=2; nmj<=n; nmj++)
        {
            j = n - nmj;
            jj = jj - (m - j);
            w[j] = zero;

            if (v[j] != zero)
            {
                // Determine a Givens rotation which eliminates the
                // j-th element of v.
                if (abs(v[nm1]) < abs(v[j]))
                {
                    cotan = v[nm1]/v[j];
                    sin = p5/sqrt(p25+p25*cotan*cotan);
                    cos = sin*cotan;
                    tau = one;
                    if (abs(cos)*giant > one) tau = one/cos;
                }
                else
                {
                    tan = v[j]/v[nm1];
                    cos = p5/sqrt(p25+p25*tan*tan);
                    sin = cos*tan;
                    tau = sin;
                }


                // Apply the transformation to v and store the information
                // necessary to recover the Givens rotation.
                v[nm1] = sin*v[j] + cos*v[nm1];
                v[j] = tau;


                // Apply the transformation to S and extend the spike in w.
                for (i=j, l=jj; i<m; i++, l++)
                {
                    temp = cos*s[l] - sin*w[i];
                    w[i] = sin*s[l] + cos*w[i];
                    s[l] = temp;
                }
            }
        }
    }


    // Add the spike from the rank 1 update to w.
    for (i=0; i<m; i++)
        w[i] = w[i] + v[nm1]*u[i];


    // Eliminate the spike.
    sing = false;
    if (n > 1)
    {
        for (j=0; j<nm1; j++)
        {
            if (w[j] != zero)
            {
                // Determine a Givens rotation which eliminates the
                // j-th element of the spike.
                if (abs(s[jj]) < abs(w[j]))
                {
                    cotan = s[jj]/w[j];
                    sin = p5/sqrt(p25+p25*cotan*cotan);
                    cos = sin*cotan;
                    tau = one;
                    if (abs(cos)*giant > one) tau = one/cos;
                }
                else
                {
                    tan = w[j]/s[jj];
                    cos = p5/sqrt(p25+p25*tan*tan);
                    sin = cos*tan;
                    tau = sin;
                }


                // Apply the transformation to S and reduce the spike in w.
                for (i=j, l=jj; i<m; i++, l++)
                {
                    temp = cos*s[l] + sin*w[i];
                    w[i] = -sin*s[l] + cos*w[i];
                    s[l] = temp;
                }


                // Store the information necessary to recover the
                // Givens rotation.
                w[j] = tau;
            }


            // Test for zero diagonal elements in the output S.
            if (s[jj] == zero) sing = true;
            jj += m - j;
        }
    }


    // Move w back into the last column of the output s.
    for (i=nm1, l=jj; i<m; i++, l++)
        s[l] = w[i];
    if (s[jj] == zero) sing = true;
}


unittest { alias r1updt!(real) rr1updt; }
