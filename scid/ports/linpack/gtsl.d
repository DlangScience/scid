/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0 
*/
module scid.ports.linpack.gtsl;


import std.math: fabs;

import scid.core.fortran;




///
void gtsl(Real)(int n, Real* c_, Real* d_, Real* e_, Real* b_, out int info)
{
      auto c = dimension(c_, n);
      auto d = dimension(d_, n);
      auto e = dimension(e_, n);
      auto b = dimension(b_, n);
//
//     dgtsl given a general tridiagonal matrix and a right hand
//     side will find the solution.
//
//     on entry
//
//        n       integer
//                is the order of the tridiagonal matrix.
//
//        c       double precision(n)
//                is the subdiagonal of the tridiagonal matrix.
//                c(2) through c(n) should contain the subdiagonal.
//                on output c is destroyed.
//
//        d       double precision(n)
//                is the diagonal of the tridiagonal matrix.
//                on output d is destroyed.
//
//        e       double precision(n)
//                is the superdiagonal of the tridiagonal matrix.
//                e(1) through e(n-1) should contain the superdiagonal.
//                on output e is destroyed.
//
//        b       double precision(n)
//                is the right hand side vector.
//
//     on return
//
//        b       is the solution vector.
//
//        info    integer
//                = 0 normal value.
//                = k if the k-th element of the diagonal becomes
//                    exactly zero.  the subroutine returns when
//                    this is detected.
//
//     linpack. this version dated 08/14/78 .
//     jack dongarra, argonne national laboratory.
//
//     no externals
//     fortran dabs
//
//     internal variables
//
      int k=1,kb=1,kp1=1,nm1=1,nm2=1;
      Real t;
//     begin block permitting ...exits to 100
//
         info = 0;
         c[1] = d[1];
         nm1 = n - 1;
         if (nm1 < 1) goto l40;
            d[1] = e[1];
            e[1] = 0.0;
            e[n] = 0.0;
//
            for (k=1; k<=nm1; k++) { // end: 30
               kp1 = k + 1;
//
//              find the largest of the two rows
//
               if (fabs(c[kp1]) < fabs(c[k])) goto l10;
//
//                 interchange row
//
                  t = c[kp1];
                  c[kp1] = c[k];
                  c[k] = t;
                  t = d[kp1];
                  d[kp1] = d[k];
                  d[k] = t;
                  t = e[kp1];
                  e[kp1] = e[k];
                  e[k] = t;
                  t = b[kp1];
                  b[kp1] = b[k];
                  b[k] = t;
 l10:          ;
//
//              zero elements
//
               if (c[k] != 0.0) goto l20;
                  info = k;
//     ............exit
                  goto l100;
 l20:          ;
               t = -c[kp1]/c[k];
               c[kp1] = d[kp1] + t*d[k];
               d[kp1] = e[kp1] + t*e[k];
               e[kp1] = 0.0;
               b[kp1] = b[kp1] + t*b[k];
 l30:       ;}
 l40:    ;
         if (c[n] != 0.0) goto l50;
            info = n;
         goto l90;
 l50:    ;
//
//           back solve
//
            nm2 = n - 2;
            b[n] = b[n]/c[n];
            if (n == 1) goto l80;
               b[nm1] = (b[nm1] - d[nm1]*b[n])/c[nm1];
               if (nm2 < 1) goto l70;
               for (kb=1; kb<=nm2; kb++) { // end: 60
                  k = nm2 - kb + 1;
                  b[k] = (b[k] - d[k]*b[k+1] - e[k]*b[k+2])/c[k];
 l60:          ;}
 l70:          ;
 l80:       ;
 l90:    ;
l100: ;
//
      return;
}


unittest
{
    alias gtsl!float fgtsl;
    alias gtsl!double dgtsl;
    alias gtsl!real rgtsl;
}
