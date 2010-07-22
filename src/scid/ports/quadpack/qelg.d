/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qelg;


import std.algorithm: min, max;
import std.math;

import scid.core.fortran;




///
void qelg(Real)(ref int n, Real* epstab_, ref Real result,
    ref Real abserr, Real* res3la_, ref int nres)
{
//***begin prologue  dqelg
//***refer to  dqagie,dqagoe,dqagpe,dqagse
//***routines called  d1mach
//***revision date  830518   (yymmdd)
//***keywords  epsilon algorithm, convergence acceleration,
//             extrapolation
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math & progr. div. - k.u.leuven
//***purpose  the routine determines the limit of a given sequence of
//            approximations, by means of the epsilon algorithm of
//            p.wynn. an estimate of the absolute error is also given.
//            the condensed epsilon table is computed. only those
//            elements needed for the computation of the next diagonal
//            are preserved.
//***description
//
//           epsilon algorithm
//           standard fortran subroutine
//           double precision version
//
//           parameters
//              n      - integer
//                       epstab(n) contains the new element in the
//                       first column of the epsilon table.
//
//              epstab - double precision
//                       vector of dimension 52 containing the elements
//                       of the two lower diagonals of the triangular
//                       epsilon table. the elements are numbered
//                       starting at the right-hand corner of the
//                       triangle.
//
//              result - double precision
//                       resulting approximation to the integral
//
//              abserr - double precision
//                       estimate of the absolute error computed from
//                       result and the 3 previous results
//
//              res3la - double precision
//                       vector of dimension 3 containing the last 3
//                       results
//
//              nres   - integer
//                       number of calls to the routine
//                       (should be zero at first call)
//
//***end prologue  dqelg
//
      Real delta1,delta2,delta3,
       epmach,epsinf,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
       oflow,res,ss,tol1,tol2,tol3;
      int i=1,ib=1,ib2=1,ie=1,indx=1,k1=1,k2=1,k3=1,limexp=1,newelm=1,num=1;
      auto epstab = dimension(epstab_, 52);
      auto res3la = dimension(res3la_, 3);
//
//           list of major variables
//           -----------------------
//
//           e0     - the 4 elements on which the computation of a new
//           e1       element in the epsilon table is based
//           e2
//           e3                 e0
//                        e3    e1    new
//                              e2
//           newelm - number of elements to be computed in the new
//                    diagonal
//           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
//           result - the element in the new diagonal with least value
//                    of error
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           oflow is the largest positive magnitude.
//           limexp is the maximum number of elements the epsilon
//           table can contain. if this number is reached, the upper
//           diagonal of the epsilon table is deleted.
//
//***first executable statement  dqelg
      epmach = Real.epsilon;
      oflow = Real.max;
      nres = nres+1;
      abserr = oflow;
      result = epstab[n];
      if (n < 3) goto l100;
      limexp = 50;
      epstab[n+2] = epstab[n];
      newelm = (n-1)/2;
      epstab[n] = oflow;
      num = n;
      k1 = n;
      for (i=1; i<=newelm; i++) { // end: l40
        k2 = k1-1;
        k3 = k1-2;
        res = epstab[k1+2];
        e0 = epstab[k3];
        e1 = epstab[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2-e1;
        err2 = fabs(delta2);
        tol2 = max(fabs(e2),e1abs)*epmach;
        delta3 = e1-e0;
        err3 = fabs(delta3);
        tol3 = max(e1abs,fabs(e0))*epmach;
        if (err2 > tol2 || err3 > tol3) goto l10;
//
//           if e0, e1 and e2 are equal to within machine
//           accuracy, convergence is assumed.
//           result = e2
//           abserr = abs(e1-e0)+abs(e2-e1)
//
        result = res;
        abserr = err2+err3;
// ***jump out of do-loop
        goto l100;
 l10:   e3 = epstab[k1];
        epstab[k1] = e1;
        delta1 = e1-e3;
        err1 = fabs(delta1);
        tol1 = max(e1abs,fabs(e3))*epmach;
//
//           if two elements are very close to each other, omit
//           a part of the table by adjusting the value of n
//
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) goto l20;
        ss = 0.1e1/delta1+0.1e1/delta2-0.1e1/delta3;
        epsinf = fabs(ss*e1);
//
//           test to detect irregular behaviour in the table, and
//           eventually omit a part of the table adjusting the value
//           of n.
//
        if(epsinf > 0.1e-3) goto l30;
 l20:   n = i+i-1;
// ***jump out of do-loop
        goto l50;
//
//           compute a new element and eventually adjust
//           the value of result.
//
 l30:   res = e1+0.1e1/ss;
        epstab[k1] = res;
        k1 = k1-2;
        error = err2+fabs(res-e2)+err3;
        if (error > abserr) goto l40;
        abserr = error;
        result = res;
 l40:;}
//
//           shift the table.
//
 l50: if (n == limexp) n = 2*(limexp/2)-1;
      ib = 1;
      if ((num/2)*2 == num) ib = 2;
      ie = newelm+1;
      for (i=1; i<=ie; i++) { // end: l60
        ib2 = ib+2;
        epstab[ib] = epstab[ib2];
        ib = ib2;
 l60:;}
      if (num == n) goto l80;
      indx = num-n+1;
      for (i=1; i<=n; i++) { // end: l70
        epstab[i]= epstab[indx];
        indx = indx+1;
 l70:;}
 l80: if (nres >= 4) goto l90;
      res3la[nres] = result;
      abserr = oflow;
      goto l100;
//
//           compute error estimate
//
 l90: abserr = fabs(result-res3la[3])+fabs(result-res3la[2])
        +fabs(result-res3la[1]);
      res3la[1] = res3la[2];
      res3la[2] = res3la[3];
      res3la[3] = result;
l100: abserr = max(abserr, 0.5e1*epmach*fabs(result));
      return;
}


unittest
{
    alias qelg!float fqelg;
    alias qelg!double dqelg;
    alias qelg!real rqelg;
}
