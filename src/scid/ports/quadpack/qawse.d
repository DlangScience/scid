/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawse;


import std.algorithm: max;
import std.math;

import scid.core.fortran;
import scid.ports.quadpack.qc25s;
import scid.ports.quadpack.qmomo;
import scid.ports.quadpack.qpsrt;




///
void qawse(Real, Func)(Func f, Real a, Real b, Real alfa, Real beta,
    int integr, Real epsabs, Real epsrel, int limit, out Real result,
    out Real abserr, out int neval, out int ier, Real* alist_,
    Real* blist_, Real* rlist_, Real* elist_, int* iord_, out int last)
{
//***begin prologue  dqawse
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1
//***keywords  automatic integrator, special-purpose,
//             algebraico-logarithmic end point singularities,
//             clenshaw-curtis method
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral i = integral of f*w over (a,b),
//            (where w shows a singular behaviour at the end points,
//            see parameter integr).
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i)).
//***description
//
//        integration of functions having algebraico-logarithmic
//        end point singularities
//        standard fortran subroutine
//        double precision version
//
//        parameters
//         on entry
//            f      - double precision
//                     function subprogram defining the integrand
//                     function f(x). the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            a      - double precision
//                     lower limit of integration
//
//            b      - double precision
//                     upper limit of integration, b.gt.a
//                     if b.le.a, the routine will end with ier = 6.
//
//            alfa   - double precision
//                     parameter in the weight function, alfa.gt.(-1)
//                     if alfa.le.(-1), the routine will end with
//                     ier = 6.
//
//            beta   - double precision
//                     parameter in the weight function, beta.gt.(-1)
//                     if beta.le.(-1), the routine will end with
//                     ier = 6.
//
//            integr - integer
//                     indicates which weight function is to be used
//                     = 1  (x-a)**alfa*(b-x)**beta
//                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
//                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
//                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
//                     if integr.lt.1 or integr.gt.4, the routine
//                     will end with ier = 6.
//
//            epsabs - double precision
//                     absolute accuracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
//            limit  - integer
//                     gives an upper bound on the number of subintervals
//                     in the partition of (a,b), limit.ge.2
//                     if limit.lt.2, the routine will end with ier = 6.
//
//         on return
//            result - double precision
//                     approximation to the integral
//
//            abserr - double precision
//                     estimate of the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            neval  - integer
//                     number of integrand evaluations
//
//            ier    - integer
//                     ier = 0 normal and reliable termination of the
//                             routine. it is assumed that the requested
//                             accuracy has been achieved.
//                     ier.gt.0 abnormal termination of the routine
//                             the estimates for the integral and error
//                             are less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                         = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit. however, if this yields no
//                             improvement, it is advised to analyze the
//                             integrand in order to determine the
//                             integration difficulties which prevent the
//                             requested tolerance from being achieved.
//                             in case of a jump discontinuity or a local
//                             singularity of algebraico-logarithmic type
//                             at one or more interior points of the
//                             integration range, one should proceed by
//                             splitting up the interval at these
//                             points and calling the integrator on the
//                             subranges.
//                         = 2 the occurrence of roundoff error is
//                             detected, which prevents the requested
//                             tolerance from being achieved.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 6 the input is invalid, because
//                             b.le.a or alfa.le.(-1) or beta.le.(-1), or
//                             integr.lt.1 or integr.gt.4, or
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                             or limit.lt.2.
//                             result, abserr, neval, rlist(1), elist(1),
//                             iord(1) and last are set to zero. alist(1)
//                             and blist(1) are set to a and b
//                             respectively.
//
//            alist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the left
//                     end points of the subintervals in the partition
//                     of the given integration range (a,b)
//
//            blist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the right
//                     end points of the subintervals in the partition
//                     of the given integration range (a,b)
//
//            rlist  - double precision
//                     vector of dimension at least limit,the first
//                      last  elements of which are the integral
//                     approximations on the subintervals
//
//            elist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the moduli of the
//                     absolute error estimates on the subintervals
//
//            iord   - integer
//                     vector of dimension at least limit, the first k
//                     of which are pointers to the error
//                     estimates over the subintervals, so that
//                     elist(iord(1)), ..., elist(iord(k)) with k = last
//                     if last.le.(limit/2+2), and k = limit+1-last
//                     otherwise form a decreasing sequence
//
//            last   - integer
//                     number of subintervals actually produced in
//                     the subdivision process
//
//***references  (none)
//***routines called  d1mach,dqc25s,dqmomo,dqpsrt
//***end prologue  dqawse
//
      Real area,area1,area12,area2,a1,
        a2,b1,b2,centre,epmach,
        errbnd,errmax,error1,erro12,error2,errsum,
        resas1,resas2,uflow;
      Real[25] ri_, rj_, rh_, rg_; 
      int iroff1,iroff2,k,maxerr,nev,nrmax;
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto elist = dimension(elist_, limit);
      auto rlist = dimension(rlist_, limit);
      auto iord = dimension(iord_, limit);
      auto rg = dimension(rg_.ptr, limit);
      auto rh = dimension(rh_.ptr, limit);
      auto ri = dimension(ri_.ptr, limit);
      auto rj = dimension(rj_.ptr, limit);
//
//            list of major variables
//            -----------------------
//
//           alist     - list of left end points of all subintervals
//                       considered up to now
//           blist     - list of right end points of all subintervals
//                       considered up to now
//           rlist(i)  - approximation to the integral over
//                       (alist(i),blist(i))
//           elist(i)  - error estimate applying to rlist(i)
//           maxerr    - pointer to the interval with largest
//                       error estimate
//           errmax    - elist(maxerr)
//           area      - sum of the integrals over the subintervals
//           errsum    - sum of the errors over the subintervals
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//           *****1    - variable for the left subinterval
//           *****2    - variable for the right subinterval
//           last      - index for subdivision
//
//
//            machine dependent constants
//            ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqawse
      epmach = Real.epsilon;
      uflow = Real.min;
//
//           test on validity of parameters
//           ------------------------------
//
      ier = 6;
      neval = 0;
      last = 0;
      rlist[1] = 0.0;
      elist[1] = 0.0;
      iord[1] = 0;
      result = 0.0;
      abserr = 0.0;
      if(b <= a || (epsabs == 0.0 && 
        epsrel < max(0.5e2*epmach,0.5e-28)) || alfa <= (-0.1e1)
        || beta <= (-0.1e1) || integr < 1 || integr > 4 || 
        limit < 2) goto l999;
      ier = 0;
//
//           compute the modified chebyshev moments.
//
      qmomo!Real(alfa,beta,ri.ptr,rj.ptr,rg.ptr,rh.ptr,integr);
//
//           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
//
      centre = 0.5*(b+a);
      qc25s!(Real, Func)(f,a,b,a,centre,alfa,beta,ri.ptr,rj.ptr,rg.ptr,rh.ptr,area1,
        error1,resas1,integr,nev);
      neval = nev;
      qc25s!(Real, Func)(f,a,b,centre,b,alfa,beta,ri.ptr,rj.ptr,rg.ptr,rh.ptr,area2,
        error2,resas2,integr,nev);
      last = 2;
      neval = neval+nev;
      result = area1+area2;
      abserr = error1+error2;
//
//           test on accuracy.
//
      errbnd = max(epsabs,epsrel*fabs(result));
//
//           initialization
//           --------------
//
      if(error2 > error1) goto l10;
      alist[1] = a;
      alist[2] = centre;
      blist[1] = centre;
      blist[2] = b;
      rlist[1] = area1;
      rlist[2] = area2;
      elist[1] = error1;
      elist[2] = error2;
      goto l20;
 l10: alist[1] = centre;
      alist[2] = a;
      blist[1] = b;
      blist[2] = centre;
      rlist[1] = area2;
      rlist[2] = area1;
      elist[1] = error2;
      elist[2] = error1;
 l20: iord[1] = 1;
      iord[2] = 2;
      if(limit == 2) ier = 1;
      if(abserr <= errbnd || ier == 1) goto l999;
      errmax = elist[1];
      maxerr = 1;
      nrmax = 1;
      area = result;
      errsum = abserr;
      iroff1 = 0;
      iroff2 = 0;
//
//            main do-loop
//            ------------
//
      for (last=3; last<=limit; last++) { //do 60 last = 3,limit
//
//           bisect the subinterval with largest error estimate.
//
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
//
        qc25s!(Real,Func)(f,a,b,a1,b1,alfa,beta,ri.ptr,rj.ptr,rg.ptr,rh.ptr,area1,
          error1,resas1,integr,nev);
        neval = neval+nev;
        qc25s!(Real,Func)(f,a,b,a2,b2,alfa,beta,ri.ptr,rj.ptr,rg.ptr,rh.ptr,area2,
          error2,resas2,integr,nev);
        neval = neval+nev;
//
//           improve previous approximations integral and error
//           and test for accuracy.
//
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if(a == a1 || b == b2) goto l30;
        if(resas1 == error1 || resas2 == error2) goto l30;
//
//           test for roundoff error.
//
        if(fabs(rlist[maxerr]-area12) < 0.1e-4*fabs(area12)
          && erro12 >= 0.99*errmax) iroff1 = iroff1+1;
        if(last > 10 && erro12 > errmax) iroff2 = iroff2+1;
 l30:   rlist[maxerr] = area1;
        rlist[last] = area2;
//
//           test on accuracy.
//
        errbnd = max(epsabs,epsrel*fabs(area));
        if(errsum <= errbnd) goto l35;
//
//           set error flag in the case that the number of interval
//           bisections exceeds limit.
//
        if(last == limit) ier = 1;
//
//
//           set error flag in the case of roundoff error.
//
        if(iroff1 >= 6 || iroff2 >= 20) ier = 2;
//
//           set error flag in the case of bad integrand behaviour
//           at interior points of integration range.
//
        if(max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*epmach)*
          (fabs(a2)+0.1e4*uflow)) ier = 3;
//
//           append the newly-created intervals to the list.
//
 l35:   if(error2 > error1) goto l40;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l50;
 l40:   alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
//
//           call subroutine dqpsrt to maintain the descending ordering
//           in the list of error estimates and select the subinterval
//           with largest error estimate (to be bisected next).
//
 l50:   qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
// ***jump out of do-loop
        if (ier != 0 || errsum <= errbnd) goto l70;
 l60: ;}
//
//           compute final result.
//           ---------------------
//
 l70: result = 0.0;
      for (k=1; k<=last; k++) { //do 80 k=1,last
        result = result+rlist[k];
 l80: ;}
      abserr = errsum;
l999: return;
}

unittest
{
    alias qawse!(float, float delegate(float)) fqawse;
    alias qawse!(double, double delegate(double)) dqawse;
    alias qawse!(double, double function(double)) dfqawse;
    alias qawse!(real, real delegate(real)) rqawse;
}
