/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qagie;


import std.algorithm: max, min;
import std.math: fabs;

import scid.common.fortran;
import scid.common.meta: Zero, One;
import scid.ports.quadpack.qk15i;
import scid.ports.quadpack.qpsrt;
import scid.ports.quadpack.qelg;




///
void qagie(Real, Func)(Func f, Real bound, int inf, Real epsabs,
    Real epsrel, int limit, out Real result, out Real abserr,
    out int neval, out int ier, Real* alist_, Real* blist_,
    Real* rlist_, Real* elist_, int* iord_, out int last)
{
//***begin prologue  dqagie
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a3a1,h2a4a1
//***keywords  automatic integrator, infinite intervals,
//             general-purpose, transformation, extrapolation,
//             globally adaptive
//***author  piessens,robert,appl. math & progr. div - k.u.leuven
//           de doncker,elise,appl. math & progr. div - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            integral   i = integral of f over (bound,+infinity)
//            or i = integral of f over (-infinity,bound)
//            or i = integral of f over (-infinity,+infinity),
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i))
//***description
//
// integration over infinite intervals
// standard fortran subroutine
//
//            f      - double precision
//                     function subprogram defining the integrand
//                     function f(x). the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            bound  - double precision
//                     finite bound of integration range
//                     (has no meaning if interval is doubly-infinite)
//
//            inf    - double precision
//                     indicating the kind of integration range involved
//                     inf = 1 corresponds to  (bound,+infinity),
//                     inf = -1            to  (-infinity,bound),
//                     inf = 2             to (-infinity,+infinity).
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
//                     in the partition of (a,b), limit.ge.1
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
//                   - ier.gt.0 abnormal termination of the routine. the
//                             estimates for result and error are less
//                             reliable. it is assumed that the requested
//                             accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit (and taking the according dimension
//                             adjustments into account). however,if
//                             this yields no improvement it is advised
//                             to analyze the integrand in order to
//                             determine the integration difficulties.
//                             if the position of a local difficulty can
//                             be determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling the
//                             integrator on the subranges. if possible,
//                             an appropriate special-purpose integrator
//                             should be used, which is designed for
//                             handling the type of difficulty involved.
//                         = 2 the occurrence of roundoff error is
//                             detected, which prevents the requested
//                             tolerance from being achieved.
//                             the error may be under-estimated.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 4 the algorithm does not converge.
//                             roundoff error is detected in the
//                             extrapolation table.
//                             it is assumed that the requested tolerance
//                             cannot be achieved, and that the returned
//                             result is the best which can be obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted that
//                             divergence can occur with any other value
//                             of ier.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                             result, abserr, neval, last, rlist(1),
//                             elist(1) and iord(1) are set to zero.
//                             alist(1) and blist(1) are set to 0
//                             and 1 respectively.
//
//            alist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the left
//                     end points of the subintervals in the partition
//                     of the transformed integration range (0,1).
//
//            blist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the right
//                     end points of the subintervals in the partition
//                     of the transformed integration range (0,1).
//
//            rlist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the integral
//                     approximations on the subintervals
//
//            elist  - double precision
//                     vector of dimension at least limit,  the first
//                     last elements of which are the moduli of the
//                     absolute error estimates on the subintervals
//
//            iord   - integer
//                     vector of dimension limit, the first k
//                     elements of which are pointers to the
//                     error estimates over the subintervals,
//                     such that elist(iord(1)), ..., elist(iord(k))
//                     form a decreasing sequence, with k = last
//                     if last.le.(limit/2+2), and k = limit+1-last
//                     otherwise
//
//            last   - integer
//                     number of subintervals actually produced
//                     in the subdivision process
//
//***references  (none)
//***routines called  d1mach,dqelg,dqk15i,dqpsrt
//***end prologue  dqagie
      Real abseps,area,area1,area12,area2,a1,
        a2,boun,b1,b2,correc,defabs,defab1,defab2,
        dres,epmach,erlarg,erlast,
        errbnd,errmax,error1,error2,erro12,errsum,ertest,oflow,resabs,
        reseps,small,uflow;
      Real[3] res3la_;
      Real[52] rlist2_;
      int id=1,ierro=1,iroff1=1,iroff2=1,iroff3=1,jupbnd=1,k=1,ksgn=1,
        ktmin=1,maxerr=1,nres=1,nrmax=1,numrl2=1;
      bool extrap,noext;
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto elist = dimension(elist_, limit);
      auto iord = dimension(iord_, limit);
      auto res3la = dimension(res3la_.ptr, 3);
      auto rlist = dimension(rlist_, limit);
      auto rlist2 = dimension(rlist2_.ptr, 52);
//
//            the dimension of rlist2 is determined by the value of
//            limexp in subroutine dqelg.
//
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
//           rlist2    - array of dimension at least (limexp+2),
//                       containing the part of the epsilon table
//                       wich is still needed for further computations
//           elist(i)  - error estimate applying to rlist(i)
//           maxerr    - pointer to the interval with largest error
//                       estimate
//           errmax    - elist(maxerr)
//           erlast    - error on the interval currently subdivided
//                       (before that subdivision has taken place)
//           area      - sum of the integrals over the subintervals
//           errsum    - sum of the errors over the subintervals
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//           *****1    - variable for the left subinterval
//           *****2    - variable for the right subinterval
//           last      - index for subdivision
//           nres      - number of calls to the extrapolation routine
//           numrl2    - number of elements currently in rlist2. if an
//                       appropriate approximation to the compounded
//                       integral has been obtained, it is put in
//                       rlist2(numrl2) after numrl2 has been increased
//                       by one.
//           small     - length of the smallest interval considered up
//                       to now, multiplied by 1.5
//           erlarg    - sum of the errors over the intervals larger
//                       than the smallest interval considered up to now
//           extrap    - logical variable denoting that the routine
//                       is attempting to perform extrapolation. i.e.
//                       before subdividing the smallest interval we
//                       try to decrease the value of erlarg.
//           noext     - logical variable denoting that extrapolation
//                       is no longer allowed (true-value)
//
//            machine dependent constants
//            ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//           oflow is the largest positive magnitude.
//
//***first executable statement  dqagie
       epmach = Real.epsilon;
//
//           test on validity of parameters
//           -----------------------------
//
      ier = 0;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      alist[1] = 0.0;
      blist[1] = 0.1e1;
      rlist[1] = 0.0;
      elist[1] = 0.0;
      iord[1] = 0;
      if(epsabs <= 0.0 && epsrel < max(0.5e2*epmach,0.5e-28))
        ier = 6;
       if(ier == 6) goto l999;
//
//
//           first approximation to the integral
//           -----------------------------------
//
//           determine the interval to be mapped onto (0,1).
//           if inf = 2 the integral is computed as i = i1+i2, where
//           i1 = integral of f over (-infinity,0),
//           i2 = integral of f over (0,+infinity).
//
      boun = bound;
      if(inf == 2) boun = 0.0;
      qk15i!(Real, Func)(f,boun,inf,0.0,1.0,result,abserr,
        defabs,resabs);
//
//           test on accuracy
//
      last = 1;
      rlist[1] = result;
      elist[1] = abserr;
      iord[1] = 1;
      dres = fabs(result);
      errbnd = max(epsabs,epsrel*dres);
      if(abserr <= 1.0e2*epmach*defabs && abserr > errbnd) ier = 2;
      if(limit == 1) ier = 1;
      if(ier != 0 || (abserr <= errbnd && abserr != resabs) || 
        abserr == 0.0) goto l130;
//
//           initialization
//           --------------
//
      uflow = Real.min_normal;
      oflow = Real.max;
      rlist2[1] = result;
      errmax = abserr;
      maxerr = 1;
      area = result;
      errsum = abserr;
      abserr = oflow;
      nrmax = 1;
      nres = 0;
      ktmin = 0;
      numrl2 = 2;
      extrap = false;
      noext = false;
      ierro = 0;
      iroff1 = 0;
      iroff2 = 0;
      iroff3 = 0;
      ksgn = -1;
      if(dres >= (0.1e1-0.5e2*epmach)*defabs) ksgn = 1;
//
//           main do-loop
//           ------------
//
      for (last=2; last<=limit; last++) { // end: 90
//
//           bisect the subinterval with nrmax-th largest error estimate.
//
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        qk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1);
        qk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2);
//
//           improve previous approximations to integral
//           and error and test for accuracy.
//
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if(defab1 == error1 || defab2 == error2)goto l15;
        if(fabs(rlist[maxerr]-area12) > 0.1e-4*fabs(area12)
         || erro12 < 0.99*errmax) goto l10;
        if(extrap) iroff2 = iroff2+1;
        if(!extrap) iroff1 = iroff1+1;
 l10:   if(last > 10 && erro12 > errmax) iroff3 = iroff3+1;
 l15:   rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
//
//           test for roundoff error and eventually set error flag.
//
        if(iroff1+iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if(iroff2 >= 5) ierro = 3;
//
//           set error flag in the case that the number of
//           subintervals equals limit.
//
        if(last == limit) ier = 1;
//
//           set error flag in the case of bad integrand behaviour
//           at some points of the integration range.
//
        if(max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*epmach)*
        (fabs(a2)+0.1e4*uflow)) ier = 4;
//
//           append the newly-created intervals to the list.
//
        if(error2 > error1) goto l20;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l30;
 l20:   alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
//
//           call subroutine dqpsrt to maintain the descending ordering
//           in the list of error estimates and select the subinterval
//           with nrmax-th largest error estimate (to be bisected next).
//
 l30:   qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
        if(errsum <= errbnd) goto l115;
        if(ier != 0) goto l100;
        if(last == 2) goto l80;
        if(noext) goto l90;
        erlarg = erlarg-erlast;
        if(fabs(b1-a1) > small) erlarg = erlarg+erro12;
        if(extrap) goto l40;
//
//           test whether the interval to be bisected next is the
//           smallest interval.
//
        if(fabs(blist[maxerr]-alist[maxerr]) > small) goto l90;
        extrap = true;
        nrmax = 2;
 l40:   if(ierro == 3 || erlarg <= ertest) goto l60;
//
//           the smallest interval has the largest error.
//           before bisecting decrease the sum of the errors over the
//           larger intervals (erlarg) and perform extrapolation.
//
        id = nrmax;
        jupbnd = last;
        if(last > (2+limit/2)) jupbnd = limit+3-last;
        for(k=id; k<=jupbnd; k++) { // end: 50
          maxerr = iord[nrmax];
          errmax = elist[maxerr];
          if(fabs(blist[maxerr]-alist[maxerr]) > small) goto l90;
          nrmax = nrmax+1;
 l50:   ;}
//
//           perform extrapolation.
//
 l60:   numrl2 = numrl2+1;
        rlist2[numrl2] = area;
        qelg!Real(numrl2,rlist2.ptr,reseps,abseps,res3la.ptr,nres);
        ktmin = ktmin+1;
        if(ktmin > 5 && abserr < 0.1e-2*errsum) ier = 5;
        if(abseps >= abserr) goto l70;
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel*fabs(reseps));
        if(abserr <= ertest) goto l100;
//
//            prepare bisection of the smallest interval.
//
 l70:   if(numrl2 == 1) noext = true;
        if(ier == 5) goto l100;
        maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = false;
        small = small*0.5;
        erlarg = errsum;
        goto l90;
 l80:   small = 0.375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[2] = area;
 l90: ;}
//
//           set final result and error estimate.
//           ------------------------------------
//
l100: if(abserr == oflow) goto l115;
      if((ier+ierro) == 0) goto l110;
      if(ierro == 3) abserr = abserr+correc;
      if(ier == 0) ier = 3;
      if(result != 0.0 && area != 0.0)goto l105;
      if(abserr > errsum)goto l115;
      if(area == 0.0) goto l130;
      goto l110;
l105: if(abserr/fabs(result) > errsum/fabs(area))goto l115;
//
//           test on divergence
//
l110: if(ksgn == (-1) && max(fabs(result),fabs(area)) <= 
       defabs*0.1e-1) goto l130;
      if(0.1e-1 > (result/area) || (result/area) > 0.1e3
        || errsum > fabs(area)) ier = 6;
      goto l130;
//
//           compute global integral sum.
//
l115: result = 0.0;
      for (k=1; k<=last; k++) { // end: 120
        result = result+rlist[k];
l120: ;}
      abserr = errsum;
l130: neval = 30*last-15;
      if(inf == 2) neval = 2*neval;
      if(ier > 2) ier=ier-1;
l999: return;
}


unittest
{
    alias qagie!(float, float delegate(float)) fqagie;
    alias qagie!(double, double delegate(double)) dqagie;
    alias qagie!(double, double function(double)) dfqagie;
    alias qagie!(real, real delegate(real)) rqagie;
}
