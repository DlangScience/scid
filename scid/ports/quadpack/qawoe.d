// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawoe;


import std.algorithm: max;
import std.math: fabs;

import scid.core.fortran;
import scid.ports.quadpack.qc25f;
import scid.ports.quadpack.qpsrt;
import scid.ports.quadpack.qelg;




///
void qawoe(Real, Func)(Func f, Real a, Real b, Real omega, int integr,
    Real epsabs, Real epsrel, int limit, int icall, int maxp1,
    out Real result, out Real abserr, out int neval, out int ier,
    out int last, Real* alist_, Real* blist_, Real* rlist_, Real* elist_,
    int* iord_, int* nnlog_, ref int momcom, Real* chebmo_)
{
//***begin prologue  dqawoe
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1
//***keywords  automatic integrator, special-purpose,
//             integrand with oscillatory cos or sin factor,
//             clenshaw-curtis method, (end point) singularities,
//             extrapolation, globally adaptive
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral
//            i = integral of f(x)*w(x) over (a,b)
//            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i)).
//***description
//
//        computation of oscillatory integrals
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
//                     upper limit of integration
//
//            omega  - double precision
//                     parameter in the integrand weight function
//
//            integr - integer
//                     indicates which of the weight functions is to be
//                     used
//                     integr = 1      w(x) = cos(omega*x)
//                     integr = 2      w(x) = sin(omega*x)
//                     if integr.ne.1 and integr.ne.2, the routine
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
//                     gives an upper bound on the number of subdivisions
//                     in the partition of (a,b), limit.ge.1.
//
//            icall  - integer
//                     if dqawoe is to be used only once, icall must
//                     be set to 1.  assume that during this call, the
//                     chebyshev moments (for clenshaw-curtis integration
//                     of degree 24) have been computed for intervals of
//                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
//                     if icall.gt.1 this means that dqawoe has been
//                     called twice or more on intervals of the same
//                     length abs(b-a). the chebyshev moments already
//                     computed are then re-used in subsequent calls.
//                     if icall.lt.1, the routine will end with ier = 6.
//
//            maxp1  - integer
//                     gives an upper bound on the number of chebyshev
//                     moments which can be stored, i.e. for the
//                     intervals of lenghts abs(b-a)*2**(-l),
//                     l=0,1, ..., maxp1-2, maxp1.ge.1.
//                     if maxp1.lt.1, the routine will end with ier = 6.
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
//                             routine. it is assumed that the
//                             requested accuracy has been achieved.
//                   - ier.gt.0 abnormal termination of the routine.
//                             the estimates for integral and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit (and taking according dimension
//                             adjustments into account). however, if
//                             this yields no improvement it is advised
//                             to analyze the integrand, in order to
//                             determine the integration difficulties.
//                             if the position of a local difficulty can
//                             be determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling the
//                             integrator on the subranges. if possible,
//                             an appropriate special-purpose integrator
//                             should be used which is designed for
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
//                             it is presumed that the requested
//                             tolerance cannot be achieved due to
//                             roundoff in the extrapolation table,
//                             and that the returned result is the
//                             best which can be obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted that
//                             divergence can occur with any other value
//                             of ier.gt.0.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             or (integr.ne.1 and integr.ne.2) or
//                             icall.lt.1 or maxp1.lt.1.
//                             result, abserr, neval, last, rlist(1),
//                             elist(1), iord(1) and nnlog(1) are set
//                             to zero. alist(1) and blist(1) are set
//                             to a and b respectively.
//
//            last  -  integer
//                     on return, last equals the number of
//                     subintervals produces in the subdivision
//                     process, which determines the number of
//                     significant elements actually in the
//                     work arrays.
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
//                     vector of dimension at least limit, the first
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
//                     elements of which are pointers to the error
//                     estimates over the subintervals,
//                     such that elist(iord(1)), ...,
//                     elist(iord(k)) form a decreasing sequence, with
//                     k = last if last.le.(limit/2+2), and
//                     k = limit+1-last otherwise.
//
//            nnlog  - integer
//                     vector of dimension at least limit, containing the
//                     subdivision levels of the subintervals, i.e.
//                     iwork(i) = l means that the subinterval
//                     numbered i is of length abs(b-a)*2**(1-l)
//
//         on entry and return
//            momcom - integer
//                     indicating that the chebyshev moments
//                     have been computed for intervals of lengths
//                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
//                     momcom.lt.maxp1
//
//            chebmo - double precision
//                     array of dimension (maxp1,25) containing the
//                     chebyshev moments
//
//***references  (none)
//***routines called  d1mach,dqc25f,dqelg,dqpsrt
//***end prologue  dqawoe
//
    Real abseps, area, area1, area12, area2, a1,
        a2, b1, b2, correc, defab1, defab2, defabs,
        domega, dres, epmach, erlarg, erlast,
        errbnd, errmax, error1, erro12, error2, errsum, ertest, oflow,
        resabs, reseps, small, uflow, width;
    int id=1, ierro=1, iroff1=1, iroff2=1, iroff3=1,
        jupbnd=1, k=1, ksgn=1, ktmin=1, maxerr=1, nev=1,
        nres=1, nrmax=1, nrmom=1, numrl2=1;
    bool extrap, noext, extall;
//
    auto alist = dimension(alist_, limit);
    auto blist = dimension(blist_, limit);
    auto rlist = dimension(rlist_, limit);
    auto elist = dimension(elist_, limit);
    auto iord = dimension(iord_, limit);
    Real[52] rlist2_;
    auto rlist2 = dimension(rlist2_.ptr, 52);
    Real[3] res3la_;
    auto res3la = dimension(res3la_.ptr, 3);
    auto chebmo = dimension(chebmo_, maxp1, 25);
    auto nnlog = dimension(nnlog_, limit);
//
//            the dimension of rlist2 is determined by  the value of
//            limexp in subroutine dqelg (rlist2 should be of
//            dimension (limexp+2) at least).
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
//           rlist2    - array of dimension at least limexp+2
//                       containing the part of the epsilon table
//                       which is still needed for further computations
//           elist(i)  - error estimate applying to rlist(i)
//           maxerr    - pointer to the interval with largest
//                       error estimate
//           errmax    - elist(maxerr)
//           erlast    - error on the interval currently subdivided
//           area      - sum of the integrals over the subintervals
//           errsum    - sum of the errors over the subintervals
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//           *****1    - variable for the left subinterval
//           *****2    - variable for the right subinterval
//           last      - index for subdivision
//           nres      - number of calls to the extrapolation routine
//           numrl2    - number of elements in rlist2. if an appropriate
//                       approximation to the compounded integral has
//                       been obtained it is put in rlist2(numrl2) after
//                       numrl2 has been increased by one
//           small     - length of the smallest interval considered
//                       up to now, multiplied by 1.5
//           erlarg    - sum of the errors over the intervals larger
//                       than the smallest interval considered up to now
//           extrap    - logical variable denoting that the routine is
//                       attempting to perform extrapolation, i.e. before
//                       subdividing the smallest interval we try to
//                       decrease the value of erlarg
//           noext     - logical variable denoting that extrapolation
//                       is no longer allowed (true  value)
//
//            machine dependent constants
//            ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//           oflow is the largest positive magnitude.
//
//***first executable statement  dqawoe
      epmach = Real.epsilon;
//
//         test on validity of parameters
//         ------------------------------
//
      ier = 0;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      alist[1] = a;
      blist[1] = b;
      rlist[1] = 0.0;
      elist[1] = 0.0;
      iord[1] = 0;
      nnlog[1] = 0;
      if ((integr != 1 && integr != 2)  ||  (epsabs <= 0.0 &&
        epsrel < 0.5e2*epmach)  ||  icall < 1  ||
        maxp1 < 1)  ier = 6;
      if (ier == 6) goto l999;
//
//           first approximation to the integral
//           -----------------------------------
//
      domega = fabs(omega);
      nrmom = 0;
      if (icall > 1) goto l5;
      momcom = 0;
  l5: qc25f!(Real, Func)(f, a, b, domega, integr, nrmom, maxp1, 0, result, abserr,
        neval, defabs, resabs, momcom, chebmo.ptr);
//
//           test on accuracy.
//
      dres = fabs(result);
      errbnd = max(epsabs, epsrel*dres);
      rlist[1] = result;
      elist[1] = abserr;
      iord[1] = 1;
      if (abserr <= 0.1e3*epmach*defabs  &&  abserr > errbnd) ier = 2;
      if (limit == 1) ier = 1;
      if (ier != 0  ||  abserr <= errbnd) goto l200;
//
//           initializations
//           ---------------
//
      uflow = Real.min_normal;
      oflow = Real.max;
      errmax = abserr;
      maxerr = 1;
      area = result;
      errsum = abserr;
      abserr = oflow;
      nrmax = 1;
      extrap = false;
      noext = false;
      ierro = 0;
      iroff1 = 0;
      iroff2 = 0;
      iroff3 = 0;
      ktmin = 0;
      small = fabs(b-a)*0.75;
      nres = 0;
      numrl2 = 0;
      extall = false;
      if (0.5*fabs(b-a)*domega > 0.2e1) goto l10;
      numrl2 = 1;
      extall = true;
      rlist2[1] = result;
 l10: if (0.25*fabs(b-a)*domega <= 0.2e1) extall = true;
      ksgn = -1;
      if(dres >= (0.1e1-0.5e2*epmach)*defabs) ksgn = 1;
//
//           main do-loop
//           ------------
//
      for (last = 2; last <= limit; last++) // end: l140
      {
//
//           bisect the subinterval with the nrmax-th largest
//           error estimate.
//
        nrmom = nnlog[maxerr]+1;
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        qc25f!(Real, Func)(f,a1,b1,domega,integr,nrmom,maxp1,0,
          area1,error1,nev,resabs,defab1,momcom,chebmo.ptr);
        neval = neval+nev;
        qc25f!(Real, Func)(f,a2,b2,domega,integr,nrmom,maxp1,1,
          area2,error2,nev,resabs,defab2,momcom,chebmo.ptr);
        neval = neval+nev;
//
//           improve previous approximations to integral
//           and error and test for accuracy.
//
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if (defab1 == error1 || defab2 == error2) goto l25;
        if (fabs(rlist[maxerr]-area12) > 0.1e-4*fabs(area12)
          || erro12 < 0.99*errmax) goto l20;
        if (extrap) iroff2 = iroff2+1;
        else iroff1 = iroff1+1;
 l20:   if (last > 10 && erro12 > errmax) iroff3 = iroff3+1;
 l25:   rlist[maxerr] = area1;
        rlist[last] = area2;
        nnlog[maxerr] = nrmom;
        nnlog[last] = nrmom;
        errbnd = max(epsabs, epsrel*fabs(area));
//
//           test for roundoff error and eventually set error flag.
//
        if (iroff1+iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
//
//           set error flag in the case that the number of
//           subintervals equals limit.
//
        if (last == limit) ier = 1;
//
//           set error flag in the case of bad integrand behaviour
//           at a point of the integration range.
//
        if (max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*epmach)
          * (fabs(a2)+0.1e4*uflow)) ier = 4;
//
//           append the newly-created intervals to the list.
//
        if(error2 > error1) goto l30;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l40;
 l30:   alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
//
//           call subroutine dqpsrt to maintain the descending ordering
//           in the list of error estimates and select the subinterval
//           with nrmax-th largest error estimate (to bisected next).
//
 l40:   qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
// ***jump out of do-loop
      if (errsum <= errbnd) goto l170;
      if (ier != 0) goto l150;
        if (last == 2 && extall) goto l120;
        if (noext) goto l140;
        if (!extall) goto l50;
        erlarg = erlarg-erlast;
        if (fabs(b1-a1) > small) erlarg = erlarg+erro12;
        if (extrap) goto l70;
//
//           test whether the interval to be bisected next is the
//           smallest interval.
//
 l50:   width = fabs(blist[maxerr]-alist[maxerr]);
        if (width > small) goto l140;
        if (extall) goto l60;
//
//           test whether we can start with the extrapolation procedure
//           (we do this if we integrate over the next interval with
//           use of a gauss-kronrod rule - see subroutine dqc25f).
//
        small = small*0.5;
        if (0.25*width*domega > 0.2e1) goto l140;
        extall = true;
        goto l130;
 l60:   extrap = true;
        nrmax = 2;
 l70:   if (ierro == 3 || erlarg <= ertest) goto l90;
//
//           the smallest interval has the largest error.
//           before bisecting decrease the sum of the errors over
//           the larger intervals (erlarg) and perform extrapolation.
//
        jupbnd = last;
        if (last > (limit/2+2)) jupbnd = limit+3-last;
        id = nrmax;
        for (k = id; k<=jupbnd; k++)     // end: l80
        {
          maxerr = iord[nrmax];
          errmax = elist[maxerr];
          if(fabs(blist[maxerr]-alist[maxerr]) > small) goto l140;
          nrmax = nrmax+1;
 l80:;  }
//
//           perform extrapolation.
//
 l90:   numrl2 = numrl2+1;
        rlist2[numrl2] = area;
        if (numrl2 < 3) goto l110;
        qelg!Real(numrl2,rlist2.ptr,reseps,abseps,res3la.ptr,nres);
        ktmin = ktmin+1;
        if (ktmin > 5 && abserr < 0.1e-2*errsum) ier = 5;
        if (abseps >= abserr) goto l100;
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs, epsrel*fabs(reseps));
// ***jump out of do-loop
        if (abserr <= ertest) goto l150;
//
//           prepare bisection of the smallest interval.
//
l100:   if (numrl2 == 1) noext = true;
        if (ier == 5) goto l150;
l110:   maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = false;
        small = small*0.5;
        erlarg = errsum;
        goto l140;
l120:   small = small*0.5;
        numrl2 = numrl2+1;
        rlist2[numrl2] = area;
l130:   ertest = errbnd;
        erlarg = errsum;
l140:;}
//
//           set the final result.
//           ---------------------
//
l150: if (abserr == oflow || nres == 0) goto l170;
      if (ier+ierro == 0) goto l165;
      if (ierro == 3) abserr = abserr+correc;
      if (ier == 0) ier = 3;
      if (result != 0.0 && area != 0.0) goto l160;
      if (abserr > errsum) goto l170;
      if (area == 0.0) goto l190;
      goto l165;
l160: if (abserr/fabs(result) > errsum/fabs(area)) goto l170;
//
//           test on divergence.
//
l165: if (ksgn == (-1) && max(fabs(result),fabs(area)) <= 
        defabs*0.1e-1) goto l190;
      if (0.1e-1 > (result/area) || (result/area) > 0.1e3
        || errsum >= fabs(area)) ier = 6;
      goto l190;
//
//           compute global integral sum.
//
l170: result = 0.0;
      for (k=1; k<=last; k++)  // end: l180
      {
        result = result+rlist[k];
l180:;}
      abserr = errsum;
l190: if (ier > 2) ier=ier-1;
l200: if (integr == 2 && omega < 0.0) result=-result;
l999: return;
}


unittest
{
    alias qawoe!(float, float delegate(float)) fqawoe;
    alias qawoe!(double, double delegate(double)) dqawoe;
    alias qawoe!(double, double function(double)) dfqawoe;
    alias qawoe!(real, real delegate(real)) rqawoe;
}
