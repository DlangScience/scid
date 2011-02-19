/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qpsrt;


import scid.core.fortran;




///
void qpsrt(Real)(int limit, int last, ref int maxerr, ref Real ermax,
    Real* elist_, int* iord_, ref int nrmax)
{
//***begin prologue  dqpsrt
//***refer to  dqage,dqagie,dqagpe,dqawse
//***routines called  (none)
//***revision date  810101   (yymmdd)
//***keywords  sequential sorting
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  this routine maintains the descending ordering in the
//            list of the local error estimated resulting from the
//            interval subdivision process. at each call two error
//            estimates are inserted using the sequential search
//            method, top-down for the largest error estimate and
//            bottom-up for the smallest error estimate.
//***description
//
//           ordering routine
//           standard fortran subroutine
//           double precision version
//
//           parameters (meaning at output)
//              limit  - integer
//                       maximum number of error estimates the list
//                       can contain
//
//              last   - integer
//                       number of error estimates currently in the list
//
//              maxerr - integer
//                       maxerr points to the nrmax-th largest error
//                       estimate currently in the list
//
//              ermax  - double precision
//                       nrmax-th largest error estimate
//                       ermax = elist(maxerr)
//
//              elist  - double precision
//                       vector of dimension last containing
//                       the error estimates
//
//              iord   - integer
//                       vector of dimension last, the first k elements
//                       of which contain pointers to the error
//                       estimates, such that
//                       elist(iord(1)),...,  elist(iord(k))
//                       form a decreasing sequence, with
//                       k = last if last.le.(limit/2+2), and
//                       k = limit+1-last otherwise
//
//              nrmax  - integer
//                       maxerr = iord(nrmax)
//
//***end prologue  dqpsrt
//
      Real errmax,errmin;
      int i=1,ibeg=1,ido=1,isucc=1,j=1,jbnd=1,jupbn=1,k=1;
      auto elist = dimension(elist_, last);
      auto iord = dimension(iord_, last);
//
//           check whether the list contains more than
//           two error estimates.
//
//***first executable statement  dqpsrt
      if (last > 2) goto l10;
      iord[1] = 1;
      iord[2] = 2;
      goto l90;
//
//           this part of the routine is only executed if, due to a
//           difficult integrand, subdivision increased the error
//           estimate. in the normal case the insert procedure should
//           start after the nrmax-th largest error estimate.
//
 l10: errmax = elist[maxerr];
      if (nrmax == 1) goto l30;
      ido = nrmax-1;
      for (i=1; i<=ido; i++) { // end: l20
        isucc = iord[nrmax-1];
// ***jump out of do-loop
        if (errmax <= elist[isucc]) goto l30;
        iord[nrmax] = isucc;
        nrmax = nrmax-1;
 l20:;}
//
//           compute the number of elements in the list to be maintained
//           in descending order. this number depends on the number of
//           subdivisions still allowed.
//
 l30: jupbn = last;
      if (last > (limit/2+2)) jupbn = limit+3-last;
      errmin = elist[last];
//
//           insert errmax by traversing the list top-down,
//           starting comparison from the element elist(iord(nrmax+1)).
//
      jbnd = jupbn-1;
      ibeg = nrmax+1;
      if (ibeg > jbnd) goto l50;
      for (i=ibeg; i<=jbnd; i++) { // end: l40
        isucc = iord[i];
// ***jump out of do-loop
        if (errmax >= elist[isucc]) goto l60;
        iord[i-1] = isucc;
 l40:;}  
 l50: iord[jbnd] = maxerr;
      iord[jupbn] = last;
      goto l90;
//
//           insert errmin by traversing the list bottom-up.
//
 l60: iord[i-1] = maxerr;
      k = jbnd;
      for (j=i; j<=jbnd; j++) { // end: l70
        isucc = iord[k];
// ***jump out of do-loop
        if (errmin < elist[isucc]) goto l80;
        iord[k+1] = isucc;
        k = k-1;
 l70:;}
      iord[i] = last;
      goto l90;
 l80: iord[k+1] = last;
//
//           set maxerr and ermax.
//
 l90: maxerr = iord[nrmax];
      ermax = elist[maxerr];
      return;
}

unittest
{
    alias qpsrt!float fqpsrt;
    alias qpsrt!double dqpsrt;
    alias qpsrt!real rqpsrt;
}
