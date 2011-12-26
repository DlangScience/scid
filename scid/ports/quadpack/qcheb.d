// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qcheb;


import scid.core.fortran;




///
void qcheb(Real)(const Real* x_, Real* fval_, Real* cheb12_, Real* cheb24_)
{
//***begin prologue  dqcheb
//***refer to  dqc25c,dqc25f,dqc25s
//***routines called  (none)
//***revision date  830518   (yymmdd)
//***keywords  chebyshev series expansion, fast fourier transform
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  this routine computes the chebyshev series expansion
//            of degrees 12 and 24 of a function using a
//            fast fourier transform method
//            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
//            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
//            where t(k,x) is the chebyshev polynomial of degree k.
//***description
//
//        chebyshev series expansion
//        standard fortran subroutine
//        double precision version
//
//        parameters
//          on entry
//           x      - double precision
//                    vector of dimension 11 containing the
//                    values cos(k*pi/24), k = 1, ..., 11
//
//           fval   - double precision
//                    vector of dimension 25 containing the
//                    function values at the points
//                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
//                    where (a,b) is the approximation interval.
//                    fval(1) and fval(25) are divided by two
//                    (these values are destroyed at output).
//
//          on return
//           cheb12 - double precision
//                    vector of dimension 13 containing the
//                    chebyshev coefficients for degree 12
//
//           cheb24 - double precision
//                    vector of dimension 25 containing the
//                    chebyshev coefficients for degree 24
//
//***end prologue  dqcheb
//
      Real alam,alam1,alam2,part1,part2, part3;
      Real[12] v_;
      int i=1,j=1;
//
      auto cheb12 = dimension(cheb12_, 13);
      auto cheb24 = dimension(cheb24_, 25);
      auto fval = dimension(fval_, 25);
      auto v = dimension(v_.ptr, 12);
      auto x = dimension(x_, 11);
//
//***first executable statement  dqcheb
      for (i=1; i<=12; i++) { // end: 10
        j = 26-i;
        v[i] = fval[i]-fval[j];
        fval[i] = fval[i]+fval[j];
 l10: ;}
      alam1 = v[1]-v[9];
      alam2 = x[6]*(v[3]-v[7]-v[11]);
      cheb12[4] = alam1+alam2;
      cheb12[10] = alam1-alam2;
      alam1 = v[2]-v[8]-v[10];
      alam2 = v[4]-v[6]-v[12];
      alam = x[3]*alam1+x[9]*alam2;
      cheb24[4] = cheb12[4]+alam;
      cheb24[22] = cheb12[4]-alam;
      alam = x[9]*alam1-x[3]*alam2;
      cheb24[10] = cheb12[10]+alam;
      cheb24[16] = cheb12[10]-alam;
      part1 = x[4]*v[5];
      part2 = x[8]*v[9];
      part3 = x[6]*v[7];
      alam1 = v[1]+part1+part2;
      alam2 = x[2]*v[3]+part3+x[10]*v[11];
      cheb12[2] = alam1+alam2;
      cheb12[12] = alam1-alam2;
      alam = x[1]*v[2]+x[3]*v[4]+x[5]*v[6]+x[7]*v[8]
        +x[9]*v[10]+x[11]*v[12];
      cheb24[2] = cheb12[2]+alam;
      cheb24[24] = cheb12[2]-alam;
      alam = x[11]*v[2]-x[9]*v[4]+x[7]*v[6]-x[5]*v[8]
        +x[3]*v[10]-x[1]*v[12];
      cheb24[12] = cheb12[12]+alam;
      cheb24[14] = cheb12[12]-alam;
      alam1 = v[1]-part1+part2;
      alam2 = x[10]*v[3]-part3+x[2]*v[11];
      cheb12[6] = alam1+alam2;
      cheb12[8] = alam1-alam2;
      alam = x[5]*v[2]-x[9]*v[4]-x[1]*v[6]
        -x[11]*v[8]+x[3]*v[10]+x[7]*v[12];
      cheb24[6] = cheb12[6]+alam;
      cheb24[20] = cheb12[6]-alam;
      alam = x[7]*v[2]-x[3]*v[4]-x[11]*v[6]+x[1]*v[8]
        -x[9]*v[10]-x[5]*v[12];
      cheb24[8] = cheb12[8]+alam;
      cheb24[18] = cheb12[8]-alam;
      for (i=1; i<=6; i++) { // end: 20
        j = 14-i;
        v[i] = fval[i]-fval[j];
        fval[i] = fval[i]+fval[j];
 l20: ;}
      alam1 = v[1]+x[8]*v[5];
      alam2 = x[4]*v[3];
      cheb12[3] = alam1+alam2;
      cheb12[11] = alam1-alam2;
      cheb12[7] = v[1]-v[5];
      alam = x[2]*v[2]+x[6]*v[4]+x[10]*v[6];
      cheb24[3] = cheb12[3]+alam;
      cheb24[23] = cheb12[3]-alam;
      alam = x[6]*(v[2]-v[4]-v[6]);
      cheb24[7] = cheb12[7]+alam;
      cheb24[19] = cheb12[7]-alam;
      alam = x[10]*v[2]-x[6]*v[4]+x[2]*v[6];
      cheb24[11] = cheb12[11]+alam;
      cheb24[15] = cheb12[11]-alam;
      for (i=1; i<=3; i++) { // end: 30
        j = 8-i;
        v[i] = fval[i]-fval[j];
        fval[i] = fval[i]+fval[j];
 l30: ;}
      cheb12[5] = v[1]+x[8]*v[3];
      cheb12[9] = fval[1]-x[8]*fval[3];
      alam = x[4]*v[2];
      cheb24[5] = cheb12[5]+alam;
      cheb24[21] = cheb12[5]-alam;
      alam = x[8]*fval[2]-fval[4];
      cheb24[9] = cheb12[9]+alam;
      cheb24[17] = cheb12[9]-alam;
      cheb12[1] = fval[1]+fval[3];
      alam = fval[2]+fval[4];
      cheb24[1] = cheb12[1]+alam;
      cheb24[25] = cheb12[1]-alam;
      cheb12[13] = v[1]-v[3];
      cheb24[13] = cheb12[13];
      alam = 0.1e1/0.6e1;
      for(i=2; i<=12; i++) { // end: 40
        cheb12[i] = cheb12[i]*alam;
 l40: ;} 
      alam = 0.5*alam;
      cheb12[1] = cheb12[1]*alam;
      cheb12[13] = cheb12[13]*alam;
      for (i=2; i<=24; i++) { // end: 50
        cheb24[i] = cheb24[i]*alam;
 l50: ;}
      cheb24[1] = 0.5*alam*cheb24[1];
      cheb24[25] = 0.5*alam*cheb24[25];
      return;
}


unittest
{
    alias qcheb!float fqcheb;
    alias qcheb!double dqcheb;
    alias qcheb!real rqcheb;
}
