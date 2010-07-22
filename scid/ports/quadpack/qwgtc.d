/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qwgtc;


///
Real qwgtc(Real)(Real x, Real c, Real p2, Real p3, Real p4, int kp)
{
//***begin prologue  dqwgtc
//***refer to dqk15w
//***routines called  (none)
//***revision date  810101   (yymmdd)
//***keywords  weight function, cauchy principal value
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  this function subprogram is used together with the
//            routine qawc and defines the weight function.
//***end prologue  dqwgtc
//
//***first executable statement  dqwgtc
      return 1.0/(x-c);
}


unittest
{
    alias qwgtc!float fqwgtc;
    alias qwgtc!double dqwgtc;
    alias qwgtc!real rqwgtc;
}
