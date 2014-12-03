/** BLAS bindings for D.

    Authors:    William V. Baxter III (with slight modifications by
                Lars Tandle Kyllingstad).
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.bindings.blas.types;


/// aliases to easily cope with "strange" blas

alias float f_float;
alias double f_double;
alias cfloat f_cfloat;
alias cdouble f_cdouble;
alias int f_int;

template isBlasType(T){
   const bool isBlasType=is(T==f_float)|| is(T==f_double) ||
       is(T==f_cfloat) || is(T==f_cdouble);
}
