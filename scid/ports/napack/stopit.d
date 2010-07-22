/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.napack.stopit;


import std.math;




//
//      ________________________________________________________
//     |                                                        |
//     |                 TEST FOR CONVERGENCE                   |
//     |                                                        |
//     |    INPUT:                                              |
//     |                                                        |
//     |         DIF   --ABSOLUTE DIFFERENCE                    |
//     |                                                        |
//     |         SIZE  --ABSOLUTE VALUE                         |
//     |                                                        |
//     |         NDIGIT--DESIRED NUMBER CORRECT DIGITS          |
//     |                                                        |
//     |         LIMIT --MAXIMUM NUMBER ITERATIONS              |
//     |                                                        |
//     |    OUTPUT:                                             |
//     |                                                        |
//     |         DIF   --POSITIVE TO CONTINUE ITERATIONS        |
//     |                                                        |
//     |         SIZE  --ABSOLUTE DIFFERENCE                    |
//     |                                                        |
//     |    BUILTIN FUNCTIONS: ABS                              |
//     |________________________________________________________|
///
void stopit(Real)(ref Real dif, ref Real size, int ndigit, int limit)
{
    
      Real e;
      static Real t;
      static int i = 0;
      dif = fabs(dif);
      size = fabs(size);
//     -----------------------------------------------
//     |*** INITIALIZATION DURING FIRST ITERATION ***|
//     -----------------------------------------------
      if (i > 0) goto l10;
      t = 10.0^^(-ndigit);
//     ------------------------------
//     |*** STOPPING CRITERION I ***|
//     ------------------------------
l10:  i++;
      e = 3*i;
      if (dif > t*size) goto l20;
      e += 1.0;
      goto l30;
//     -------------------------------
//     |*** STOPPING CRITERION II ***|
//     -------------------------------
l20:  if (i < limit) goto l40;
      e += 2.0;
l30:  dif = -dif;
      i = 0;
l40:  size = e;
}


unittest
{
    alias stopit!float fstopit;
    alias stopit!double dstopit;
    alias stopit!real rstopit;
}
