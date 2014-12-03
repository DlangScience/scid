/** Various fundamental constants.  Physical constants are given in SI units.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
    Macros:
        SUB = <sub>$0</sub>
*/
module scid.constants;




// ==================== PHYSICAL CONSTANTS ====================

/// Boltzmann constant, $(I k$(SUB B))  [J/K]
enum real boltzmannConstant = 1.38064_88e-23L;

/// Planck constant, &#x210E; [J s]
enum real planckConstant = 6.62606_896e-34L;

/// Reduced Planck constant, &#x210F; &#x2261; &#x210E;/2&#x03C0;  [J s]
enum real hBar = 1.05457_1628e-34L;




// ==================== MATHEMATICAL CONSTANTS ====================

/// Pi squared, &#x03C0;&#x00B2;
enum real piSquared = 0x9.de9e64df22ef2d2p+0L;

/// Two times pi, 2&#x03C0;
enum real twoPi = 0xc.90fdaa22168c235p-1L;

/// Euler-Mascheroni constant, &#x03B3;
enum real eulerGamma = 0x9.3c467e37db0c7a5p-4L;

