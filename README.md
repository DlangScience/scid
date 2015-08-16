[![Build Status](https://travis-ci.org/DlangScience/scid.svg?branch=master)](https://travis-ci.org/DlangScience/scid)
SciD
====

SciD is a collection of numerical routines and bindings written in and for
the D programming language.  It contains:

  * `scid.calculus`: Numerical integration (quadrature) and differentiation.
  * `scid.constants`: Fundamental constants of mathematics and Nature.
  * `scid.functions`: Mathematical special functions.
  * `scid.linalg`: Linear algebra (i.e., nice interfaces to a few LAPACK functions).
  * `scid.matrix`: LAPACK-compatible matrix view of ordinary arrays.
  * `scid.nonlinear`: Methods for nonlinear equation solving.

There are several ways to get SciD:

  * It is available as a [Dub package](http://code.dlang.org/packages/scid).
  * Users of Debian, Ubuntu, Linux Mint and similar operating systems may
    install it from the [D-APT package repository](http://d-apt.sourceforge.net/).
  * You can clone the [GitHub repository](https://github.com/DlangScience/scid)
    and build it yourself.

Please submit questions and bug reports using the
[GitHub issue tracker](https://github.com/DlangScience/scid/issues).

API documentation
-----------------

[Here](https://dlangscience.github.io/scid/api/)

Requirements
------------
All you need in order to use SciD is:

  * A somewhat up-to-date D compiler
  * The LAPACK library (and, by extension, a BLAS library)

Contributing
------------
Contributions are very welcome!  Please submit your work as
[pull requests on GitHub](https://github.com/DlangScience/scid/pulls).
By submitting a pull request, you implicitly accept that your contribution
is placed under the same licence as SciD.  See LICENCE.txt for details.

While there is no formal coding style guide for SciD, new code is expected to
follow the same style as existing code (which is more or less conventional
in the D community anyway).

Use four spaces for indentation (we don't need no stinkin' tabs!) , and UNIX
line endings (a single newline).
