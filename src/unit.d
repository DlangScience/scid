#!/usr/local/bin/rdmd --shebang --force -w -d -unittest -L-lgfortran -L-lblas -L-llapack

// Run unittests.
module unit;


import std.stdio;


import scid.bindings.blas.blas;
import scid.bindings.blas.dblas;
import scid.bindings.blas.types;

import scid.ports.intde.intde1;
import scid.ports.intde.intde2;

import scid.ports.minpack.dogleg;
import scid.ports.minpack.enorm;
import scid.ports.minpack.fdjac1;
import scid.ports.minpack.hybrd;
import scid.ports.minpack.qform;
import scid.ports.minpack.qrfac;
import scid.ports.minpack.r1mpyq;
import scid.ports.minpack.r1updt;

import scid.ports.napack.addchg;
import scid.ports.napack.quasi;
import scid.ports.napack.stopit;

import scid.ports.quadpack.qag;
import scid.ports.quadpack.qage;
import scid.ports.quadpack.qagi;
import scid.ports.quadpack.qagie;
import scid.ports.quadpack.qagp;
import scid.ports.quadpack.qagpe;
import scid.ports.quadpack.qags;
import scid.ports.quadpack.qagse;
import scid.ports.quadpack.qawc;
import scid.ports.quadpack.qawce;
import scid.ports.quadpack.qawf;
import scid.ports.quadpack.qawfe;
import scid.ports.quadpack.qawo;
import scid.ports.quadpack.qawoe;
import scid.ports.quadpack.qaws;
import scid.ports.quadpack.qawse;
import scid.ports.quadpack.qc25c;
import scid.ports.quadpack.qc25f;
import scid.ports.quadpack.qc25s;
import scid.ports.quadpack.qcheb;
import scid.ports.quadpack.qelg;
import scid.ports.quadpack.qng;
import scid.ports.quadpack.qk15;
import scid.ports.quadpack.qk15i;
import scid.ports.quadpack.qk15w;
import scid.ports.quadpack.qk21;
import scid.ports.quadpack.qk31;
import scid.ports.quadpack.qk41;
import scid.ports.quadpack.qk51;
import scid.ports.quadpack.qk61;
import scid.ports.quadpack.qmomo;
import scid.ports.quadpack.qpsrt;
import scid.ports.quadpack.qwgtc;
import scid.ports.quadpack.qwgtf;
import scid.ports.quadpack.qwgts;

import scid.core.fortran;
import scid.core.memory;
import scid.core.meta;
import scid.core.testing;
import scid.core.traits;

import scid.internal.calculus.integrate_qng;

import scid.calculus;
import scid.exception;
import scid.linalg;
import scid.matrix;
import scid.nonlinear;
import scid.types;
import scid.util;


void main()
{
    version(unittest)  {
        auto rep = checkReport();
        writefln("\n%s checks performed, %s succeeded and %s failed",
            rep.total, rep.success, rep.fail);
    }
    else  { writeln("Compile with -unittest to run unit tests."); }
}

