module scid.lapack;

version( nodeps ) {
	alias naive_ lapack;
} else {
	static import scid.bindings.lapack.dlapack;
	alias scid.bindings.lapack.dlapack lapack;
}
