set PATH=C:\D\dmd2\\windows\bin;C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\\bin;%PATH%

echo bindings\blas\blas.d >Release\scid.build.rsp
echo bindings\blas\dblas.d >>Release\scid.build.rsp
echo bindings\blas\types.d >>Release\scid.build.rsp
echo bindings\lapack\dlapack.d >>Release\scid.build.rsp
echo bindings\lapack\lapack.d >>Release\scid.build.rsp
echo common\fortran.d >>Release\scid.build.rsp
echo common\memory.d >>Release\scid.build.rsp
echo common\meta.d >>Release\scid.build.rsp
echo common\testing.d >>Release\scid.build.rsp
echo common\traits.d >>Release\scid.build.rsp
echo internal\calculus\integrate_qk.d >>Release\scid.build.rsp
echo internal\calculus\integrate_qng.d >>Release\scid.build.rsp
echo ports\intde\intde1.d >>Release\scid.build.rsp
echo ports\intde\intde2.d >>Release\scid.build.rsp
echo ports\linpack\gtsl.d >>Release\scid.build.rsp
echo ports\minpack\dogleg.d >>Release\scid.build.rsp
echo ports\minpack\enorm.d >>Release\scid.build.rsp
echo ports\minpack\fdjac1.d >>Release\scid.build.rsp
echo ports\minpack\hybrd.d >>Release\scid.build.rsp
echo ports\minpack\qform.d >>Release\scid.build.rsp
echo ports\minpack\qrfac.d >>Release\scid.build.rsp
echo ports\minpack\r1mpyq.d >>Release\scid.build.rsp
echo ports\minpack\r1updt.d >>Release\scid.build.rsp
echo ports\napack\addchg.d >>Release\scid.build.rsp
echo ports\napack\quasi.d >>Release\scid.build.rsp
echo ports\napack\stopit.d >>Release\scid.build.rsp
echo ports\quadpack\qag.d >>Release\scid.build.rsp
echo ports\quadpack\qage.d >>Release\scid.build.rsp
echo ports\quadpack\qagi.d >>Release\scid.build.rsp
echo ports\quadpack\qagie.d >>Release\scid.build.rsp
echo ports\quadpack\qagp.d >>Release\scid.build.rsp
echo ports\quadpack\qagpe.d >>Release\scid.build.rsp
echo ports\quadpack\qags.d >>Release\scid.build.rsp
echo ports\quadpack\qagse.d >>Release\scid.build.rsp
echo ports\quadpack\qawc.d >>Release\scid.build.rsp
echo ports\quadpack\qawce.d >>Release\scid.build.rsp
echo ports\quadpack\qawf.d >>Release\scid.build.rsp
echo ports\quadpack\qawfe.d >>Release\scid.build.rsp
echo ports\quadpack\qawo.d >>Release\scid.build.rsp
echo ports\quadpack\qawoe.d >>Release\scid.build.rsp
echo ports\quadpack\qaws.d >>Release\scid.build.rsp
echo ports\quadpack\qawse.d >>Release\scid.build.rsp
echo ports\quadpack\qc25c.d >>Release\scid.build.rsp
echo ports\quadpack\qc25f.d >>Release\scid.build.rsp
echo ports\quadpack\qc25s.d >>Release\scid.build.rsp
echo ports\quadpack\qcheb.d >>Release\scid.build.rsp
echo ports\quadpack\qelg.d >>Release\scid.build.rsp
echo ports\quadpack\qk15.d >>Release\scid.build.rsp
echo ports\quadpack\qk15i.d >>Release\scid.build.rsp
echo ports\quadpack\qk15w.d >>Release\scid.build.rsp
echo ports\quadpack\qk21.d >>Release\scid.build.rsp
echo ports\quadpack\qk31.d >>Release\scid.build.rsp
echo ports\quadpack\qk41.d >>Release\scid.build.rsp
echo ports\quadpack\qk51.d >>Release\scid.build.rsp
echo ports\quadpack\qk61.d >>Release\scid.build.rsp
echo ports\quadpack\qmomo.d >>Release\scid.build.rsp
echo ports\quadpack\qng.d >>Release\scid.build.rsp
echo ports\quadpack\qpsrt.d >>Release\scid.build.rsp
echo ports\quadpack\qwgtc.d >>Release\scid.build.rsp
echo ports\quadpack\qwgtf.d >>Release\scid.build.rsp
echo ports\quadpack\qwgts.d >>Release\scid.build.rsp
echo main.d >>Release\scid.build.rsp
echo calculus.d >>Release\scid.build.rsp
echo constants.d >>Release\scid.build.rsp
echo exception.d >>Release\scid.build.rsp
echo functions.d >>Release\scid.build.rsp
echo linalg.d >>Release\scid.build.rsp
echo matrix.d >>Release\scid.build.rsp
echo nonlinear.d >>Release\scid.build.rsp
echo types.d >>Release\scid.build.rsp
echo util.d >>Release\scid.build.rsp
echo matclosure.d >>Release\scid.build.rsp
echo expressions.d >>Release\scid.build.rsp
echo vecclosure.d >>Release\scid.build.rsp
echo scalclosure.d >>Release\scid.build.rsp

dmd -release -X -Xf"Release\scid.json" -of"Release\scid.exe" -deps="Release\scid.dep" -map "Release\scid.map" -L/NOMAP @Release\scid.build.rsp
if errorlevel 1 goto reportError
if not exist "Release\scid.exe" (echo "Release\scid.exe" not created! && goto reportError)

goto noError

:reportError
echo Building Release\scid.exe failed!

:noError
