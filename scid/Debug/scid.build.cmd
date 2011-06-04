set PATH=C:\D\dmd2\\windows\bin;C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\\bin;%PATH%
set DMD_LIB=;C:\Users\Cristi\Documents\Code\D\scid\deps

echo bindings\blas\blas.d >Debug\scid.build.rsp
echo bindings\blas\dblas.d >>Debug\scid.build.rsp
echo bindings\blas\types.d >>Debug\scid.build.rsp
echo bindings\lapack\dlapack.d >>Debug\scid.build.rsp
echo bindings\lapack\lapack.d >>Debug\scid.build.rsp
echo common\fortran.d >>Debug\scid.build.rsp
echo common\memory.d >>Debug\scid.build.rsp
echo common\meta.d >>Debug\scid.build.rsp
echo common\testing.d >>Debug\scid.build.rsp
echo common\traits.d >>Debug\scid.build.rsp
echo internal\calculus\integrate_qk.d >>Debug\scid.build.rsp
echo internal\calculus\integrate_qng.d >>Debug\scid.build.rsp
echo ports\intde\intde1.d >>Debug\scid.build.rsp
echo ports\intde\intde2.d >>Debug\scid.build.rsp
echo ports\linpack\gtsl.d >>Debug\scid.build.rsp
echo ports\minpack\dogleg.d >>Debug\scid.build.rsp
echo ports\minpack\enorm.d >>Debug\scid.build.rsp
echo ports\minpack\fdjac1.d >>Debug\scid.build.rsp
echo ports\minpack\hybrd.d >>Debug\scid.build.rsp
echo ports\minpack\qform.d >>Debug\scid.build.rsp
echo ports\minpack\qrfac.d >>Debug\scid.build.rsp
echo ports\minpack\r1mpyq.d >>Debug\scid.build.rsp
echo ports\minpack\r1updt.d >>Debug\scid.build.rsp
echo ports\napack\addchg.d >>Debug\scid.build.rsp
echo ports\napack\quasi.d >>Debug\scid.build.rsp
echo ports\napack\stopit.d >>Debug\scid.build.rsp
echo ports\quadpack\qag.d >>Debug\scid.build.rsp
echo ports\quadpack\qage.d >>Debug\scid.build.rsp
echo ports\quadpack\qagi.d >>Debug\scid.build.rsp
echo ports\quadpack\qagie.d >>Debug\scid.build.rsp
echo ports\quadpack\qagp.d >>Debug\scid.build.rsp
echo ports\quadpack\qagpe.d >>Debug\scid.build.rsp
echo ports\quadpack\qags.d >>Debug\scid.build.rsp
echo ports\quadpack\qagse.d >>Debug\scid.build.rsp
echo ports\quadpack\qawc.d >>Debug\scid.build.rsp
echo ports\quadpack\qawce.d >>Debug\scid.build.rsp
echo ports\quadpack\qawf.d >>Debug\scid.build.rsp
echo ports\quadpack\qawfe.d >>Debug\scid.build.rsp
echo ports\quadpack\qawo.d >>Debug\scid.build.rsp
echo ports\quadpack\qawoe.d >>Debug\scid.build.rsp
echo ports\quadpack\qaws.d >>Debug\scid.build.rsp
echo ports\quadpack\qawse.d >>Debug\scid.build.rsp
echo ports\quadpack\qc25c.d >>Debug\scid.build.rsp
echo ports\quadpack\qc25f.d >>Debug\scid.build.rsp
echo ports\quadpack\qc25s.d >>Debug\scid.build.rsp
echo ports\quadpack\qcheb.d >>Debug\scid.build.rsp
echo ports\quadpack\qelg.d >>Debug\scid.build.rsp
echo ports\quadpack\qk15.d >>Debug\scid.build.rsp
echo ports\quadpack\qk15i.d >>Debug\scid.build.rsp
echo ports\quadpack\qk15w.d >>Debug\scid.build.rsp
echo ports\quadpack\qk21.d >>Debug\scid.build.rsp
echo ports\quadpack\qk31.d >>Debug\scid.build.rsp
echo ports\quadpack\qk41.d >>Debug\scid.build.rsp
echo ports\quadpack\qk51.d >>Debug\scid.build.rsp
echo ports\quadpack\qk61.d >>Debug\scid.build.rsp
echo ports\quadpack\qmomo.d >>Debug\scid.build.rsp
echo ports\quadpack\qng.d >>Debug\scid.build.rsp
echo ports\quadpack\qpsrt.d >>Debug\scid.build.rsp
echo ports\quadpack\qwgtc.d >>Debug\scid.build.rsp
echo ports\quadpack\qwgtf.d >>Debug\scid.build.rsp
echo ports\quadpack\qwgts.d >>Debug\scid.build.rsp
echo main.d >>Debug\scid.build.rsp
echo calculus.d >>Debug\scid.build.rsp
echo constants.d >>Debug\scid.build.rsp
echo exception.d >>Debug\scid.build.rsp
echo functions.d >>Debug\scid.build.rsp
echo linalg.d >>Debug\scid.build.rsp
echo matrix.d >>Debug\scid.build.rsp
echo nonlinear.d >>Debug\scid.build.rsp
echo types.d >>Debug\scid.build.rsp
echo util.d >>Debug\scid.build.rsp
echo matclosure.d >>Debug\scid.build.rsp
echo expressions.d >>Debug\scid.build.rsp
echo vecclosure.d >>Debug\scid.build.rsp
echo scalclosure.d >>Debug\scid.build.rsp

dmd -g -debug -X -Xf"Debug\scid.json" -of"Debug\scid.exe_cv" -deps="Debug\scid.dep" -map "Debug\scid.map" -L/NOMAP C:\Users\Cristi\Documents\Code\D\scid\deps\blaslapackdll.lib @Debug\scid.build.rsp
if errorlevel 1 goto reportError
if not exist "Debug\scid.exe_cv" (echo "Debug\scid.exe_cv" not created! && goto reportError)
echo Converting debug information...
"C:\Program Files (x86)\VisualD\cv2pdb\cv2pdb.exe" -D2 "Debug\scid.exe_cv" "Debug\scid.exe"
if errorlevel 1 goto reportError
if not exist "Debug\scid.exe" (echo "Debug\scid.exe" not created! && goto reportError)

goto noError

:reportError
echo Building Debug\scid.exe failed!

:noError
