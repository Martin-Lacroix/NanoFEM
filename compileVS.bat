:: Sets the environment variables and 64 bits compiler

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\Common7\Tools\vsdevcmd.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"

:: Compile the code and delete objects

cl /I C:\ProgramData\Alglib\src /EHsc /O2 ^
/DAE_CPU=AE_INTEL /DAE_OS=AE_WINDOWS /DAE_MKL ^
C:\ProgramData\Alglib\src\*.cpp source\*.cpp alglib316_64mkl.lib
del *.obj