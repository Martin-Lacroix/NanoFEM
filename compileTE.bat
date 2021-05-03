:: Sets the environment variables and 64 bits compiler

call C:\Windows\System32\cmd.exe /E:ON /K ""C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019"

:: Compile the code and delete objects

icx /I C:\ProgramData\Alglib\src /EHsc /O3 ^
/DAE_CPU=AE_INTEL /DAE_OS=AE_WINDOWS /DAE_MKL ^
C:\ProgramData\Alglib\src\*.cpp source\*.cpp alglib316_64mkl.lib
del *.obj