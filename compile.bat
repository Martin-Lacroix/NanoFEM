g++ -I C:\ProgramData\Alglib\src -L C:\Users\ORBBE\Desktop\Plugin\include\ ^
source\*.cpp C:\ProgramData\Alglib\src\*.cpp -o Mechanics.exe -O3 ^
-DAE_CPU=AE_INTEL -DAE_OS=AE_WINDOWS -pthread -DAE_MKL alglib*64mkl.lib

::g++ -O3 -I C:\ProgramData\Alglib\src -L C:\Users\ORBBE\Desktop\Plugin\include\ ^
::source\*.cpp C:\ProgramData\Alglib\src\*.cpp -o Mechanics.exe