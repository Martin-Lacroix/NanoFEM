:: g++ -I C:\ProgramData\Alglib\src -o Mechanics.exe -O3 ^
:: -DAE_OS=AE_WINDOWS -DAE_CPU=AE_INTEL -pthread -static -static-libgcc -static-libstdc++ ^
:: C:\ProgramData\Alglib\src\*.cpp source\*.cpp


g++ -I C:\ProgramData\Alglib\src -o Mechanics.exe -O3 ^
-DAE_OS=AE_WINDOWS -DAE_CPU=AE_INTEL -DAE_MKL ^
-pthread -static -static-libgcc -static-libstdc++ ^
C:\ProgramData\Alglib\src\*.cpp source\*.cpp alglib316_64mkl.lib