g++ -msse2 -I C:\ProgramData\Alglib\src -o Mechanics.exe -O3 ^
-DAE_OS=AE_WINDOWS -DAE_CPU=AE_INTEL -pthread ^
-static -static-libgcc -static-libstdc++ ^
C:\ProgramData\Alglib\src\*.cpp source\*.cpp -lm