g++ -msse2 -I C:\Local\Alglib\src -o Mechanics.exe -O3 ^
-DAE_OS=AE_WINDOWS -DAE_CPU=AE_INTEL -pthread ^
-static -static-libgcc -static-libstdc++ ^
C:\Local\Alglib\src\*.cpp source\*.cpp -lm
