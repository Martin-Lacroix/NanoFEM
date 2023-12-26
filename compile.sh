g++ -msse2 -I /opt/alglib/src -o Mechanics -O3 \
-pthread -static -static-libgcc -static-libstdc++ \
/opt/alglib/src/*.cpp source/*.cpp -lm