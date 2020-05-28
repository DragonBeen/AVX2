CC = gcc
CFLAGS =-O3 -mavx2
CFLAGS1=-O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99
CFLAGS2=-O0 -mavx2

all :
	$(CC) $(CFLAGS) -c LizarMong_KEM_Comfort.c main.c randombytes.c fips202.c fips202x4.c xe5.c
	$(CC) $(CFLAGS) -o LizarMong_KEM_Comfort LizarMong_KEM_Comfort.o main.o randombytes.o fips202.o fips202x4.o xe5.o keccak4x/KeccakP-1600-times4-SIMD256.o
	
run : all
	./LizarMong_KEM_Comfort

clean :
	rm -f *.o
	rm -f LizarMong_KEM_Comfort

new :
	make clean
	make all
	./LizarMong_KEM_Comfort
