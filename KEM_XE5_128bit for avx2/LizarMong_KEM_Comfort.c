#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

#include "LizarMong_KEM.h"
#include "randombytes.h"
#include "fips202.h"
#include "xef.h"
#include "fips202x4.h"

int Keygen(unsigned char *pk, unsigned char *sk){
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N*2]={0,};
	unsigned char seed_a[SEED_LEN];
	uint16_t sk_s[HS];
	int i, j;
////////// Gen poly a ////////// 
	randombytes(seed_a, SEED_LEN);	
	shake256x4(pk_a, pk_a+(LWE_N/4), pk_a+2*(LWE_N/4), pk_a+3*(LWE_N/4), LWE_N/4, seed_a, seed_a+8, seed_a+16, seed_a+24, SEED_LEN/4);

////////// Gen poly s and s_idx  //////////
	memset(sk, 0, LWE_N+LWE_N/8);
	unsigned char seed_s[HS*3], temp;
	unsigned int sk_random_idx;
	int hw=0, count = 0, neg_start=0, back_position = HS;

	randombytes(seed_s, HS*3);
	while (hw < HS) {
		sk_random_idx = seed_s[count++]; 
		sk_random_idx <<= 8;
		sk_random_idx ^= seed_s[count++];
		temp = (sk_random_idx & 0x02) - 1;
		sk_random_idx >>= 2;
		sk_random_idx &= (LWE_N - 1);
		if (sk[sk_random_idx] == 0) {
			sk[sk_random_idx] = temp;
			hw++;
			if (sk[sk_random_idx]==0x01){sk_s[neg_start++] = sk_random_idx;}
			if (sk[sk_random_idx]==0xff){sk_s[--back_position] = sk_random_idx;}
		}
		if (count >= HS*3 - 2) {
			randombytes(seed_s, HS*3);
			count = 0;
		}
	}
	if (hw != HS) { // fault detecting
		return 3;
	}

/* ADD compare PKE algorithm */

////////// Gen u and Concat sk = (sk || u) /////////
	unsigned char u[LWE_N/8];
	randombytes(u, LWE_N/8);
	memcpy(sk+LWE_N, u, LWE_N/8);

/*		END of the ADD		*/

////////// Initialize b as an error polynomial e ////////// 
	unsigned char b0, b1, tmp2[LWE_N/4];
	randombytes(tmp2,LWE_N/4);
	__m256i x256, y256, z256, w256;
	__m256i hex_01 = _mm256_set1_epi8(0x01);
	
	for(j=0; j<LWE_N/4; j+=16){
	// step 1
	x256 = _mm256_loadu_si256((__m256i_u*)&tmp2[j]);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_loadu_si256((__m256i_u*) &tmp2[j+8]);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+0], z256);
	// step 2
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+8], z256);
	// step 3
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+16], z256);
	// step 4
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+24], z256);
	// step 5
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+32], z256);
	// step 6
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+40], z256);
	// step 7
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+48], z256);
	// step 8
	x256 = _mm256_srli_epi16(x256, 1);
	y256 = _mm256_and_si256(x256, hex_01);
	w256 = _mm256_srli_epi16(w256, 1);
	z256 = _mm256_and_si256(w256, hex_01);

	z256 = _mm256_sub_epi8(y256, z256);
	_mm256_storeu_si256((__m256i*)&pk_b[(j*4)+56], z256);

	}

	
	if (j != (LWE_N/4)) { // fault detecting
		return 3;
	}

////////// mult a*s and add e(pk_b) ////////// 
	uint16_t deg;

 for (i = 0; i < neg_start; ++i){
        deg = sk_s[i];
        for ( j = 0; j < LWE_N; j+=32) {
        x256 = _mm256_loadu_si256((__m256i*) &pk_b[deg+j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
        z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&pk_b[deg+j], z256);
        }

    }
    for (i = neg_start; i < HS; ++i){
        deg = sk_s[i];
        for ( j = 0; j < LWE_N; j+=32) {
         x256 = _mm256_loadu_si256((__m256i*) &pk_b[deg+j]);
         y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
         z256 = _mm256_add_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&pk_b[deg+j], z256);
        }
    }



    //for (j = 0; j < LWE_N; ++j) {pk_b[j] -= pk_b[LWE_N + j];}

    for (j = 0; j < LWE_N; j+=32){
        x256 = _mm256_loadu_si256((__m256i*) &pk_b[j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_b[LWE_N+j]);
        z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&pk_b[j], z256);
    }


////////// Concat seed_genA || pk_b //////////
	memcpy(pk, seed_a, SEED_LEN);
	memcpy(pk+SEED_LEN, pk_b, LWE_N);

	return 0;
}



int Enc(unsigned char *c, unsigned char *shared_k, const unsigned char *pk){ 
	int i, j;
	unsigned char c1h_a[LWE_N*2]={0,};
	unsigned char c1h_b[LWE_N*2]={0,};	
	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL; 


////////// Generate a random polynomial delta ////////// 
	unsigned char delta[size_of_delta] = { 0, };
	int sum = 0;
	randombytes(delta, size_of_delta);
	for (i = 0; i < size_of_delta; ++i) {
		sum += delta[i];
	}
	if (sum == 0) { // fault detecting
		return 3;
	}

////////// Set r = H(delta) and Gen r_idx  ////////// 
	unsigned char r[LWE_N]={0,}, temp;
	uint16_t r_idx[HR];
	unsigned int r_random_idx; 
	int hw=0, count = 0, neg_start = 0, back_position = HR;

	hash = calloc(HR*3, sizeof(unsigned char));
	shake256(hash, HR*3, delta, size_of_delta);	
	//shake256x4(hash, hash+((HR*3)/4), hash+2*((HR*3)/4), hash+3*((HR*3)/4), (HR*3)/4, delta, delta+((size_of_delta)/4), delta+2*((size_of_delta)/4), delta+3*((size_of_delta)/4), size_of_delta/4);

	//shake256x4(pk_a, pk_a+(LWE_N/4), pk_a+2*(LWE_N/4), pk_a+3*(LWE_N/4), LWE_N/4, seed_a, seed_a+8, seed_a+16, seed_a+24, SEED_LEN/4);
	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		temp = (r_random_idx & 0x02) - 1;
		r_random_idx >>= 2;
		r_random_idx = r_random_idx & (LWE_N - 1);  
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = temp;
			hw++;
			if (r[r_random_idx] == 0x01){r_idx[neg_start++] = r_random_idx;}
			if (r[r_random_idx] == 0xff){r_idx[--back_position] = r_random_idx;}
		}
		if (count >= (HR*3 - 2)) { 
			shake256(hash, HR*3, hash, HR*3);
			count = 0;
		}
	}
	
////////// Encoding delta using Error Correcting Code ////////// 
	unsigned char delta_hat[LWE_N / 8]={0,};

	memcpy(delta_hat, delta, size_of_delta);
	xe5_234_compute(delta_hat);				 

////////// Parse seed_a||pk_b from pk and Make pk_a ////////// 
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N];
	unsigned char seed_a[SEED_LEN];

	memcpy(seed_a, pk, SEED_LEN);
	shake256x4(pk_a, pk_a+(LWE_N/4), pk_a+2*(LWE_N/4), pk_a+3*(LWE_N/4), LWE_N/4, seed_a, seed_a+8, seed_a+16, seed_a+24, SEED_LEN/4);

	memcpy(pk_b, pk+SEED_LEN, LWE_N);

////////// Initialize c1h_b as q/2 * delta_hat //////////  		
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;	
		}
	}


////////// Compute a * r and b * r, and then add to c1h_a and c1h_b, respectively. ////////// 
	uint16_t deg;

	__m256i x256, y256, z256;
 for (i = 0; i < neg_start; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
        x256 = _mm256_loadu_si256((__m256i*) &c1h_a[deg+j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
        z256 = _mm256_add_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_a[deg+j], z256);
	}
    }
    for (i = neg_start; i < HR; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
         x256 = _mm256_loadu_si256((__m256i*) &c1h_a[deg+j]);
         y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
         z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_a[deg+j], z256);
        }
    }
 for (i = 0; i < neg_start; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
        x256 = _mm256_loadu_si256((__m256i*) &c1h_b[deg+j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_b[j]);
        z256 = _mm256_add_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_b[deg+j], z256);
	}
    }
    for (i = neg_start; i < HR; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
         x256 = _mm256_loadu_si256((__m256i*) &c1h_b[deg+j]);
         y256 = _mm256_loadu_si256((__m256i*) &pk_b[j]);
         z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_b[deg+j], z256);
        }
    }

	__m256i two256 = _mm256_set1_epi8(0x02);
	__m256i eight256 = _mm256_set1_epi8(0x08);
	__m256i fc256 = _mm256_set1_epi8(0xfc);
	__m256i f0256 = _mm256_set1_epi8(0xf0);

    for (j = 0; j < LWE_N; j+=32){
        __m256i x256 = _mm256_loadu_si256((__m256i*) &c1h_a[j]);
        __m256i y256 = _mm256_loadu_si256((__m256i*) &c1h_a[LWE_N+j]);
        __m256i z256 = _mm256_sub_epi8(x256, y256);
	z256 = _mm256_add_epi8(z256, two256);
	z256 = _mm256_and_si256(z256, fc256);
        _mm256_storeu_si256((__m256i*)&c[j], z256);
        x256 = _mm256_loadu_si256((__m256i*) &c1h_b[j]);
        y256 = _mm256_loadu_si256((__m256i*) &c1h_b[LWE_N+j]);
        z256 = _mm256_sub_epi8(x256, y256);
	z256 = _mm256_add_epi8(z256, eight256);
	z256 = _mm256_and_si256(z256, f0256);
        _mm256_storeu_si256((__m256i*)&c[LWE_N+j], z256);
    }
	
////////// Send c1h_a and c1h_b from mod q to mod p and mod k ////////// 



/* ADD compare PKE algorithm */

////////// G(c1,delta) ////////// 
	hash_t=calloc((LWE_N+LWE_N+LWE_N/8), sizeof(unsigned char));
	memcpy(hash_t, c, (LWE_N+LWE_N));	
	memcpy(hash_t+(LWE_N+LWE_N), delta_hat, LWE_N/8);
	sha3_256(shared_k, hash_t, LWE_N+LWE_N+LWE_N/8);
	free(hash_t);

/*		END of the ADD		*/

	free(hash);
	return 0;
}




int Dec(unsigned char *shared_k, unsigned char *c, const unsigned char *sk, const unsigned char *pk){
	int res = 0;
	int i, j;

	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL;

	unsigned char delta_hat[LWE_N / 8]={0,};
	unsigned char delta[size_of_delta]={0,};
	unsigned char c1h_a[LWE_N*2] = { 0, };
	unsigned char c1h_b[LWE_N*2] = { 0, };
	unsigned char decomp_delta[LWE_N*2]={0,};

	uint16_t deg;
	__m256i x256, y256, z256;
////////// Initialize decomp_delta(=c2) as c1h_a(=c1) //////////
	memcpy(decomp_delta, c+LWE_N, LWE_N);
	memcpy(c1h_a, c, LWE_N);

//////// Omit the task of changing mod_k to mod_p. because the data_type is unsigned_char. //////// 


////////// Gen s_idx //////////
	uint16_t sk_s[HS];
	int neg_start = 0, back_position = HS;	
	for (i = 0; i < LWE_N; ++i) {
		if (sk[i] == 0x01){sk_s[neg_start++] = i;}
		if (sk[i] == 0xff){sk_s[--back_position] = i;}

			
	}
////////// Compute delta (c2 + c1 * s) ////////// 

	for (i = 0; i < neg_start; ++i){
		deg = sk_s[i];
		for ( j = 0; j < LWE_N; j+=32) {
			x256 = _mm256_loadu_si256((__m256i*) &decomp_delta[deg+j]);
			y256 = _mm256_loadu_si256((__m256i*) &c1h_a[j]);
			z256 = _mm256_add_epi8(x256, y256);
			_mm256_storeu_si256((__m256i*)&decomp_delta[deg+j], z256);
		}
	}
	for (i = neg_start; i < HS; ++i){
		deg = sk_s[i];
		for ( j = 0; j < LWE_N; j+=32) {
			x256 = _mm256_loadu_si256((__m256i*) &decomp_delta[deg+j]);
			y256 = _mm256_loadu_si256((__m256i*) &c1h_a[j]);
			z256 = _mm256_sub_epi8(x256, y256);
			_mm256_storeu_si256((__m256i*)&decomp_delta[deg+j], z256);
		}
	}
	__m256i hex_40 = _mm256_set1_epi8(0x40);
	for (j = 0; j < LWE_N; j+=32){
		x256 = _mm256_loadu_si256((__m256i*) &decomp_delta[j]);
		y256 = _mm256_loadu_si256((__m256i*) &decomp_delta[LWE_N+j]);
		z256 = _mm256_sub_epi8(x256, y256);
		z256 = _mm256_add_epi8(z256, hex_40);
		_mm256_storeu_si256((__m256i*)&decomp_delta[j], z256);
	}

//////// Compute delta_hat' = 2/p * delta //////// 
	for (i = 0; i < LWE_N; ++i) {
		//decomp_delta[i] += 0x40;
		 decomp_delta[i] >>= _8_LOG_T;
	}
////////// Set delta_hat' ////////// 
	for (i = 0; i < LWE_N/8; ++i) {
		for (j = 0; j < 8; ++j) {
			uint8_t a = (decomp_delta[8 * i + j]) << j;
			delta_hat[i] ^= a;  
		}
	}
	
////////// Decoding delta_hat using Error Correcting Code ////////// 
	xe5_234_compute(delta_hat);
	xe5_234_fixerr(delta_hat);
	memcpy(delta, delta_hat, size_of_delta);

////////// Set r = H(delta) and Gen r_idx ////////// 
	unsigned char r[LWE_N] = { 0, }, temp;
	unsigned int r_idx[HR];
	unsigned int r_random_idx; 
	int hw=0, count = 0;
	neg_start = 0;
	back_position = HR;
	hash = calloc(HR*3, sizeof(unsigned char));
	shake256(hash, HR*3, delta, size_of_delta);	
	//shake256x4(hash, hash+((HR*3)/4), hash+2*((HR*3)/4), hash+3*((HR*3)/4), (HR*3)/4, delta, delta+((size_of_delta)/4), delta+2*((size_of_delta)/4), delta+3*((size_of_delta)/4), size_of_delta/4);
	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		temp = (r_random_idx & 0x02) - 1;
		r_random_idx >>= 2;
		r_random_idx = r_random_idx & (LWE_N - 1); // (seed1) || (seed2) 
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = temp;
			hw++;
			if (r[r_random_idx] == 0x01){r_idx[neg_start++] = r_random_idx;}
			if (r[r_random_idx] == 0xff){r_idx[--back_position] = r_random_idx;}
		}
		if (count >= (HR*3 - 2)) { 
			shake256(hash, HR*3, hash, HR*3);
			count = 0;
		}
	}
	
////////// Encoding delta using Error Correcting Code //////////
	memset(delta_hat, 0, LWE_N/8);
	memcpy(delta_hat, delta, size_of_delta);
	xe5_234_compute(delta_hat);

////////// Parse seed_a||pk_b from pk & Make pk_a ////////// 
	unsigned char pk_a[LWE_N], pk_b[LWE_N], seed_a[SEED_LEN];

	memcpy(seed_a, pk, SEED_LEN);
	shake256x4(pk_a, pk_a+(LWE_N/4), pk_a+2*(LWE_N/4), pk_a+3*(LWE_N/4), LWE_N/4, seed_a, seed_a+8, seed_a+16, seed_a+24, SEED_LEN/4);

	memcpy(pk_b, pk+SEED_LEN, LWE_N);

////////// Initialize c1h_b as q/2 * delta_hat //////////  	
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;
		}
	}

////////// Compute a * r and b * r, and then add to c2h_a and c2h_b, respectively. ///////// 
	memset(c1h_a, 0, LWE_N*2);


 for (i = 0; i < neg_start; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
        x256 = _mm256_loadu_si256((__m256i*) &c1h_a[deg+j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
        z256 = _mm256_add_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_a[deg+j], z256);
	}
    }
    for (i = neg_start; i < HR; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
         x256 = _mm256_loadu_si256((__m256i*) &c1h_a[deg+j]);
         y256 = _mm256_loadu_si256((__m256i*) &pk_a[j]);
         z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_a[deg+j], z256);
        }
    }
 for (i = 0; i < neg_start; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
        x256 = _mm256_loadu_si256((__m256i*) &c1h_b[deg+j]);
        y256 = _mm256_loadu_si256((__m256i*) &pk_b[j]);
        z256 = _mm256_add_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_b[deg+j], z256);
	}
    }
    for (i = neg_start; i < HR; ++i){
        deg = r_idx[i];
        for ( j = 0; j < LWE_N; j+=32) {
         x256 = _mm256_loadu_si256((__m256i*) &c1h_b[deg+j]);
         y256 = _mm256_loadu_si256((__m256i*) &pk_b[j]);
         z256 = _mm256_sub_epi8(x256, y256);
        _mm256_storeu_si256((__m256i*)&c1h_b[deg+j], z256);
        }
    }

	__m256i two256 = _mm256_set1_epi8(0x02);
	__m256i eight256 = _mm256_set1_epi8(0x08);
	__m256i fc256 = _mm256_set1_epi8(0xfc);
	__m256i f0256 = _mm256_set1_epi8(0xf0);


    for (j = 0; j < LWE_N; j+=32){
        __m256i x256 = _mm256_loadu_si256((__m256i*) &c1h_a[j]);
        __m256i y256 = _mm256_loadu_si256((__m256i*) &c1h_a[LWE_N+j]);
        __m256i z256 = _mm256_sub_epi8(x256, y256);
	z256 = _mm256_add_epi8(z256, two256);
	z256 = _mm256_and_si256(z256, fc256);
        _mm256_storeu_si256((__m256i*)&c1h_a[j], z256);

        x256 = _mm256_loadu_si256((__m256i*) &c1h_b[j]);
        y256 = _mm256_loadu_si256((__m256i*) &c1h_b[LWE_N+j]);
        z256 = _mm256_sub_epi8(x256, y256);
	z256 = _mm256_add_epi8(z256, eight256);
	z256 = _mm256_and_si256(z256, f0256);
        _mm256_storeu_si256((__m256i*)&c1h_b[j], z256);
    }


////////// check c == c_hat ? ////////// 

	for (i = 0; i < LWE_N; ++i) {
		if ((c[i] != ((c1h_a[i] + 0x02) & 0xfc)) || (c[LWE_N + i] != ((c1h_b[i] + 0x08) & 0xf0))){
			// G(c1, u)
			hash_t=calloc((LWE_N+LWE_N+LWE_N/8), sizeof(unsigned char));
			memcpy(hash_t, c, LWE_N+LWE_N);
			memcpy(hash_t+LWE_N+LWE_N, sk+LWE_N, LWE_N/8);
			
			
			sha3_256(shared_k, hash_t, LWE_N+LWE_N+LWE_N/8);

			free(hash_t);
			free(hash);
			return res = 2;
		}
	}

////////// G(c1, delta) ////////// 
	hash_t=calloc((LWE_N+LWE_N+LWE_N/8), sizeof(unsigned char));
	memcpy(hash_t, c, LWE_N+LWE_N);
	memcpy(hash_t+LWE_N+LWE_N, delta_hat, LWE_N/8);
	sha3_256(shared_k, hash_t, LWE_N+LWE_N+LWE_N/8);

	free(hash_t);
	free(hash);

	return res;
}
