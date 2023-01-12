/**
 * MIT License
 * Copyright (c) 2023 Robert1037
 * sha256fast.c v2.0 x64
 * Last modified: 2023-01-12
 **/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define M_LEN 64 //each M block is 512 bits == 64 bytes
#define MAX_LEN (M_LEN << 4)
#define Rn(X, n) ((X << (32 - n)) | (X >> n))
#define Sn(X, n) (X >> n)
#define Ch(X, Y, Z) ((X & Y) ^ (~X & Z))
#define Maj(X, Y, Z) ((X & Y) ^ (X & Z) ^ (Y & Z))
#define Sigma_E0(X) (Rn(X, 2) ^ Rn(X, 13) ^ Rn(X, 22)) //Σ0(X) = R2(X) ⊕ R13(X) ⊕ R22(X)
#define Sigma_E1(X) (Rn(X, 6) ^ Rn(X, 11) ^ Rn(X, 25)) //Σ1(X) = R6(X) ⊕ R11(X) ⊕ R25(X)
#define Sigma_o0(X) (Rn(X, 7) ^ Rn(X, 18) ^ Sn(X, 3)) //σ0(X) = R7(X) ⊕ R18(X) ⊕ S3(X)
#define Sigma_o1(X) (Rn(X, 17) ^ Rn(X, 19) ^ Sn(X, 10)) //σ1(X) = R17(X) ⊕ R19(X) ⊕ S10(X)
int main()
{
    FILE *fp;
    clock_t t1, t2;
    unsigned int K[48] = {
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };
    register unsigned int *Kp, *W, *M, a, b, c, d, e, f, g, h, T1,
        H0 = 0x6a09e667,  H1 = 0xbb67ae85,  H2 = 0x3c6ef372,  H3 = 0xa54ff53a,
        H4 = 0x510e527f,  H5 = 0x9b05688c,  H6 = 0x1f83d9ab,  H7 = 0x5be0cd19;
    register char *name, *ch, tmp1, tmp2;
    printf("filename or string :");
    for (ch = name = (char*)calloc(MAX_LEN, 1);  (*ch = getchar()) != '\n';  ch++)
        ;
    *ch = '\0';
    t1 = clock(); //start
    if(fp = fopen(name, "rb")) {
        fseek(fp, 0, SEEK_END);
        a = ftell(fp);
        b = a & 63; // a % 64;
        b = a + ((b < 56) ? (M_LEN - b) : (M_LEN | (M_LEN - b)));
        if (!(ch = (char*)calloc(b + 254, 1))) // 254 == 62 * sizeof(int), it's for W.
            return 0;
        rewind(fp);
        if ((c = fread(ch, 1, a, fp)) != a)
            return 0;
        fclose(fp);
    } else {
        a = ch - name;
        b = a & 63; // a % 64;
        b = a + ((b < 56) ? (M_LEN - b) : (M_LEN | (M_LEN - b)));
        ch = (char*)realloc(name, b + 254);
    }
    ch[a] = 0x80;
    name = ch + (a | 3); // ch + 3 + a - (a % 4);
    do {
        tmp1 = *name;  *name = name[-3];  name--;
        tmp2 = *name;  *name = name[-1];  name--;
        *name = tmp2;  name--;
        *name = tmp1;  name--;
    } while (ch < name);
    M = (unsigned int*)ch;
    b >>= 2;
    W = M + b;
    W[-2] = a >> 29;
    W[-1] = a << 3;
    do { // 1 -> N blocks of M
        a = H0;  b = H1;  c = H2;  d = H3;
        e = H4;  f = H5;  g = H6;  h = H7;
        //do { // j : 0 -> 15
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x428a2f98; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x71374491; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0xb5c0fbcf; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0xe9b5dba5; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x3956c25b; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x59f111f1; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x923f82a4; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0xab1c5ed5; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0xd807aa98; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x12835b01; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x243185be; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x550c7dc3; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x72be5d74; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x80deb1fe; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0x9bdc06a7; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + 0xc19bf174; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  M++;
        //} while (Kp != &K[16]);
        Kp = K;
        do { // j : 16 -> 63
            *W = Sigma_o1(W[-2]) + W[-7] + Sigma_o0(W[-15]) + W[-16];
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + *Kp;
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d);
            W++;  Kp++;
        } while (Kp != &K[46]);
        T1 = h + Sigma_E1(e) + Ch(e, f, g) + Sigma_o1(W[-2]) + W[-7] + Sigma_o0(W[-15]) + W[-16] + 0xbef9a3f7;
        h = g;  g = f;  f = e;  e = d + T1;
        d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d);
        T1 = h + Sigma_E1(e) + Ch(e, f, g) + Sigma_o1(W[-1]) + W[-6] + Sigma_o0(W[-14]) + W[-15] + 0xc67178f2;
        H1 += a;  H2 += b;  H3 += c;  H4 += d + T1;
        H5 += e;  H6 += f;  H7 += g;  H0 += T1 + Sigma_E0(a) + Maj(a, b, c);
        W -= 62;
    } while (M != W);
    t2 = clock(); //finish
    printf("%.8x%.8x%.8x%.8x%.8x%.8x%.8x%.8x", H0, H1, H2, H3, H4, H5, H6, H7);
    printf("\ntotal cost time: %lf s\n", (double)(t2 - t1)/CLOCKS_PER_SEC);
#ifdef _WIN64
    system("pause");
#endif
    return 0;
}