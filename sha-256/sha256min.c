/**
 * MIT License
 * Copyright (c) 2023 Robert1037
 * sha256min.c v1.0 win10 x64
 * Last modified: 2023-01-12
 **/
#include <stdio.h>
#include <stdlib.h>
#define M_LEN 64 //each M block is 512 bits == 64 bytes
#define Rn(X, n) ((X << (32 - n)) | (X >> n))
#define Sn(X, n) (X >> n)
#define Ch(X, Y, Z) ((X & Y) ^ (~X & Z))
#define Maj(X, Y, Z) ((X & Y) ^ (X & Z) ^ (Y & Z))
#define Sigma_E0(X) (Rn(X, 2) ^ Rn(X, 13) ^ Rn(X, 22)) //Σ0(X) = R2(X) ⊕ R13(X) ⊕ R22(X)
#define Sigma_E1(X) (Rn(X, 6) ^ Rn(X, 11) ^ Rn(X, 25)) //Σ1(X) = R6(X) ⊕ R11(X) ⊕ R25(X)
#define Sigma_o0(X) (Rn(X, 7) ^ Rn(X, 18) ^ Sn(X, 3)) //σ0(X) = R7(X) ⊕ R18(X) ⊕ S3(X)
#define Sigma_o1(X) (Rn(X, 17) ^ Rn(X, 19) ^ Sn(X, 10)) //σ1(X) = R17(X) ⊕ R19(X) ⊕ S10(X)
unsigned int K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};
int main(int argc, char **argv)
{
    FILE *fp = fopen(argv[1], "rb");
    if (!fp)
        return 0;
    register unsigned int *Kp, *W, *M, a, b, c, d, e, f, g, h, T1,
        H0 = 0x6a09e667,  H1 = 0xbb67ae85,  H2 = 0x3c6ef372,  H3 = 0xa54ff53a,
        H4 = 0x510e527f,  H5 = 0x9b05688c,  H6 = 0x1f83d9ab,  H7 = 0x5be0cd19;
    register char *chp, *ch, tmp;
    fseek(fp, 0, SEEK_END);
    a = ftell(fp);
    b = a & 63; // a % 64;
    b = a + ((b < 56) ? (M_LEN - b) : (M_LEN | (M_LEN - b)));
    if (!(ch = (char*)calloc(b + 256, 1))) // 256 == 64 * sizeof(int), it's for W.
        return 0;
    rewind(fp);
    fread(ch, 1, a, fp);
    fclose(fp);
    ch[a] = 0x80;
    chp = ch + b - 12;
    do {
        tmp = chp[0];  chp[0] = chp[3];  chp[3] = tmp;
        tmp = chp[1];  chp[1] = chp[2];  chp[2] = tmp;
        chp -= 4;
    } while (ch <= chp);
    M = (unsigned int*)ch;
    b >>= 2;
    W = M + b;
    W[-2] = a >> 29;
    W[-1] = a << 3;
    do { // 1 -> N blocks of M
        a = H0;  b = H1;  c = H2;  d = H3;
        e = H4;  f = H5;  g = H6;  h = H7;
        Kp = W + 16;
        do { // j : 0 -> 15
            *W = *M;
            W++;  M++;
        } while (Kp != W);
        Kp += 48;
        do { // j : 16 -> 63
            *W = Sigma_o1(W[-2]) + W[-7] + Sigma_o0(W[-15]) + W[-16];
            W++;
        } while (Kp != W);
        Kp = K;  W -= 64;
        do {
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + *Kp; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  Kp++;
        } while (Kp <= &K[63]);
        H0 += a;  H1 += b;  H2 += c;  H3 += d;
        H4 += e;  H5 += f;  H6 += g;  H7 += h;
        W -= 64;
    } while (M != W);
    printf("%.8x%.8x%.8x%.8x%.8x%.8x%.8x%.8x\n", H0, H1, H2, H3, H4, H5, H6, H7);
    return 0;
}