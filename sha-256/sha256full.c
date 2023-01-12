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
static unsigned int K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};
void Hash(register unsigned int *M, register unsigned int *W)
{
    register unsigned int *Kp, a, b, c, d, e, f, g, h, T1,
        H0 = 0x6a09e667,  H1 = 0xbb67ae85,  H2 = 0x3c6ef372,  H3 = 0xa54ff53a,
        H4 = 0x510e527f,  H5 = 0x9b05688c,  H6 = 0x1f83d9ab,  H7 = 0x5be0cd19;
    do { // 1 -> N blocks of M
        a = H0;  b = H1;  c = H2;  d = H3;
        e = H4;  f = H5;  g = H6;  h = H7;
        Kp = K;
        do { // j : 0 -> 15
            *W = *M;          
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + *Kp; //T1 ← h + Σ1(e) + Ch(e, f, g) + Wj + Kj
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d); //T2 ← Σ0(a) + Maj(a, b, c)
            W++;  Kp++;  M++;
        } while (Kp != &K[16]);
        do { // j : 16 -> 63
            //Wj ← σ1(Wj−2) + Wj−7 + σ0(Wj−15) + Wj−16
            *W = Sigma_o1(W[-2]) + W[-7] + Sigma_o0(W[-15]) + W[-16];
            T1 = h + Sigma_E1(e) + Ch(e, f, g) + *W + *Kp;
            h = g;  g = f;  f = e;  e = d + T1;
            d = c;  c = b;  b = a;  a = T1 + Sigma_E0(b) + Maj(b, c, d);
            W++;  Kp++;
        } while (Kp != &K[63]);
        T1 = h + Sigma_E1(e) + Ch(e, f, g) + *Kp + Sigma_o1(W[-2]) + W[-7] + Sigma_o0(W[-15]) + W[-16];
        H1 += a;  H2 += b;  H3 += c;  H4 += d + T1;
        H5 += e;  H6 += f;  H7 += g;  H0 += T1 + Sigma_E0(a) + Maj(a, b, c);
        W -= 63;
    } while (M != W);
    printf("hash: %.8x%.8x%.8x%.8x%.8x%.8x%.8x%.8x", H0, H1, H2, H3, H4, H5, H6, H7);
}
int main()
{
    printf("*************** sha-256 ***************\nPlease input a file name.\n\
If the file dose not exist,\nthe inputted string will be computed.\n");
    FILE *fp;
    clock_t t1, t2, t3;
    unsigned int *M, size, len;
    register char *name, *ch, tmp1, tmp2;
    for ( ; ; free(M)) {
        printf("(string length < %d, file size < 2 GB)\n(ENTER without input to quit)\n:",
            MAX_LEN - M_LEN);
        for (ch = name = (char*)calloc(MAX_LEN, 1);  (*ch = getchar()) != '\n';  ch++)
            ;
        if (ch == name) //quit
            return 0;
        *ch = '\0';
        t1 = clock(); //convertion start
        if(fp = fopen(name, "rb")) {
            fseek(fp, 0, SEEK_END);
            size = ftell(fp);
            printf("file: %s exists\nlength: %u bytes\n", name, size);
            len = size & 63; // size % 64;
            len = size + ((len < 56) ? (M_LEN - len) : (M_LEN | (M_LEN - len)));
            ch = (char*)calloc(len + 256, 1); // 256 == 64 * sizeof(int), it's for W.
            if (!ch) {
                printf("ERROR: calloc failed!\n");
                exit(0);
            }
            rewind(fp);
            fread(ch, 1, size, fp);
            fclose(fp);
            free(name);
        } else {
            size = ch - name;
            printf("file dose not exist\nlength: %u bytes\n", size);
            len = size & 63; // size % 64;
            len = size + ((len < 56) ? (M_LEN - len) : (M_LEN | (M_LEN - len)));
            ch = (char*)realloc(name, len + 256); // 256 == 64 * sizeof(int), it's for W.
        }
        ch[size] = 0x80;
        name = ch + (size | 3); // ch + 3 + size - (size % 4);
        do {
            tmp1 = *name;  *name = name[-3];  name--;
            tmp2 = *name;  *name = name[-1];  name--;
            *name = tmp2;  name--;
            *name = tmp1;  name--;
        } while (ch < name);
        M = (unsigned int*)ch;
        len >>= 2;
        M[len - 2] = size >> 29;
        M[len - 1] = size << 3;
        t2 = clock(); //computation start
        Hash(M, M + len);
        t3 = clock(); //finish
        printf("\nconvertion  cost time: %ldms\ncomputation cost time: %ldms\ntotal cost time: %ldms\n",
            t2 - t1, t3 - t2, t3 - t1);
    }
}