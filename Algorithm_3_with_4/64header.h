#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/time.h>
#include "generic_utils.h"
#include "spqlios/lagrangehalfc_impl.h"
#include "poc_64types.h"
#include <chrono>

void torus64PolynomialMultNaive_plain_aux(Torus64* __restrict result, const int* __restrict poly1, const Torus64* __restrict poly2, const int N) {
    const int _2Nm1 = 2*N-1;
    Torus64 ri;
    for (int i=0; i<N; i++) {
        ri=0;
        for (int j=0; j<=i; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
    for (int i=N; i<_2Nm1; i++) {
        ri=0;
        for (int j=i-N+1; j<N; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
}


// A and B of size = size
// R of size = 2*size-1
void Karatsuba64_aux(Torus64* R, const int* A, const Torus64* B, const int size, const char* buf){
    const int h = size / 2;
    const int sm1 = size-1;

    //we stop the karatsuba recursion at h=4, because on my machine,
    //it seems to be optimal
    if (h<=4) {
        torus64PolynomialMultNaive_plain_aux(R, A, B, size);
        return;
    }

    //we split the polynomials in 2
    int* Atemp = (int*) buf; buf += h*sizeof(int);
    Torus64* Btemp = (Torus64*) buf; buf+= h*sizeof(Torus64);
    Torus64* Rtemp = (Torus64*) buf; buf+= size*sizeof(Torus64);
    //Note: in the above line, I have put size instead of sm1 so that buf remains aligned on a power of 2

    for (int i = 0; i < h; ++i) Atemp[i] = A[i] + A[h+i];
    for (int i = 0; i < h; ++i) Btemp[i] = B[i] + B[h+i];

    // Karatsuba recursivly
    Karatsuba64_aux(R, A, B, h, buf); // (R[0],R[2*h-2]), (A[0],A[h-1]), (B[0],B[h-1])
    Karatsuba64_aux(R+size, A+h, B+h, h, buf); // (R[2*h],R[4*h-2]), (A[h],A[2*h-1]), (B[h],B[2*h-1])
    Karatsuba64_aux(Rtemp, Atemp, Btemp, h, buf);
    R[sm1]=0; //this one needs to be set manually
    for (int i = 0; i < sm1; ++i) Rtemp[i] -= R[i] + R[size+i];
    for (int i = 0; i < sm1; ++i) R[h+i] += Rtemp[i];
}




// poly1, poly2 and result are polynomials mod X^N+1
void torus64PolynomialMultKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomiala* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; //that's large enough to store every tmp variables (2*2*N*8)

    // Karatsuba
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N2-1; ++i) result->coefs[i] = R[i] - R[N2+i];
    result->coefs[N2-1] = R[N2-1];

    delete[] R;
    delete[] buf;
}





void torus64PolynomialMultAddKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomiala* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; //that's large enough to store every tmp variables (2*2*N*8)

    // Karatsuba
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N2-1; ++i) result->coefs[i] += R[i] - R[N2+i];
    result->coefs[N2-1] += R[N2-1];

    delete[] R;
    delete[] buf;
}


void tLwe64EncryptZero(TLweSample64* cipher, const double stdev, const Globals* env){
    const int N = env->N;
    const int k= env->k;
    for (int j = 0; j < N; ++j) cipher->b->coefs[j] = random_gaussian64(0, stdev);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j) cipher->a[i].coefs[j] = random_int64();
    }
    for (int i = 0; i < k; ++i) torus64PolynomialMultAddKaratsuba_lvl2(cipher->b, &env->tlwekey[i], &cipher->a[i], env);
}

void tLwe64Encrypt(TLweSample64* cipher,const Torus64Polynomial* mess, const double stdev, const Globals* env){
    const int N = env->N;
    // const int k= env->k;

   tLwe64EncryptZero(cipher, stdev, env);

    for (int32_t j = 0; j < N; ++j)
       cipher->b->coefs[j] += mess->coefs[j];
}

void tLwe64EncryptZero_debug(TLweSample64* cipher, const double stdev, const Globals* env){
    const int N = env->N;
    const int k= env->k;

    
    for (int j = 0; j < N; ++j)
        cipher->b->coefs[j] = random_gaussian64(0, stdev);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j)
            cipher->a[i].coefs[j] = random_int64();
    }

  
    for (int i = 0; i < k; ++i)
        torus64PolynomialMultAddKaratsuba_lvl2(cipher->b, &env->tlwekey[i], &cipher->a[i], env);

   
}

void tLwe64Encrypt_debug(TLweSample64* cipher, const Torus64Polynomial* mess, const double stdev, const Globals* env){
    const int N = env->N;


    tLwe64EncryptZero_debug(cipher, stdev, env);

    for (int32_t j = 0; j < N; ++j)
        cipher->b->coefs[j] += mess->coefs[j];

}

void tLwe64Phase_lvl2(Torus64Polynomial* phase, const TLweSample64* cipher, const Globals* env){
    const int N = env->N;
    const int k= env->k;

    //since we only have AddMult, we compute the opposite of the phase
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -cipher->b->coefs[j];
    }



    for (int i = 0; i < k; ++i) {
        torus64PolynomialMultAddKaratsuba_lvl2(phase, &env->tlwekey[i], &cipher->a[i], env);
    }


    //and we negate the result
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -phase->coefs[j];
    }

}
/*
void circuitPrivKS(TLweSampleFFTa* result, const int u, const LweSample64* x, const Globals* env) {
    const int kslen = env->t;
    const int k=env->k;
    const int N = env->N; // N_lvl1 = n_lvl1
    const int basebit = env->basebit;
    const int base = 1<<basebit;       // base=2 in [CGGI16]
    const int mask = base - 1;
    const int64_t prec_offset = UINT64_C(1)<<(64-(1+basebit*kslen)); //precision ILA: revoir

    // clear result
    for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            result->a[i].values[j] = 0;
        }
    }

    // Private Key Switching
    for (int i = 0; i <= N; ++i) {
        const uint64_t aibar = x->a[i] + prec_offset;

        for (int j = 0; j < kslen; ++j) {
            const uint64_t aij = (aibar>>(64-(j+1)*basebit)) & mask;

            if (aij != 0){
                for (int q = 0; q <= k; ++q) {
                        for(int p=0;p<N;++p){
        result->a[q].values[p] -= env->privKS[u][i][j][aij].a[q].values[p]; }

//result->a[q].coefs[p] -= env->privKS[u][i][j][aij].a[q].coefs[p];
                }
            }
        }
    }

}
*/


void circuitPrivKS(TLweSample64* result, const int u, const LweSample64* x, const Globals* env) {
    const int kslen = env->t;
    const int k=env->k;
    const int N = env->N; // N_lvl1 = n_lvl1
    const int basebit = env->basebit;
    const int base = 1<<basebit;       // base=2 in [CGGI16]
    const int mask = base - 1;
    const int64_t prec_offset = UINT64_C(1)<<(64-(1+basebit*kslen)); //precision ILA: revoir

    // clear result
    for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            result->a[i].coefs[j] = 0;
        }
    }

    // Private Key Switching
    for (int i = 0; i <= N; ++i) {
        const uint64_t aibar = x->a[i] + prec_offset;

        for (int j = 0; j < kslen; ++j) {
            const uint64_t aij = (aibar>>(64-(j+1)*basebit)) & mask;

            if (aij != 0){
                for (int q = 0; q <= k; ++q) {
                        for(int p=0;p<N;++p){
        //result->a[q].values[p] -= env->privKS[u][i][j][aij].a[q].values[p]; }

result->a[q].coefs[p] -= env->privKS[u][i][j][aij].a[q].coefs[p];}
                }
            }
        }
    }

}

/////////////////////////////////////////////////////////////////
//     FFT SPECIFIC SECTION                                    //
/////////////////////////////////////////////////////////////////
/*
LagrangeHalfCPolynomiala* new_LagrangeHalfCPolynomiala_array(int nbelts, int N) {
    return new_array1<LagrangeHalfCPolynomiala>(nbelts,N);
}

void delete_LagrangeHalfCPolynomial_array(int nbelts, LagrangeHalfCPolynomiala* data) {
    delete_array1<LagrangeHalfCPolynomiala>(data);
}
*/

#ifdef USE_FFT
void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomiala* result, const IntPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048.execute_reverse_int(result->values, source->coefs);
}


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomiala* result, const Globals* env) {
    const int N = env->N;
    for (int i=0; i<N; i++)
    result->values[i] = 0;
}
void LagrangeHalfCPolynomialAddTo_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const Globals* env) {
    const int N = env->N;

    for (int i=0; i<N; i++) {
        double ra = a->values[i];
       // double ia = a->values[Ns2+i];
       // double rb = b->values[i];
       // double ib = b->values[Ns2+i];
        result->values[i] += ra;
       // result->values[i+Ns2] += ra+ia;
    }

   // LagrangeHalfCPolynomialAddTo(result->values, a->values, Ns2);
}
void LagrangeHalfCPolynomialSubTo_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const Globals* env) {
    const int N = env->N;

    for (int i=0; i<N; i++) {
        double ra = a->values[i];
       // double ia = a->values[Ns2+i];
       // double rb = b->values[i];
       // double ib = b->values[Ns2+i];
        result->values[i] -= ra;
       // result->values[i+Ns2] += ra+ia;
    }

   // LagrangeHalfCPolynomialAddTo(result->values, a->values, Ns2);
}
void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const LagrangeHalfCPolynomiala* b, const Globals* env) {
    const int Ns2 = env->N/2;
    /*
    for (int i=0; i<Ns2; i++) {
        double ra = a->values[i];
        double ia = a->values[Ns2+i];
        double rb = b->values[i];
        double ib = b->values[Ns2+i];
        result->values[i] += ra*rb-ia*ib;
        result->values[i+Ns2] += ra*ib+ia*rb;
    }
    */
    LagrangeHalfCPolynomialAddMulASM(result->values, a->values, b->values, Ns2);
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048.execute_direct_torus64(result->coefs, source->values);
}

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomiala* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048.execute_reverse_torus64(result->values, source->coefs);
}


#else
//these are fake and slow versions of the FFT, that use Karatsuba instead
void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomiala* result, const IntPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    result->setIntPoly(source, 2048);
}


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomiala* result, const Globals* env) {
    assert(env->N==2048);
    result->setZeroTorus64Poly(2048);
}

void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const LagrangeHalfCPolynomiala* b, const Globals* env) {
    assert(env->N==2048);
    assert(result->torus64Poly!=0);
assert(a->intPoly!=0);
    assert(b->torus64Poly!=0);
    torus64PolynomialMultAddKaratsuba_lvl2(result->torus64Poly, a->intPoly, b->torus64Poly, env);
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    assert(source->torus64Poly!=0);
    for (int i=0; i<2048; i++) result->coefs[i]=source->torus64Poly->coefs[i];
}

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomiala* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N==2048);
    result->setTorus64Poly(source, 2048);
}
#endif



void tLwe64NoiselessTrivial(TLweSample64* cipher, const Torus64Polynomial* mess, const Globals* env){
    const int N = env->N;
    const int k= env->k;

     for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            cipher->a[i].coefs[j] = 0;
        }
    }
 for (int j = 0; j < N; ++j) {

    cipher->b->coefs[j] = mess->coefs[j];
  }
}

void int_to_bin_digit(unsigned int in, int count, int64_t* out)
{
        unsigned int mask =1U << (count-1);
                int k;
                        for (k=0;k< count; k++){
                                                out[k]=(in & mask) ? 1 : 0;
                                                                       in <<=1;


                                                                                        }

}

void tGswTorus64PolynomialDecompH(IntPolynomiala* result, const Torus64Polynomial* sample, const Globals* env){
            const int N = env->N;
            const int l = env->l;
            const int Bgbit = env->bgbit;
            // peut etre tout cela dans le env
            const uint64_t Bg = UINT64_C(1)<<Bgbit;
            const uint64_t mask = Bg-1;
            const int64_t halfBg = Bg/2;
            uint64_t* buf = env->torusDecompBuf;
            const uint64_t offset = env->torusDecompOffset;

    //First, add offset to everyone
    for (int j = 0; j < N; ++j) buf[j]=sample->coefs[j]+offset;

    //then, do the decomposition (in parallel)
    for (int p = 0; p < l; ++p) {
        const int decal = (64-(p+1)*Bgbit);
         int* res_p = result[p].coefs; // res is a int (ok 32)
        for (int j = 0; j < N; ++j) {
            uint64_t temp1 = (buf[j] >> decal) & mask;
            res_p[j] = temp1 - halfBg;
        }
    }
}

        void tGsw64DecompH(IntPolynomiala* result, const TLweSample64* sample, const Globals* env){
    const int l = env->l;
    const int k=env->k;
    for (int i = 0; i <= k; ++i) tGswTorus64PolynomialDecompH(result+(i*l), &sample->a[i], env);
}

        void tGswExternMulToTLwe1(TLweSample64 *accum, const TGswSample64 *sample, const Globals *env) {

    const int32_t N = env->N;
    const int32_t k= env->k;
    const int l = env->l;
    const int32_t kpl = (k+1)*l;
    //TODO: improve this new/delete

 IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
    tGsw64DecompH(decomp, accum, env);
    // tLweClear(accum, par);
 for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            accum->a[i].coefs[j] = 0;
        }
    }

    for (int32_t i = 0; i < kpl; i++) {
      //   tLweAddMulRTo(accum, &dec[i], &sample->all_sample[i], par);
      //  }

     for (int j = 0; j <= k; ++j) torus64PolynomialMultAddKaratsuba_lvl2(accum->a+j, &decomp[i], &sample->allsamples[i].a[j], env);
    }

}


        void CMux(TLweSample64 *result, const TGswSample64 *eps, const TLweSample64 *c0, TLweSample64 *c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

         IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N); Torus64Polynomial* phase = new Torus64Polynomial(N);

         TLweSampleFFTa* accFFT = new TLweSampleFFTa(N);
         TGswSampleFFTa* epsFFT= new TGswSampleFFTa(l,N);

         for (int i=0;i<kpl;i++)
                for (int q=0;q<=k;q++)
                  TorusPolynomial64_ifft_lvl2(&epsFFT->allsamples[i].a[q],&eps->allsamples[i].a[q],  env);

 //     tLweSubTo(c1,c0, params2);//c1=c1-c0
         for (int q = 0; q <= k; ++q)
         for (int j = 0; j < N; ++j) c1->a[q].coefs[j] -= c0->a[q].coefs[j];

 //tGswFFTExternMulToTLwe(c1, eps, params1);//c1=c1*eps

        tGsw64DecompH(decomp, c1, env);
        for (int p = 0; p < kpl; ++p) IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);
        // accFFT initialization
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

        // external product FFT
auto start = std::chrono::high_resolution_clock::now();
         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &epsFFT->allsamples[p].a[q], env);

auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double,std::milli> execution_time = end-start;
std::cout <<"one external product w/o conversion from FFT: " << execution_time.count()<<" ms "<<std::endl;

        // conversion from FFT
        for (int q = 0; q <= k; ++q) TorusPolynomial64_fft_lvl2(c1->a+q,accFFT->a+q, env);


//       tLweAddTo(c1, c0, params2);//c1=c1+c0

        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) c1->a[q].coefs[j] += c0->a[q].coefs[j];
for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].coefs[j] = c1->a[q].coefs[j];

}



void CMuxFFT(TLweSample64 *result, const TGswSampleFFTa *eps, const TLweSample64 *c0, TLweSample64 *c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

        IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N);
        TLweSampleFFTa* accFFT= new TLweSampleFFTa(N);






        for (int q = 0; q <= k; ++q)
         for (int j = 0; j < N; ++j) c1->a[q].coefs[j] -= c0->a[q].coefs[j];



        tGsw64DecompH(decomp, c1, env);


        for (int p = 0; p < kpl; ++p)

        IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);




        // external product FFT


         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);


        // conversion from FFT

     for (int q = 0; q <= k; ++q) TorusPolynomial64_fft_lvl2(c1->a+q,accFFT->a+q, env);



        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) c1->a[q].coefs[j] += c0->a[q].coefs[j];


        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].coefs[j] = c1->a[q].coefs[j];

        delete_array1<IntPolynomiala>(decomp);
        delete_array1<LagrangeHalfCPolynomiala>(decompFFT);
        delete accFFT;



}
void CMuxFFTdb(TLweSampleFFTa *result, const TGswSampleFFTa *eps, const Torus64 c0, const Torus64 c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

        IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N);
        TLweSampleFFTa* accFFT= new TLweSampleFFTa(N);
        TLweSampleFFTa* tempFFT= new TLweSampleFFTa(N);
        TLweSample64 *temp = new TLweSample64(N);
        Torus64 cn= 0;
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

        for (int j = 0; j < N; ++j) {
                temp->a[0].coefs[j]=0;
                temp->a[1].coefs[j] =0;

                                      }
           temp->a[1].coefs[0]=c0;// Set data c0 as a noiseless TRLWE sample

        for (int q = 0; q <= k; ++q)

        TorusPolynomial64_ifft_lvl2(tempFFT->a+q,temp->a+q, env);// convert noiseless TRLWE sample c0 to fft form




         cn= c1-c0;//c1=c1-c0;


         for (int j = 0; j < N; ++j)
                temp->a[1].coefs[0] =cn;


        tGsw64DecompH(decomp, temp , env);


        for (int p = 0; p < kpl; ++p)

        IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);




        // external product FFT


         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);


        // conversion from FFT

   for (int q = k; q <= k; ++q) LagrangeHalfCPolynomialAddTo_lvl2(accFFT->a+q,tempFFT->a+q, env); //c1=c1+c0


        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].values[j] = accFFT->a[q].values[j];


        delete_array1<IntPolynomiala>(decomp);
        delete_array1<LagrangeHalfCPolynomiala>(decompFFT);
        delete accFFT;
        delete tempFFT;
        delete temp;


}


        void CMuxFFTa(TLweSampleFFTa *result, const TGswSampleFFTa *eps, const TLweSampleFFTa *c0, TLweSampleFFTa *c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

        IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N);
        TLweSampleFFTa* accFFT= new TLweSampleFFTa(N);
        TLweSample64* acc = new TLweSample64(N);

        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

//       TGswSampleFFTa* epsFFT= new TGswSampleFFTa(l,N);

//      tLweSubTo(c1,c0, params2);//c1=c1-c0



        for (int q = 0; q <= k; ++q)
        {
        LagrangeHalfCPolynomialSubTo_lvl2(c1->a+q,c0->a+q, env);//c1=c1-c0
        TorusPolynomial64_fft_lvl2(acc->a+q,c1->a+q, env);
        }







        tGsw64DecompH(decomp, acc, env);


        for (int p = 0; p < kpl; ++p)

        IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);




        // external product FFT





         for (int p = 0; p < kpl; ++p)
  for (int q = 0; q <= k; ++q)                                                          LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);



          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddTo_lvl2(accFFT->a+q,c0->a+q, env); //c1=c1+c0


        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].values[j] = accFFT->a[q].values[j];


        delete_array1<IntPolynomiala>(decomp);
        delete_array1<LagrangeHalfCPolynomiala>(decompFFT);
        delete accFFT;
        delete acc;

}

        void CMuxDecompFFT(TLweSample64* c0 ,const TGswSampleFFTa *eps, const  LagrangeHalfCPolynomiala* decompFFT,     const Globals* env){
         int32_t N= env->N;
         int32_t k= env->k;
         int32_t l=env->l;
         int32_t kpl= (k+1)*l;

        TLweSampleFFTa* accFFTa= new TLweSampleFFTa(N);
        TLweSample64 *c11 = new TLweSample64(N);
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFTa->a+q, env);
         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q)          LagrangeHalfCPolynomialAddMul_lvl2(accFFTa->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);

        for (int q = 0; q <= k; ++q) TorusPolynomial64_fft_lvl2(c11->a+q,accFFTa->a+q, env);


        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) c0->a[q].coefs[j] += c11->a[q].coefs[j];


        delete accFFTa;
        delete c11;
}

        void CMuxDecompFFTa(TLweSampleFFTa *c0 ,const TGswSampleFFTa *eps, const  LagrangeHalfCPolynomiala* decompFFT,  const Globals* env){
         int32_t N= env->N;
         int32_t k= env->k;
         int32_t l=env->l;
         int32_t kpl= (k+1)*l;

        TLweSampleFFTa* accFFTa= new TLweSampleFFTa(N);
        //TLweSampleFFTa *c11 = new TLweSampleFFTa(N);
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFTa->a+q, env);

        for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFTa->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);


        for (int q = 0; q <= k; ++q)

        LagrangeHalfCPolynomialAddTo_lvl2(c0->a+q,accFFTa->a+q, env);


        delete accFFTa;

        }












        void shift(TLweSample64 *result, int j,TLweSample64 *sample,int N){

//int N = env->N;
//int k=env->k;


        for (int i=0; i<N; i++){
          if (i+j<N){
  result->a[0].coefs[j+i] = sample->a[0].coefs[i];
   result->a[1].coefs[j+i] = sample->a[1].coefs[i];
          }
          else
           {
   result->a[0].coefs[(i+j)%N]= -sample->a[0].coefs[i];
   result->a[1].coefs[(i+j)%N] = -sample->a[1].coefs[i];

            }


         }
        }

        void tGsw64Encrypt(TGswSample64* cipher, const int mess, const double stdev, const Globals* env){
            const int l = env->l;
            const int Bgbit = env->bgbit;
            const int k=env->k;
         for (int bloc = 0; bloc <= k; ++bloc) {
         for (int i = 0; i < l; ++i) {
            // encryption of 0
            tLwe64EncryptZero(&cipher->samples[bloc][i], stdev, env);
            // add mess*h[i]
            cipher->samples[bloc][i].a[bloc].coefs[0] += mess * (UINT64_C(1) << (64-(i+1)*Bgbit));
                                                }
                                     }
          }


    void tGsw64Encrypt_poly(TGswSample64* cipher, const IntPolynomiala* mess, const double stdev, const Globals* env) {
    const int l = env->l;
    const int Bgbit = env->bgbit;
    const int k = env->k;

    for (int bloc = 0; bloc <= k; ++bloc) {
        for (int i = 0; i < l; ++i) {
           
            tLwe64EncryptZero(&cipher->samples[bloc][i], stdev, env);

            for (int j = 0; j < env->N; ++j) {
                cipher->samples[bloc][i].a[bloc].coefs[j] += mess->coefs[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
            }
        }
    }
}


    void tGsw64Encrypt_poly_2(TGswSample64* cipher, const IntPolynomiala* mess, const double stdev, const Globals* env){
    const int N = env->N;
    const int l = env->l;
    const int Bgbit = env->bgbit;
    const int k = env->k;

    //for (int bloc = 0; bloc <= k; ++bloc) {
        for (int i = 0; i < l; ++i) {
           
            tLwe64EncryptZero_debug(&cipher->samples[0][i], stdev, env);
            tLwe64EncryptZero_debug(&cipher->samples[1][i], stdev, env);


           
            for (int j = 0; j < N; ++j) {
                cipher->samples[0][i].a[0].coefs[j] += mess->coefs[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
                cipher->samples[1][i].a[1].coefs[j] += mess->coefs[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
            }
        }
    //}
}



void tLweExtractLweSampleIndex64(LweSample64* result, const TLweSample64* x, const int32_t index, const Globals *env) {
    const int32_t N = env->N;
    const int32_t k = env->k;
    assert(env->smalln == k*N);

    for (int32_t i=0; i<k; i++) {
      for (int32_t j=0; j<=index; j++)
        result->a[i*N+j] = x->a[i].coefs[index-j];
      for (int32_t j=index+1; j<N; j++)
        result->a[i*N+j] = -x->a[i].coefs[N+index-j];
    }
    result->a[N] = x->a[k].coefs[index];
}



Torus64 lwe64Phase_lvl2(const LweSample64* cipher, const Globals* env) {
    const int n = env->N;
    Torus64 res = *cipher->b;
    for (int i = 0; i < n; ++i) {
        res -= cipher->a[i]*env->lwekey[i];
    }
    return res;
}


void packing_algorithm2(TLweSample64* rlweResult, const TGswSample64** ksk, const TLweSample64* rlweInput, const int32_t index, const Globals* env) {

    const int k = env->k;
    const int N = env->N;
    const int l = env->l;


    // Sample Exraction
    LweSample64* extract_result = new LweSample64(N); 
    tLweExtractLweSampleIndex64(extract_result, rlweInput, index, env);



    Torus64Polynomial* a_poly = new Torus64Polynomial(N);
    for (int i = 0; i < N; ++i) {
        a_poly->coefs[i] = extract_result->a[i];  
    }


    Torus64 b_scalar = extract_result->a[N];  

   
    IntPolynomiala* decomp_a_scalars = new_array1<IntPolynomiala>(l, N);
    tGswTorus64PolynomialDecompH(decomp_a_scalars, a_poly, env);

    
    Torus64Polynomial** temp_result_a = new Torus64Polynomial*[l];
    Torus64Polynomial** temp_result_b = new Torus64Polynomial*[l];
    for (int i = 0; i < l; ++i) {
        temp_result_a[i] = new Torus64Polynomial(N);
        temp_result_b[i] = new Torus64Polynomial(N);
    }

   
    Torus64Polynomial** sum_result_a = new Torus64Polynomial*[N];
    Torus64Polynomial** sum_result_b = new Torus64Polynomial*[N];
    for (int i = 0; i < N; ++i) {
        sum_result_a[i] = new Torus64Polynomial(N);
        sum_result_b[i] = new Torus64Polynomial(N);

       
        for (int j = 0; j < N; ++j) {
            sum_result_a[i]->coefs[j] = 0;
            sum_result_b[i]->coefs[j] = 0;
        }
    }

   
    for (int idx = 0; idx < N; ++idx) {
       
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < N; ++j) {
                temp_result_a[i]->coefs[j] = 0;
                temp_result_b[i]->coefs[j] = 0;
            }
        }

        
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < N; ++j) {
                temp_result_a[i]->coefs[j] = (decomp_a_scalars+i)->coefs[idx] * ksk[idx]->allsamples[i].a[0].coefs[j];
                temp_result_b[i]->coefs[j] = (decomp_a_scalars+i)->coefs[idx] * ksk[idx]->allsamples[i].a[1].coefs[j];
            }
        }


       
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < N; ++j) {
                sum_result_a[idx]->coefs[j] += temp_result_a[i]->coefs[j];
                sum_result_b[idx]->coefs[j] += temp_result_b[i]->coefs[j];
            }
        }
    }

    
    for (int j = 0; j < N; ++j) {
        rlweResult->a[0].coefs[j] = 0;  
        rlweResult->a[1].coefs[j] = 0;  
    }


    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            rlweResult->a[0].coefs[j] -= sum_result_a[i]->coefs[j];
            rlweResult->a[1].coefs[j] -= sum_result_b[i]->coefs[j];
        }
    }

    
    rlweResult->a[1].coefs[0] = b_scalar + rlweResult->a[1].coefs[0];


    
    delete extract_result;
    delete a_poly;

    for (int i = 0; i < l; ++i) {
    delete temp_result_a[i];
    delete temp_result_b[i];
    }
    delete[] temp_result_a;
    delete[] temp_result_b;

    for (int i = 0; i < N; ++i) {
        delete sum_result_a[i];
        delete sum_result_b[i];
    }
    delete[] sum_result_a;
    delete[] sum_result_b;
    delete_array1<IntPolynomiala>(decomp_a_scalars);

}


void KSKGen_RGSW(TGswSample64* ksk, const IntPolynomiala* info_sk, const Globals* env) {
   
    const int l = env->l;            
    const int Bgbit = env->bgbit;  
    const int N = env->N;           
    const double stdev = pow(2., -55); 

    
    std::vector<uint64_t> gadget_vector(l);
    for (int i = 0; i < l; ++i) {
        
        gadget_vector[i] = (UINT64_C(1) << (64 - (i + 1) * Bgbit));
    }

  
    for (int i = 0; i < l; ++i) {  

        tLwe64EncryptZero(&ksk->samples[0][i], stdev, env);

       
        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[0].coefs[j] += info_sk->coefs[j] * gadget_vector[i];
            ksk->samples[0][i].a[1].coefs[j] += info_sk->coefs[j] * gadget_vector[i];
        }
    }

    std::cout << "Key switching key generation completed." << std::endl;
}

void KSKGen_RGSW_2_debug(TGswSample64* ksk, const IntPolynomiala* info_sk, const Globals* env) {
    
    const int l = env->l;            
    const int Bgbit = env->bgbit;     
    const int N = env->N;             
    const double stdev = pow(2., -55); 

    std::vector<uint64_t> gadget_vector(l);
    for (int i = 0; i < l; ++i) {
        gadget_vector[i] = (UINT64_C(1) << (64 - (i + 1) * Bgbit));
    }



    for (int i = 0; i < l; ++i) {  

      
        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].b->coefs[j] = random_gaussian64(0, stdev);
        }

       
        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[0].coefs[j] = random_int64();
        }

       
        torus64PolynomialMultAddKaratsuba_lvl2(ksk->samples[0][i].b, env->tlwekey, &ksk->samples[0][i].a[0], env);


        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[1].coefs[j] += info_sk->coefs[j] * gadget_vector[i];

        }

        
    }

   
}




void unpacking_algorithm4(TGswSample64* result, const TLweSample64** rlweInputs, const TGswSample64* convk, const Globals* env) {

    const int l = env->l;            
    const int N = env->N;           


    for (int k = 0; k <= env->k; ++k) {
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < N; ++j) {
                result->samples[k][i].a[0].coefs[j] = rlweInputs[i]->a[0].coefs[j];
                result->samples[k][i].a[1].coefs[j] = rlweInputs[i]->a[1].coefs[j];
            }
        }
    }

   
    for (int i = 0; i < l; ++i) {
        // C[i] <- ExternalProd(C_i, convk)
        tGswExternMulToTLwe1(&result->samples[0][i], convk, env);

       
    }

    
}



void Alg3_XPowerShift(TLweSample64* result, const TLweSample64* input, int shift, const Globals* env) {
    const int N = env->N; 

    for (int j = 0; j < N; ++j) {
        int new_index = (j * shift) % N;  
        Torus64 sign = ((j * shift) / N) % 2 == 0 ? 1 : -1; 

        result->a[0].coefs[new_index] = sign * input->a[0].coefs[j]; 
        result->b->coefs[new_index] = sign * input->b->coefs[j]; 
    }
}

void Alg3_XPowerShift_sk(IntPolynomiala* result, const IntPolynomiala* input, int shift, const Globals* env) {
    const int N = env->N;

    for (int j = 0; j < N; ++j) {
        int new_index = (j * shift) % N;  
        Torus64 sign = ((j * shift) / N) % 2 == 0 ? 1 : -1;  

        result->coefs[new_index] = sign * input->coefs[j]; 
    }
}

void Alg3_XPowerShift_TorusPoly(Torus64Polynomial* result, const Torus64Polynomial* input, int shift, const Globals* env) {
    const int N = env->N; 

    for (int j = 0; j < N; ++j) {
        int new_index = (j * shift) % N;  
        Torus64 sign = ((j * shift) / N) % 2 == 0 ? 1 : -1;  

        result->coefs[new_index] = sign * input->coefs[j]; 
    }
}

void packing_algorithm3(TLweSample64* result, const TLweSample64* rlweInput, const TGswSample64** ksk, const Globals* env) {
    const int k = env->k;
    const int N = env->N;
    const int l = env->l;
    const int logN = log2(N);


    TGswSampleFFTa** kskFFT = new TGswSampleFFTa*[logN];
    for (int i = 0; i < logN; ++i) {
        kskFFT[i] = new TGswSampleFFTa(l, N);
        for (int p = 0; p < l; ++p)
            for (int q = 0; q <= k; ++q)
                TorusPolynomial64_ifft_lvl2(&kskFFT[i]->allsamples[p].a[q], &ksk[i]->allsamples[p].a[q], env);
    }

    for (int i = 0; i <= k; ++i) {
        for (int j = 0; j < N; ++j) {
            result->a[i].coefs[j] = rlweInput->a[i].coefs[j];
        }
    }


    for (int iter = 0; iter < logN; ++iter) {
        int shift = (N / (1 << iter)) + 1;

        TLweSample64* c_prime = new TLweSample64(N);
        Alg3_XPowerShift(c_prime, result, shift, env);



        Torus64Polynomial* c_prime_a0 = new Torus64Polynomial(N);
        Torus64Polynomial* c_prime_a1 = new Torus64Polynomial(N);

        for (int j = 0; j < N; ++j) {
            c_prime_a0->coefs[j] = c_prime->a[0].coefs[j];
            c_prime_a1->coefs[j] = c_prime->a[1].coefs[j];
        }

 
        IntPolynomiala* decomp_c_prime_a0 = new_array1<IntPolynomiala>(l, N);
        tGswTorus64PolynomialDecompH(decomp_c_prime_a0, c_prime_a0, env);

  
        LagrangeHalfCPolynomiala* decompFFT_c_prime_a0 = new_array1<LagrangeHalfCPolynomiala>(l, N);

        for (int p = 0; p < l; ++p) {
            IntPolynomial_ifft_lvl2(decompFFT_c_prime_a0 + p, decomp_c_prime_a0 + p, env);
        }

        TLweSampleFFTa* accFFT = new TLweSampleFFTa(N);


        for(int q=0; q<=k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

        for (int p = 0; p < l; ++p) {
            LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a, decompFFT_c_prime_a0 + p, &kskFFT[iter]->allsamples[p].a[0], env);
            LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a + 1, decompFFT_c_prime_a0 + p, &kskFFT[iter]->allsamples[p].a[1], env);
        }



        TLweSample64* acc = new TLweSample64(N);
        for (int q = 0; q <= k; ++q)
            TorusPolynomial64_fft_lvl2(acc->a + q, accFFT->a + q, env);


        TLweSample64* temp_result = new TLweSample64(N);
        for (int j = 0; j < N; ++j){
            temp_result->a[0].coefs[j] = -acc->a[0].coefs[j];
            temp_result->b->coefs[j] = c_prime_a1->coefs[j] - acc->b->coefs[j];
        }



        for (int i = 0; i <= k; ++i)
            for (int j = 0; j < N; ++j)
                result->a[i].coefs[j] += temp_result->a[i].coefs[j];



        delete_array1<IntPolynomiala>(decomp_c_prime_a0);
        delete_array1<LagrangeHalfCPolynomiala>(decompFFT_c_prime_a0);

        delete accFFT;
        delete acc;
        delete temp_result;
        delete c_prime_a0;
        delete c_prime_a1;
        delete c_prime;
    }


    for (int i = 0; i < logN; ++i) {
        delete kskFFT[i];
    }
    delete[] kskFFT;

}





void left_shift_by_one(TLweSample64 *output, TLweSample64 *input, int N){
    for(int i=0; i<N; i++){
        output->a[0].coefs[i] = input->a[0].coefs[i+1];
        output->a[1].coefs[i] = input->a[1].coefs[i+1];
    }
    output->a[0].coefs[N-1] = -input->a[0].coefs[0];
    output->a[1].coefs[N-1] = -input->a[1].coefs[0];
}

void left_shift_by_one_poly(IntPolynomiala *output, IntPolynomiala *input, int N){
    for(int i=0; i<N; i++){
        output->coefs[i] = input->coefs[i+1];
    }

}