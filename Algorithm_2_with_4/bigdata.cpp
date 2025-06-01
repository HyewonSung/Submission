#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "generic_utils.h"
#include "spqlios/lagrangehalfc_impl.h"
#include <chrono>
#include "64header.h"

#include "poc_64types.h"

Random* global_random = new Random();
const int Globals::k=1;
const int Globals::N=2048;
const int Globals::t=12;
const int Globals::smalln=N;
const int Globals::bgbit=16;
const int Globals::l=5;
const int Globals::basebit=4;
using namespace std;


Globals::Globals(){





        // Offset
        torusDecompOffset = 0;
        for (int i = 0; i <= l; ++i) torusDecompOffset |= (UINT64_C(1)<<(63-i*bgbit));

        // Buffer
        torusDecompBuf = new uint64_t[N];

        //secret keys (and their polynomial interpretation)

       lwekey = new int[N];
        for (int i = 0; i < N; ++i) lwekey[i] = random_bit();

        tlwekey = new IntPolynomiala(N);
        for (int i = 0; i < N; ++i) tlwekey->coefs[i] = lwekey[i];
        // level 2
        in_key = new int64_t[N+1];
        for (int i = 0; i < N; ++i) in_key[i] = lwekey[i];
        in_key[N] = -1;

        //temp = new TLweSample64(N);



  }

  int64_t array_max(Torus64Polynomial *a, int size)
  {


        int64_t max = abs(a->coefs[0]);

        for(int i=1;i< size;i++) {

        int64_t temp= abs(a->coefs[i]);
      if(temp > max) max=temp;
        }

        return max;

  }

 ///////////////////////////////////****** MAIN ******/////////////////////////////////////////


int main(){

//double alpha =pow(2.,-55);
  double alpha1=pow(2.,-55);
//4.26e-13;//pow(2.,-55);


  uint64_t n=1<<20; //DB
  uint64_t pack=1; // data element size: 384 Byte
  uint64_t m=n/pack;

  uint64_t mm=1<<10; // CMUX per 2^10 data elts
  uint64_t tot =m/mm;// after 10 levels, run CMUX with tot(m/2^10) TRLWE samples

  cout << "the number of DB: "<< n << endl;
  cout << "the number of data in one plaintext(#packing): "<< pack << endl;
  cout << "the number of packed DB: "<< m  << endl;
  cout << "each size of data element: "<< "384Byte" << endl;
  cout << "small cmux gate with " << mm << "data elements" << endl;




  int index=rand()%m;
  cout << "-------------index:"<<index<<"------------------" << endl;
   //  rand()%m;//m-1;
  cout << "generating key switching materials (not depending on DB size) may take some time ... " << endl;


  Globals* env = new Globals();



        int32_t k = env->k;
        int32_t N = env->N;
        int32_t l=env->l;

  cout<< "the ciphertext polynomial degree: "<< N << endl;


  int32_t kpl =(k+1)*l;

  int Bgbit=env->bgbit;
  double log_2m=ceil(log(m)/log(2));
  Torus32 log2m = static_cast<int>(log_2m);
  double log_2mm=ceil(log(mm)/log(2));
  Torus32 log2mm = static_cast<int>(log_2mm);





  int MM1=8;
  uint64_t MM=int64_t(1)<<MM1; //plaintext modulus : 2^12

  double log_2t=ceil(log(MM)/log(2));
  uint64_t logt = static_cast<int>(log_2t);

  uint64_t db=logt*N*(m/64);

  uint64_t mdb = max(db,mm);
  uint64_t mct = max(mm,tot);



  TLweSample64 *enc_dataa =new TLweSample64(N);
  TLweSampleFFTa *enc_temp = new_array1<TLweSampleFFTa>(mct,N);
  TLweSampleFFTa *enc_data = new_array1<TLweSampleFFTa>(mct,N);











  Torus64Polynomial *data= new Torus64Polynomial(db);

  int64_t bit[log2m];

  int_to_bin_digit(index,log2m, bit);

        for(uint32_t j=0;j<db;j++)
          data->coefs[j]=(int64_t(1)<<56)*(17*j*rand()%MM); // database as plaintext......?!?!?!?! with a little trick..........OMG!!!!


   Torus64Polynomial *indexdata= new Torus64Polynomial(N);

          for(int32_t j=1;j<N;j++)
          indexdata->coefs[j]=0;

        indexdata->coefs[0]=data->coefs[index];


// client  is making a query.
//////////////////////////////////////////////////////////////////////////
std::cout << "Now client is making a query" << std::endl;

// Step 1: l개의 Torus64Polynomial 동적 할당 (message)
std::cout << "[Step 1] Allocating and initializing message..." << std::endl;
Torus64Polynomial** message = new Torus64Polynomial*[l];
for (int i = 0; i < l; i++) {
    message[i] = new Torus64Polynomial(N);
    for (int j = 0; j < N; j++) message[i]->coefs[j] = 0;
}

// Step 2: bit 정보를 이용해 message 설정
std::cout << "[Step 2] Setting up message coefficients using bit information..." << std::endl;
for (int i = 0; i < l; i++) {
    for (int32_t j = log2m - 1; j >= 0; --j) {
        message[i]->coefs[log2m - 1 - j] = bit[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
    }
}

// Step 3: l개의 TLweSample64 암호문 동적 할당 (cipher)
std::cout << "[Step 3] Encrypting messages into cipher..." << std::endl;
TLweSample64** cipher = new TLweSample64*[l];
for (int i = 0; i < l; i++) {
    cipher[i] = new TLweSample64(N);
    tLwe64Encrypt_debug(cipher[i], message[i], alpha1, env);
}

// Step 4: info_sk 생성 (N 개)
std::cout << "[Step 4] Creating info_sk..." << std::endl;
IntPolynomiala** info_sk = new IntPolynomiala*[N];
for (int iter = 0; iter < N; ++iter) {
    info_sk[iter] = new IntPolynomiala(N);
    for (int j = 0; j < N; ++j) {
        info_sk[iter]->coefs[j] = 0;
    }
    info_sk[iter]->coefs[0] = env->tlwekey->coefs[iter];
}

// Step 5: ksk 생성 (N 개)
std::cout << "[Step 5] Generating key switching keys (ksk)..." << std::endl;
TGswSample64** ksk = new TGswSample64*[N];
for (int i = 0; i < N; i++) {
    ksk[i] = new TGswSample64(env->l, N);
    KSKGen_RGSW_2_debug(ksk[i], info_sk[i], env);
}

// Step 6: cipher_prime 생성
std::cout << "[Step 6] Running packing algorithm2..." << std::endl;

TLweSample64*** cipher_prime = new TLweSample64**[l];
auto start_total = std::chrono::high_resolution_clock::now();

double total_time = 0.0;  // 전체 실행 시간 누적
int count = 0;             // 실행 횟수 카운트

for (int i = 0; i < l; i++) {
    cipher_prime[i] = new TLweSample64*[log2m];
    for (int j = 0; j < log2m; j++) {
        cipher_prime[i][j] = new TLweSample64(N);
        int32_t index = j;

        // 실행 시간 측정
        auto start = std::chrono::high_resolution_clock::now();
        packing_algorithm2(cipher_prime[i][j], const_cast<const TGswSample64**>(ksk), cipher[i], index, env);
        auto end = std::chrono::high_resolution_clock::now();

        // 개별 실행 시간 저장
        std::chrono::duration<double, std::milli> duration = end - start;
        total_time += duration.count();  // 총 시간 누적
        count++;                         // 실행 횟수 증가
    }
}

auto end_total = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> total_duration = end_total - start_total;

// 평균 시간 계산
double avg_time = (count > 0) ? (total_time / count) : 0.0;

std::cout << "Packing (Algorithm2) done." << std::endl;
std::cout << "Total execution time: " << total_duration.count() << " ms" << std::endl;
std::cout << "Average execution time per packing: " << avg_time << " ms" << std::endl;

// Step 7: minus_sk 생성
std::cout << "[Step 7] Generating negative secret key..." << std::endl;
IntPolynomiala* minus_sk = new IntPolynomiala(N);
for (int i = 0; i < N; i++) {
    minus_sk->coefs[i] = -env->tlwekey->coefs[i];
}

// Step 8: convk 생성 및 암호화
std::cout << "[Step 8] Encrypting convk..." << std::endl;
TGswSample64* convk = new TGswSample64(l, N);
tGsw64Encrypt_poly_2(convk, minus_sk, pow(2., -55), env);

// Step 9: extract 생성 및 unpacking_algorithm4 실행
std::cout << "[Step 9] Running unpacking_algorithm4...(RGSW ciphertexts)" << std::endl;
TGswSample64* extract = new_array1<TGswSample64>(log2m, l, N);
TGswSampleFFTa* extFFT = new_array1<TGswSampleFFTa>(log2m, l, N);

auto start0 = std::chrono::high_resolution_clock::now();
for (int j = 0; j < log2m; j++) {
    TLweSample64** rlweInputs = new TLweSample64*[l];
    for (int i = 0; i < l; i++) {
        rlweInputs[i] = cipher_prime[i][j];
    }
    unpacking_algorithm4(&extract[j], const_cast<const TLweSample64**>(rlweInputs), convk, env);
    delete[] rlweInputs;
}
std::cout << "Unpacking (Algorithm4) done" << std::endl;

auto end0 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> execution_time0 = end0 - start0;
std::cout << "Query unpacking step takes: " << execution_time0.count() << " ms" << std::endl;




// 메모리 해제
/*
std::cout << "[Cleanup] Freeing allocated memory..." << std::endl;
for (int i = 0; i < l; i++) {
    delete message[i];
    delete message_scaled_by_N[i];
    delete cipher[i];

    for (int j = 0; j < log2m; j++) {
        delete rot_cipher[i][j];
        delete cipher_prime[i][j];
    }
    delete[] rot_cipher[i];
    delete[] cipher_prime[i];
}
delete[] message;
delete[] message_scaled_by_N;
delete[] cipher;
delete[] rot_cipher;
delete[] cipher_prime;

for (int j = 0; j < log2m; j++) {

    delete ksk[j];
}

delete[] ksk;
delete minus_sk;
delete convk;

std::cout << "Memory cleanup completed!" << std::endl;
 */

////////////////////////////////////////////////////////////////////////////////



  printf("STep2 starts\n\n");




  auto start = std::chrono::high_resolution_clock::now();

  auto start20 = std::chrono::high_resolution_clock::now();
  for (int s=0;s<log2m;s++){
         for (int i=0;i<kpl;i++)
                for (int q=0;q<=k;q++)
                     TorusPolynomial64_ifft_lvl2(&extFFT[s].allsamples[i].a[q],&extract[s].allsamples[i].a[q],  env);
            }

  auto end20 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time20 = end20-start20;
  std::cout <<"query:conversion to FFT: " << execution_time20.count()<<" ms "<<std::endl;

  int32_t j=0;
   int32_t z =0;
   int32_t mk=mm;

  while(z<tot){

          while(j<mm)
         {

         CMuxFFTdb(&enc_temp[j],&extFFT[0],data->coefs[z*mk+j],data->coefs[z*mk+j+1], env);
// CMuxFFT(&enc_temp[j],&extFFT[0],&enc_data[j],&enc_data[j+1], env);


        j=j+2;
       }
        //printf("%d\n",z*mk);

        j=0;

        mm=mm/2;


        for (int32_t i=1; i<log2mm;++i)
        {



          while(j<mm)
           {
           CMuxFFTa(&enc_temp[j],&extFFT[i],&enc_temp[j*2],&enc_temp[(j+1)*2], env);


            j=j+2;



           }


        j=0;
        mm=mm/2;

        }



   for(int q=0;q<=k;++q)

    {for(int c=0;c<N;++c)
     enc_data[z].a[q].values[c]=enc_temp[0].a[q].values[c];
    }
    mm=mk;
    j=0;
    z=z+1;




  }



j=0;



 while(j<tot)
        {

 //      CMuxDecompFFTa(&enc_temp[j],&extFFT[0],decompFFT+(kpl*j), env);
         CMuxFFTa(&enc_temp[j],&extFFT[log2mm],&enc_data[j],&enc_data[j+1], env);


           j=j+2;

          }

        j=0;
        tot=tot/2;


        for (int32_t i=log2mm+1; i<log2m;++i)
        {



          while(j<tot)
           {
           CMuxFFTa(&enc_temp[j],&extFFT[i],&enc_temp[j*2],&enc_temp[(j+1)*2], env);


            j=j+2;



           }


        j=0;
        tot=tot/2;

 }






  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time = end-start;
  std::cout <<"main computation step(running m-1 cmux gates) takes: " << execution_time.count()<<" ms "<<std::endl;




//After the second step the output ciphertext is packed including the desired item.








  for (int q = 0; q <= k; ++q)

                TorusPolynomial64_fft_lvl2(&enc_dataa->a[q],&enc_temp[0].a[q], env);

   Torus64Polynomial *decrypt = new Torus64Polynomial(N);
  auto start11 = std::chrono::high_resolution_clock::now();


        tLwe64Phase_lvl2(decrypt,enc_dataa, env);


  auto end11 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time11 = end11-start11;
  std::cout <<" answer decode takes: "<< execution_time11.count()<<" ms "<<std::endl;

  // 🔹 decrypt와 indexdata의 처음 10개 coefficient 출력
std::cout << "🔹 decrypt coefficients (first 10): ";
for (int i = 0; i < 10; ++i) {
    std::cout << decrypt->coefs[i] << " ";
}
std::cout << std::endl;

std::cout << "🔹 indexdata coefficients (first 10): ";
for (int i = 0; i < 10; ++i) {
    std::cout << indexdata->coefs[i] << " ";
}
std::cout << std::endl;

  Torus64Polynomial *error = new Torus64Polynomial(N);

   for (int i=0;i<N;++i){

        error->coefs[i]=decrypt->coefs[i]-indexdata->coefs[i];
    }
  int64_t temp1 =  array_max(error,N);
        //infinite norm of add error

  double bit_ea = ceil(log(temp1)/log(2));
  cout<< "output noise budget: " << 64-logt-bit_ea-1 << endl;




  int aa=0;
  cout << "Test TLweSymDecrypt " << endl;
    for (int32_t i = 0; i < N; ++i) {
         if (abs(decrypt->coefs[i]-indexdata->coefs[i]) >int64_t(1)<<(63-MM1))

   aa+=1;
    }
        if (aa>0)
    cout<< "decryption failure? = "<< "Yes" << endl;
        else
    cout<< "decryption failure? = "<< "No" << endl;







}