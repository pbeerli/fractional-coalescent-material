
#include <immintrin.h>

int main()
{
  double xx1[]={1,0.1,0.2,0.3};
  double xx2[]={1,0.5,0.6,0.7};
  double xx3[]={0,0,0,0};
  double p1[4][4] = {{0.99, 0.0009, 0.002, 0.0008},{0.0007,0.99,0.0001,0.0002},\
		     {0.008, 0.003, 0.98,0.0008},{0.0007,0.0005,0.0004,0.97}};
  double p2[4][4] = {{0.89, 0.0009, 0.002, 0.0008},{0.0007,0.89,0.0001,0.0002},\
		     {0.008, 0.003, 0.88,0.0008},{0.0007,0.0005,0.0004,0.87}};
  __m256d x1,x2,x3;
  x1 = _mm256_load_pd(xx1);
  x2 = _mm256_load_pd(xx2);
  __m256d pp11,pp12,pp13,pp14,pp21,pp22,pp23,pp24;;
  pp11 = _mm256_load_pd((double *) (&p1[0])); 
  pp12 = _mm256_load_pd((double *) (&p1[1])); 
  pp13 = _mm256_load_pd((double *) (&p1[2])); 
  pp14 = _mm256_load_pd((double *) (&p1[3])); 
  //
  pp21 = _mm256_load_pd((double *) (&p2[0])); 
  pp22 = _mm256_load_pd((double *) (&p2[1])); 
  pp23 = _mm256_load_pd((double *) (&p2[2])); 
  pp24 = _mm256_load_pd((double *) (&p2[3])); 

  __m256d y0 = _mm256_mul_pd(x1,pp11);
  __m256d y1 = _mm256_mul_pd(x1,pp12);
  __m256d y2 = _mm256_mul_pd(x1,pp13);
  __m256d y3 = _mm256_mul_pd(x1,pp14);
  __m256d z0 = _mm256_mul_pd(x2,pp21);
  __m256d z1 = _mm256_mul_pd(x2,pp22);
  __m256d z2 = _mm256_mul_pd(x2,pp23);
  __m256d z3 = _mm256_mul_pd(x2,pp24);

  __m256d sy01 = _mm256_hadd_pd(y0,y1);//a1,b1,a2,b2
  __m256d sy02 = _mm256_hadd_pd(y2,y3);//c1,d1,c2,d2
  __m256d yblend = _mm256_blend_pd(sy01, sy02, 0b1100);//a1,b1,c2,d2
  __m256d yperm = _mm256_permute2f128_pd(sy01, sy02, 0x21);//a2,b2,c1,d1
  __m256d h1 =  _mm256_add_pd(yperm, yblend);//a,b,c,d

  __m256d sz01 = _mm256_hadd_pd(z0,z1);//a1,b1,a2,b2
  __m256d sz02 = _mm256_hadd_pd(z2,z3);//c1,d1,c2,d2
  __m256d zblend = _mm256_blend_pd(sz01, sz02, 0b1100);//a1,b1,c2,d2
  __m256d zperm = _mm256_permute2f128_pd(sz01, sz02, 0x21);//a2,b2,c1,d1
  __m256d h2 =  _mm256_add_pd(zperm, zblend);//a,b,c,d

  x3 = _mm256_mul_pd(h1,h2);
  // store into the standard pointer location
  _mm256_store_pd(xx3,x3);

  return 0;
}
