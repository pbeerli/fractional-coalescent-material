
#include <immintrin.h>

double erfscalar(double z)
{
  t = 1.0 / (1.0 + 0.5 * fabs(z));
  // use Horner's method
  ans = 1.0 - t * exp( -z*z -  1.26551223 +
      t * ( 1.00002368 +
      t * ( 0.37409196 + 
      t * ( 0.09678418 + 
      t * (-0.18628806 + 
      t * ( 0.27886807 + 
      t * (-1.13520398 + 
      t * ( 1.48851587 + 
      t * (-0.82215223 + 
      t * ( 0.17087277))))))))))
    if (z >= 0.0)
      return ans;
    else
      return -ans;
}


double erfpara(double z)
{
  t = 1.0 / (1.0 + 0.5 * fabs(z));
  // use Horner's method
  ans = 1.0 - t * exp( -z*z -  1.26551223 +
      t * ( 1.00002368 +
      t * ( 0.37409196 + 
      t * ( 0.09678418 + 
      t * (-0.18628806 + 
      t * ( 0.27886807 + 
      t * (-1.13520398 + 
      t * ( 1.48851587 + 
      t * (-0.82215223 + 
      t * ( 0.17087277))))))))))
    if (z >= 0.0)
      return ans;
    else
      return -ans;
}


int main()
{
  for (i=0;i<10000;i++)
    {
      double x = erf(y);
      if (x>maxx)
	maxx= x;
    }
  
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
