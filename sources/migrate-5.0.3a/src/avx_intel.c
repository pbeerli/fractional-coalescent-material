
#define bool int
#define false 0
#define true 1
#define OSXSAVEFlag (1UL<<27)
#define AVXFlag     ((1UL<<28)|OSXSAVEFlag)
#define FMAFlag     ((1UL<<12)|AVXFlag|OSXSAVEFlag)
#define CLMULFlag   ((1UL<< 1)|AVXFlag|OSXSAVEFlag)
#define VAESFlag    ((1UL<<25)|AVXFlag|OSXSAVEFlag)

#include <stdio.h>
#include <stdint.h>

typedef uint32_t U32;
 
void cpuid(unsigned info, unsigned *eax, unsigned *ebx, unsigned *ecx, unsigned *edx)
{
  __asm__(
	  "cpuid;"                                            /* assembly code */
	  :"=a" (*eax), "=b" (*ebx), "=c" (*ecx), "=d" (*edx) /* outputs */
	  :"a" (info)                                         /* input: info into eax */
	   /* clobbers: none */
	  );
}

bool SimdDetectFeature(U32 idFeature)
{
	unsigned int EAX, EBX, ECX, EDX;
	printf("passed flag %i\n",idFeature);
	cpuid(0, &EAX, &EBX, &ECX, &EDX);
	if (EBX == 0x756e6547  && EDX == 0x49656e69 &&  ECX == 0x6c65746e)
	  printf("intel chip ");
	printf("(EBX == 0x756e6547  && EDX == 0x49656e69 &&  ECX == 0x6c65746e)\n");
	printf("%#010x %#010x %#010x\n",EBX,EDX,ECX);
	cpuid(1, &EAX, &EBX, &ECX, &EDX);
	printf("eax=%i: %#010x %#010x %#010x %#010x [%#010x %#10x]\n", 
	       1, EAX, EBX, ECX, EDX, ECX & idFeature, idFeature);
	if((ECX & idFeature) != idFeature)
	  return false;
	return true;
}


int main()
{
  if(SimdDetectFeature(AVXFlag)==true)
    printf("AVX is available\n");
  else
    printf("AVX is NOT available\n");

  return 0;
}
