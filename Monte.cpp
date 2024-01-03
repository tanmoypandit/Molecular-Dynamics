#include <iostream>
#include <cstdlib>
#include <ctime>
#include<fstream>
#include<cmath>
using namespace std;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (idum <= 0) {
		if (-(idum) < 1) idum=1;
		else idum = -(idum);
		idum2=(idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float fun(float x)
{
    float y;

    y= 1/(1+x*x);

    return(y);



}

int main() {
    
    long seed=-8;
    float fseed;
	int i,j,counter=0,N=100000;
	float a,b,h=0.01,x,y,area;
    cout<<"Enter the lower limit:\n";
    cin>>a;
cout<<"Enter the upper limit:\n";
    cin>>b;
	for(i=0;i<N;i++)
	{
		srand(time(NULL));
		seed=rand();
        fseed=ran2(seed);
    x=a+fseed*(b-a) ;
	srand(time(NULL));
		seed=rand();
        fseed=ran2(seed);    
     y=fun(a)+(fun(b)-fun(a))*fseed;

     if(y<=fun(x))
     {

         counter++;
     }

}

area=counter*(b-a)*fabs(fun(b)-fun(a))/N;
cout<<area;

}
