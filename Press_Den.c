
/*HERE WE SIMULATED THE DYNAMICS OF BROWNIAN PARTICLE IN A HARMONIC OSCILLATOR POTENTIAL
WITH 4/(R^12) REPULSIVE INTERACTION*/

#include<stdlib.h>

#include<math.h>
#include<stdio.h>

#define NOP 784             /*NUMBER OF PARTICLES*/
#define GAMMA 1.0            /*FRICTIONAL CONSTANT*/
#define SIGMA  1.1           /*NOMINAL PARTICLE DIAMETER*/
#define TIME 2.5       /*TOTAL TIME FOR EVOLUTION*/
#define TEMP 1.0             /*TEMPERATURE OF THE HEAT BATH*/
#define Kb   1.0            /*BOLTZMAN CONSTANT*/
#define deltaT  0.00001       /*TIME steps*/
#define BETA 1/(Kb*TEMP)    /*INVERSE OF ENERGY OF THE HEAT BATH*/
#define D 1/(BETA*GAMMA)    /*DIFFUSION CONSTANT*/
#define T_START   1.5            /*TIME FROM WHICH PRESSURE CALCULATION STARTED*/
#define R_CUT 2.5
#define DIMENSION 2.0


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
float gasdev(long idum)
{
	float ran1(long idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

struct particles
{
    float x;             /*X COORDINATES OF THE PARTICLES*/
    float y;             /*Y COORDINATES OF THE PARTICLES*/
    float x_euler;       /*X COORDINATES USED DURING EULER UPDATES*/
    float y_euler;       /*Y COORDINATEDS*/
    int n_number;        /*NUMBER OF NEIGHBOR PARTICLES*/
    float max_disp;      /*MAXIMUM DISPLACEMENT OF THE PARTICLES STARTING FROM A PARTICULAR TIME*/
    float force1_x;      /*FORCE ALONG THE X-DIRECTION BY THE NEIGHBORING PARTICLES BEFORE EULER UPDATE*/
    float force1_y;      /*FORCE ALONG THE Y-DIRECTION BY THE NEIGHBORING PARTICLES BEFORE EULER UPDATE*/
    float force2_x;      /*FORCE CALCULATED FROM COORDINATES OF EULER UPDATE*/
    float force2_y;      /*FORCE CALCULATED FROM COORDINATES OF EULER UPDATE*/
};


float WCA_POT(float r)
{
float f;
f=12.0/pow(r,13.0);
return (f);

}

void main()
{



    FILE *f5;
    f5=fopen("DENSITY_PRESSURE_VARIANCE.DAT","w");
    long seed=1;
    float max=0,min=0, ENERGY=0.0,dr,t,a,h,Lx,Ly,dx,dy;
    int i,j,counter=0,ensemble;
    struct particles arr[NOP];
    float VOLUME,PRESSURE=0.0,DENSITY,VIR=0.0,VAR=0.0,PRESSURE_SQUARED=0.0,force_x=0.0,force_y=0.0,force=0.0;


for(DENSITY=0.9;DENSITY<=.9;DENSITY+=0.01)
    {

       a=sqrt(2/(sqrt(3)*DENSITY));
       h=sqrt(3)*a/2;

       a=floorf(a*100)/100;
       h=floorf(h*100)/100;

        Lx=sqrt(NOP)*a;
        Ly=sqrt(NOP)*h;


        VOLUME=Lx*Ly;


        /*---------------------------------------------------*/
        /* DISTRIBUTING THE PARTICLES IN A TRIANGULAR LATTICE*/
        /*---------------------------------------------------*/
    for(i=0;i<sqrt(NOP);i++)
    {
     for(j=0;j<sqrt(NOP);j++)
     {
        if(i%2==0)
        {
        arr[counter].x=j*a-Lx/2.0;
        arr[counter].y=h*i-Ly/2.0;

        arr[counter].x=floorf( arr[counter].x*100)/100;
        arr[counter].y=floorf( arr[counter].y*100)/100;

        }
       else if(i%2!=0)
        {
        arr[counter].x=(j+0.5)*a-Lx/2.0;
        arr[counter].y=i*h-Ly/2.0;
        arr[counter].x=floorf( arr[counter].x*100)/100;
        arr[counter].y=floorf( arr[counter].y*100)/100;
        }
        counter++;
        }
    }
counter=0;

        /*------------------------------------*/
        /*EVOLUTION OF THE SYSTEM THROUGH TIME*/
        /*------------------------------------*/
counter=0;
for(t=0;t<=TIME;t+=deltaT)
   {
       /*INITIALIZATION*/
 for(i=0;i<NOP;i++)
 {
     arr[i].force1_x=0;
     arr[i].force1_y=0;
     arr[i].force2_x=0;
     arr[i].force2_y=0;
 }
        /*-----------------------------------------------------------------------------------------*/
        /*FINDING THE FORCE ALONG THE X AND Y DIRECTION BY ALL  THE  PARTICLES--BEFORE EULER UPDATE*/
        /*-----------------------------------------------------------------------------------------*/
    for(i=0;i<NOP-1;i++)
    {
    for(j=i+1;j<NOP;j++)
    {
        dx=arr[i].x-arr[j].x;
        dx=dx-Lx*rint(dx/Lx);

        dy=arr[i].y-arr[j].y;
        dy=dy-Ly*rint(dy/Ly);

        dr=sqrt(pow(dx,2)+pow(dy,2));

        force_x=dx*12/pow(dr,14);
        force_y=dy*12/pow(dr,14);

        if(dr<R_CUT)
        {
        arr[i].force1_x+=force_x;
        arr[i].force1_y+=force_y;

        arr[j].force1_x+=-force_x;
        arr[j].force1_y+=-force_y;
        }
       }
    }
        /*---------------------------*/
        /*EULER UPDATE OF COORDINATES*/
        /*---------------------------*/
       for(i=0;i<NOP;i++)
       {

        arr[i].x_euler=arr[i].x+D*BETA*arr[i].force1_x*deltaT+gasdev(seed)*sqrt(2*D*deltaT);
        arr[i].x_euler=arr[i].x_euler-Lx*rint(arr[i].x_euler/Lx);

        arr[i].y_euler=arr[i].y+D*BETA*arr[i].force1_y*deltaT+gasdev(seed)*sqrt(2*D*deltaT);
        arr[i].y_euler=arr[i].y_euler-Ly*rint(arr[i].y_euler/Ly);

       }
        /*-----------------------------------------------------------------------------------------------*/
        /*FINDING THE FORCE ALONG THE X AND Y DIRECTION BY THE NEIGHBORING PARTICLES--AFTER EULER UPDATE*/
        /*-----------------------------------------------------------------------------------------------*/


for(i=0;i<NOP-1;i++)
{
    for(j=i+1;j<NOP;j++)
    {
        dx=arr[i].x_euler-arr[j].x_euler;
        dx=dx-Lx*rint(dx/Lx);

        dy=arr[i].y_euler-arr[j].y_euler;
        dy=dy-Ly*rint(dy/Ly);

        dr=sqrt(pow(dx,2)+pow(dy,2));

        force_x=dx*12/pow(dr,14);
        force_y=dy*12/pow(dr,14);


        if(dr<R_CUT)
        {

        arr[i].force2_x+=force_x;
        arr[i].force2_y+=force_y;

         arr[j].force2_x+=-force_x;
        arr[j].force2_y+=-force_y;
        }

    }

}

        /*----------------------------------*/
        /*STRATONOVICH UPDATE OF COORDINATES*/
        /*----------------------------------*/
for(i=0;i<NOP;i++)
       {

        arr[i].x=arr[i].x+D*BETA*0.5*(arr[i].force1_x+arr[i].force2_x)*deltaT+gasdev(seed)*sqrt(2*D*deltaT);
        arr[i].x=arr[i].x-Lx*rint(arr[i].x/Lx);

        arr[i].y=arr[i].y+D*BETA*0.5*(arr[i].force1_y+arr[i].force2_y)*deltaT+gasdev(seed)*sqrt(2*D*deltaT);
        arr[i].y=arr[i].y-Ly*rint(arr[i].y/Ly);

       }




if(t>=T_START)
{
   /*---------------------*/
   /*CALCULATION OF VIRIAL*/
   /*---------------------*/
   for(i=0;i<NOP-1;i++)
   {
       for(j=i+1;j<NOP;j++)
       {
        dx=arr[i].x-arr[j].x;
        dx=dx-Lx*rint(dx/Lx);

        dy=arr[i].y-arr[j].y;
        dy=dy-Ly*rint(dy/Ly);

        dr=sqrt(pow(dx,2)+pow(dy,2));

        force=12.0/pow(dr,13.0);

        if(dr<R_CUT)
           VIR+=dr*force;
       }
   }
   /*--------------------*/
   /*PRESSURE CALCULATION*/
   /*--------------------*/

   PRESSURE+=(DENSITY/BETA)+VIR/(DIMENSION*VOLUME);
   PRESSURE_SQUARED+=pow((DENSITY/BETA)+VIR/(DIMENSION*VOLUME),2);
   VIR=0.0;
}
}//TIME LOOP
   VAR=sqrt(PRESSURE_SQUARED/((TIME-T_START)/deltaT)-pow(PRESSURE/((TIME-T_START)/deltaT),2));
   fprintf(f5,"%f %f %f\n",DENSITY,PRESSURE/((TIME-T_START)/deltaT),VAR);

VIR=0.0;
PRESSURE=0.0;
PRESSURE_SQUARED=0.0;


    }//DENSITY LOOP

}//MAIN












