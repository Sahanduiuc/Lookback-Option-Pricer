#include <stdio.h>
#include <stdlib.h>

#include <time.h> 
#include <math.h>
#include "pnl/pnl_random.h" 
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_fft.h"






/////////////////////////////procedures from wienerhopf.c//////////////////////////////
static void fft_real (PnlVect * a,int fft_size, int sign)
{
  int    i, opposite;
  double last, second;

  switch (sign)
    {
    case 1 : /* backward */
      second = pnl_vect_get (a, 1);
      opposite = 1;
      for ( i=1 ; i<fft_size - 1 ; i++ )
        {
          pnl_vect_set (a, i, opposite * pnl_vect_get (a, i + 1));
          opposite = - opposite;
        }
      pnl_vect_set (a , fft_size - 1, second);

      pnl_real_ifft_inplace (a->array, fft_size);
      break;
    case -1 : /* forward */
      pnl_real_fft_inplace (a->array, fft_size);
      last = pnl_vect_get (a, fft_size - 1);
      opposite = -1;
      for ( i=fft_size -1 ; i>1 ; i-- )
        {
          pnl_vect_set (a, i, opposite * pnl_vect_get (a, i-1));
          opposite = - opposite;
        }
      pnl_vect_set (a , 1, last);
      break;
    }
}

static int expectation(long int kmax,PnlVect * vv1, PnlVect * bl1)
  
{
  long int i;
  double t1,t2;

  
  vv1->array[0]= vv1->array[0]*bl1->array[0];
  vv1->array[1]= vv1->array[1]*bl1->array[1];

  for(i = 1; i <= kmax-1; i++)
    {
      t1 = vv1->array[2*i];
      t2 = vv1->array[2*i+1];
      vv1->array[2*i]= t1*bl1->array[2*i]-t2*bl1->array[2*i+1];
      vv1->array[2*i+1]= t2*bl1->array[2*i]+t1*bl1->array[2*i+1];
    }
  fft_real(vv1, 2*kmax, 1);
  return 1;
}

//=================================================================
static int findcoefnew(int model, double mu, double sigma, double lm1, double lp1,
		    double num, double nup,
		    double cnum,double cnup, double q, double r1,
		    double T, double h, long int kmax,
		    double er, long int Nt,
		    PnlVect * al1)

{
  long int Nx;
  PnlVect *bl1;
  PnlVect *alin1;
  double sg2;
	double cpl;
	double cml;

  long int k;
  long int i;
  double mod;

  double xi,xip,xim,anp,anm,sxi,nup1,num1;
  long int Nmax=2*kmax;
  double lpm1;

  Nx=Nmax; /*number of space points*/

  /*Memory allocation for space grid*/

  bl1 = pnl_vect_create(Nx+1);
  alin1 = pnl_vect_create(Nx);

  ///////////////////////// this part depends on the model
if(model==1/*TSL*/)
{

  lpm1=q+pow(lp1,nup)*cnup+pow(-lm1,num)*cnum+fabs(mu)/h;

  nup1=nup/2.0;
  num1=num/2.0;
  LET(bl1,0)=q;
  xi=-M_PI/h;
  xip=pow(xi*xi+lp1*lp1,nup1);
  xim=pow(xi*xi+lm1*lm1,num1);
  anp=nup*atan(xi/lp1);
  anm=num*atan(xi/lm1);
  LET(bl1,1)=lpm1-cnup*xip*cos(anp)-cnum*xim*cos(anm)+fabs(mu)/h;
  LET(bl1,2*kmax)=-cnup*xip*sin(anp)-cnum*xim*sin(anm);
  sxi=xi/kmax;
  xi=sxi;

  i=1;
  do
    {
      xip=pow(xi*xi+lp1*lp1,nup1);
      xim=pow(xi*xi+lm1*lm1,num1);
      anp=nup*atan(xi/lp1);
      anm=num*atan(xi/lm1);
      LET(bl1,2*i)=lpm1-cnup*xip*cos(anp)-cnum*xim*cos(anm)-fabs(mu)*cos(xi*h)/h;
      LET(bl1,2*i+1)=-cnup*xip*sin(anp)-cnum*xim*sin(anm)-mu*sin(xi*h)/h;
      i++;
      xi=xi+sxi;
 }while(i<kmax);
}// END TSL
if(model==2/*NIG*/)
{
       lpm1=q-sqrt(-lp1*lm1)*cnup;

	LET(bl1,0)=q;
	xi=-M_PI/h;
	xip=xi*xi-lp1*lm1;
	xim=-(lp1+lm1)*xi;
	anp=0.5*atan(xim/xip);
	LET(bl1,1)=lpm1+cnum*pow(xip*xip+xim*xim, 0.25)*cos(anp);
	LET(bl1,2*kmax)=cnup*pow(xip*xip+xim*xim, 0.25)*sin(anp);
	sxi=xi/kmax;
	xi=sxi;

       i=1;
	do
	{
		xip=xi*xi-lp1*lm1;
		xim=-(lp1+lm1)*xi;
		anp=0.5*atan(xim/xip);
		LET(bl1,2*i)=lpm1+cnum*pow(xip*xip+xim*xim, 0.25)*cos(anp);
		LET(bl1,2*i+1)=cnup*pow(xip*xip+xim*xim, 0.25)*sin(anp);
		i++;
		xi=xi+sxi;
	}while(i<kmax);
}// END NIG
if(model == 3/*VGP*/)
{
       lpm1=q-log(-lp1*lm1)*cnup;

	LET(bl1,0)=q;
	xi=-M_PI/h;
	xip=xi*xi-lp1*lm1;
	xim=-(lp1+lm1)*xi;
	anp=atan(xim/xip);
	LET(bl1,1)=lpm1+cnup*log(xip*xip+xim*xim)/2;
	LET(bl1,2*kmax)=cnup*anp-mu*xi;
	sxi=xi/kmax;
	xi=sxi;

     i=1;
	do
	{
		xip=xi*xi-lp1*lm1;
		xim=-(lp1+lm1)*xi;
		anp=atan(xim/xip);
		LET(bl1,2*i)=lpm1+cnup*log(xip*xip+xim*xim)/2;
		LET(bl1,2*i+1)=cnup*anp-mu*xi;
		i++;
		xi=xi+sxi;
	}while(i<kmax);
}// END VGP
if(model == 4/*KOU*/)
{
    sg2=sigma*sigma/2;
	cpl=cnup*lp1;
	cml=cnum*lm1;

	LET(bl1,0)=q;
	xi=-M_PI/h;
	xip=xi*xi+lp1*lp1;
	xim=xi*xi+lm1*lm1;
	LET(bl1,1)=q+(sg2+cnup/xip+cnum/xim)*xi*xi;
	LET(bl1,2*kmax)=xi*(cpl/xip+cml/xim)-mu*xi;
	sxi=xi/kmax;
	xi=sxi;

	i=1;
	do
	{
		xip=xi*xi+lp1*lp1;
		xim=xi*xi+lm1*lm1;
		LET(bl1,2*i)=q+(sg2+cnup/xip+cnum/xim)*xi*xi;
		LET(bl1,2*i+1)=xi*(cpl/xip+cml/xim)-mu*xi;
		i++;
		xi=xi+sxi;
	}while(i<kmax);

}//   END KOU
if(model == 6/*BS or Heston*/)
{
	// num==var
	double sg2=sigma/2;
	if (sg2>0.00001)
	{	LET(bl1,0)=q;
		xi=-M_PI/h;
		xip=xi*xi;
		LET(bl1,1)=q+sg2*xip;
		LET(bl1,2*kmax)=-mu*xi;
		sxi=xi/kmax;
		xi=sxi;

		i=1;
		do
		{
		xip=xi*xi;
		LET(bl1,2*i)=q+sg2*xip;
		LET(bl1,2*i+1)=-mu*xi;
		i++;
		xi=xi+sxi;
		}while(i<kmax);
	}
	else
		{
		lpm1=q+fabs(mu)/h;
		LET(bl1,0)=q;
		xi=-M_PI/h;
		xip=xi*xi;
		LET(bl1,1)=lpm1+sg2*xip+fabs(mu)/h;
		LET(bl1,2*kmax)=0;
		sxi=xi/kmax;
		xi=sxi;

		i=1;
		do
		{
		xip=xi*xi;
		LET(bl1,2*i)=lpm1+sg2*xip-fabs(mu)*cos(xi*h)/h;
		LET(bl1,2*i+1)=-mu*sin(xi*h)/h;
		i++;
		xi=xi+sxi;
		}while(i<kmax);
	}
}	//END Heston
  //////////////////////////////////////////////////////

  LET(alin1,0)=q/GET(bl1,0);
  mod=GET(bl1,2*kmax)*GET(bl1,2*kmax)+GET(bl1,1)*GET(bl1,1);
  LET(alin1,1)=q*GET(bl1,1)/mod;

  i=1;
  do
    {   mod=GET(bl1,2*i+1)*GET(bl1,2*i+1)+GET(bl1,2*i)*GET(bl1,2*i);
      LET(alin1,2*i+1)=-q*GET(bl1,2*i+1)/mod;
      LET(alin1,2*i)=q*GET(bl1,2*i)/mod;
      i++;
    }while(i<kmax);


  k=0;
  while(k<2*kmax)
    {
      k++;
      LET(al1,k)=GET(alin1,k-1);
 }
  pnl_vect_free(&alin1);
  pnl_vect_free( &bl1);
  return 1;
}
//--------------------------------------------findfactor---------------
static int findfactor(double np, double nm, double h, long int kmax,  double lp1, double lm1, PnlVect *alin1, PnlVect *tp1, PnlVect *tm1)
{
	 PnlVect *tb1, *tbp1, *tbm1, *ltp1, *ltm1;

 long int i;

  double xi, xip, xim, anp, anm, sxi, abp, abm, nup1, num1;
 double t1, t2;
 double angle;
  double mod;

long int Nmax=2*kmax;

  tbp1= pnl_vect_create(Nmax+2);
 ltp1= pnl_vect_create(Nmax+2);
 ltm1= pnl_vect_create(Nmax+2);
 tb1= pnl_vect_create(Nmax+2);
 tbm1= pnl_vect_create(Nmax+2);

nup1 = 0.;
num1 = 0.;

 if (np==nm)
   {
     nup1=np/2.0;
     num1=nm/2.0;
   }
 if (np>nm)
   {
     nup1=np;
	num1=0.0;
   }
 if (np<nm)
   {
     nup1=0.0;
	num1=nm;
   }

     abp=pow(lp1,nup1);
     abm=pow(-lm1,num1);

     LET(ltp1,0)=1;
     LET(ltm1,0)=1;
     xi=-M_PI/h;

     xip=pow(xi*xi+lp1*lp1,nup1/2.0)/abp;
     xim=pow(xi*xi+lm1*lm1,num1/2.0)/abm;
     anp=nup1*atan(xi/lp1);
     anm=num1*atan(xi/lm1);
     LET(ltp1,1)=xip*cos(anp);
     LET(ltm1,1)=xim*cos(anm);
     LET(ltp1,2*kmax)=xip*sin(anp);
     LET(ltm1,2*kmax)=xim*sin(anm);
     sxi=xi/kmax;
     xi=sxi;

    i=1;
 do
 {
  xip=pow(xi*xi+lp1*lp1,nup1/2.0)/abp;
  xim=pow(xi*xi+lm1*lm1,num1/2.0)/abm;
  anp=nup1*atan(xi/lp1);
  anm=num1*atan(xi/lm1);
  LET(ltp1,2*i)=xip*cos(anp);
  LET(ltp1,2*i+1)=xip*sin(anp);
  LET(ltm1,2*i)=xim*cos(anm);
  LET(ltm1,2*i+1)=xim*sin(anm);
  i++;
  xi=xi+sxi;
 }while(i<kmax);

 LET(alin1,1)=GET(alin1,1)*(GET(ltp1,1)*GET(ltm1,1)-GET(ltp1,2*kmax)*GET(ltm1,2*kmax));
 for(i = 1; i <= kmax-1; i++)
   {
     t1 = GET(alin1,2*i);
     t2 = GET(alin1,2*i+1);
     LET(alin1,2*i)= t1*GET(ltm1,2*i)-t2*GET(ltm1,2*i+1);
     LET(alin1,2*i+1)= t2*GET(ltm1,2*i)+t1*GET(ltm1,2*i+1);
   }
 for(i = 1; i <= kmax-1; i++)
   {
        t1 = GET(alin1,2*i);
        t2 = GET(alin1,2*i+1);
        LET(alin1,2*i)= t1*GET(ltp1,2*i)-t2*GET(ltp1,2*i+1);
        LET(alin1,2*i+1)= t2*GET(ltp1,2*i)+t1*GET(ltp1,2*i+1);
   }

 LET(ltp1,0)=1/GET(ltp1,0);
 LET(ltm1,0)=1/GET(ltm1,0);
 LET(ltp1,1)=GET(ltp1,1)/(GET(ltp1,1)*GET(ltp1,1)+GET(ltp1,2*kmax)*GET(ltp1,2*kmax));
 LET(ltm1,1)=GET(ltm1,1)/(GET(ltm1,1)*GET(ltm1,1)+GET(ltm1,2*kmax)*GET(ltm1,2*kmax));

i=1;
 do
   {   mod=GET(ltp1,2*i+1)*GET(ltp1,2*i+1)+GET(ltp1,2*i)*GET(ltp1,2*i);
     LET(ltp1,2*i+1)=-GET(ltp1,2*i+1)/mod;
     LET(ltp1,2*i)=GET(ltp1,2*i)/mod;

     mod=GET(ltm1,2*i+1)*GET(ltm1,2*i+1)+GET(ltm1,2*i)*GET(ltm1,2*i);
     LET(ltm1,2*i+1)=-GET(ltm1,2*i+1)/mod;
     LET(ltm1,2*i)=GET(ltm1,2*i)/mod;
     i++;
   }while(i<kmax);

 LET(tb1,0)=log(GET(alin1,0));
 LET(tb1,1)=log(GET(alin1,1));

 i=1;
 do
   {   mod=GET(alin1,2*i+1)*GET(alin1,2*i+1)+GET(alin1,2*i)*GET(alin1,2*i);
     LET(tb1,2*i)=0.5*log(mod);
     if (GET(alin1,2*i)==0)
       {
	 if (GET(alin1,2*i+1)>0)
	   {LET(tb1,2*i+1)=M_PI/2.0;}
	 else
	   {LET(tb1,2*i+1)=-M_PI/2.0;}
       }
     else
       { angle=atan(GET(alin1,2*i+1)/GET(alin1,2*i));
	 if (GET(alin1,2*i)>0) {LET(tb1,2*i+1)=angle;}
	 if (GET(alin1,2*i)<0) {if (GET(alin1,2*i+1)<0) {LET(tb1,2*i+1)=angle-M_PI;}
	   else {LET(tb1,2*i+1)=angle+M_PI;}
	 }
       }
     i++;
   }while(i<kmax);


 fft_real(tb1, 2*kmax, 1);

  i=1;
 LET(tbp1,0)=0;
 LET(tbm1,0)=0;
 do
 { LET(tbp1,0)=GET(tbp1,0)-GET(tb1,i);
   LET(tbp1,i)=GET(tb1,i);
   LET(tbm1,i)=0;
   i++;
 }while(i<kmax);
 do
   { LET(tbm1,i)=GET(tb1,i);
     LET(tbp1,i)=0;
     LET(tbm1,0)=GET(tbm1,0)-GET(tb1,i);
     i++;
   }while(i<2*kmax);

 fft_real(tbp1, 2*kmax, -1);
 fft_real(tbm1, 2*kmax, -1);

 LET(tp1,0)=exp(GET(tbp1,0));
 LET(tp1,1)=exp(GET(tbp1,1));
 LET(tm1,0)=exp(GET(tbm1,0));
 LET(tm1,1)=exp(GET(tbm1,1));

 i=1;
 do
   {   mod=exp(GET(tbp1,2*i));
     LET(tp1,2*i)=mod*cos(GET(tbp1,2*i+1));
     LET(tp1,2*i+1)=mod*sin(GET(tbp1,2*i+1));
     mod=exp(GET(tbm1,2*i));
     LET(tm1,2*i)=mod*cos(GET(tbm1,2*i+1));
     LET(tm1,2*i+1)=mod*sin(GET(tbm1,2*i+1));
     i++;
   }while(i<kmax);

 LET(tp1,0)= GET(tp1,0)*GET(ltp1,0);
 LET(tp1,1)= GET(tp1,1)*GET(ltp1,1);

 LET(tm1,0)= GET(tm1,0)*GET(ltm1,0);
 LET(tm1,1)= GET(tm1,1)*GET(ltm1,1);

  for(i = 1; i <= kmax-1; i++)
    {
      t1 = GET(tm1,2*i);
      t2 = GET(tm1,2*i+1);
      LET(tm1,2*i)= t1*GET(ltm1,2*i)-t2*GET(ltm1,2*i+1);
      LET(tm1,2*i+1)= t2*GET(ltm1,2*i)+t1*GET(ltm1,2*i+1);
    }
  for(i = 1; i <= kmax-1; i++)
    {
      t1 = GET(tp1,2*i);
      t2 = GET(tp1,2*i+1);
      LET(tp1,2*i)= t1*GET(ltp1,2*i)-t2*GET(ltp1,2*i+1);
      LET(tp1,2*i+1)= t2*GET(ltp1,2*i)+t1*GET(ltp1,2*i+1);
    }
   for(i = 1; i <= kmax-1; i++)
    {
      t1 = GET(tp1,2*i);
      t2 = GET(tp1,2*i+1);
      LET(tp1,2*i)= t1*GET(ltp1,2*i)-t2*GET(ltp1,2*i+1);
      LET(tp1,2*i+1)= t2*GET(ltp1,2*i)+t1*GET(ltp1,2*i+1);
    }

 pnl_vect_free(&tbp1);
 pnl_vect_free(&ltp1);
 pnl_vect_free(&ltm1);
 pnl_vect_free(&tb1);
 pnl_vect_free(&tbm1);

 return 1;
}//end findfactor

//=================================================================
static int findcoef(int model, double mu, double sigma, double lm1, double lp1,
		    double num, double nup,
		    double cnum,double cnup, double q, double r1,
		    double T, double h, long int kmax,
		    double er, long int Nt,
		    PnlVect * al1)
  
{
  PnlVect *bl1;
  long int k;
  long int i;
  double mod;
  PnlVect *alin1;
  double xi,xip,xim,anp,anm,sxi,nup1,num1;
  double lpm1;
 
  double sg2;
  double cpl;     //  !!!!!!
  double cml;

  /*Memory allocation for space grid*/
  
  bl1=pnl_vect_create(2*kmax+1);
 
  ///////////////////////// this part depends on model
if(model==1/*TSL*/)
{
  
  lpm1=q+pow(lp1,nup)*cnup+pow(-lm1,num)*cnum;  
  
  nup1=nup/2;
  num1=num/2;
  bl1->array[0]=q;
  xi=-M_PI/h;
  xip=pow(xi*xi+lp1*lp1,nup1);
  xim=pow(xi*xi+lm1*lm1,num1);
  anp=nup*atan(xi/lp1);
  anm=num*atan(xi/lm1);
  bl1->array[1]=lpm1-cnup*xip*cos(anp)-cnum*xim*cos(anm); 
  bl1->array[2*kmax]=-cnup*xip*sin(anp)-cnum*xim*sin(anm)-mu*xi;;
  sxi=xi/kmax;
  xi=sxi;
  
  i=1;
  do
    {  
      xip=pow(xi*xi+lp1*lp1,nup1);
      xim=pow(xi*xi+lm1*lm1,num1);
      anp=nup*atan(xi/lp1);
      anm=num*atan(xi/lm1);
      bl1->array[2*i]=lpm1-cnup*xip*cos(anp)-cnum*xim*cos(anm);
      bl1->array[2*i+1]=-cnup*xip*sin(anp)-cnum*xim*sin(anm)-mu*xi;;
      i++;
      xi=xi+sxi;
 }while(i<kmax);
}// END TSL
if(model==2/*NIG*/)
{  
       lpm1=q-sqrt(-lp1*lm1)*cnup;
   
	bl1->array[0]=q; 
	xi=-M_PI/h;
	xip=xi*xi-lp1*lm1;
	xim=-(lp1+lm1)*xi;
	anp=0.5*atan(xim/xip);
	bl1->array[1]=lpm1+cnum*pow(xip*xip+xim*xim, 0.25)*cos(anp);
	bl1->array[2*kmax]=cnup*pow(xip*xip+xim*xim, 0.25)*sin(anp);
	sxi=xi/kmax;
	xi=sxi;
	
       i=1;
	do
	{   
		xip=xi*xi-lp1*lm1;
		xim=-(lp1+lm1)*xi;
		anp=0.5*atan(xim/xip);
		bl1->array[2*i]=lpm1+cnum*pow(xip*xip+xim*xim, 0.25)*cos(anp);	
		bl1->array[2*i+1]=cnup*pow(xip*xip+xim*xim, 0.25)*sin(anp);
		i++;
		xi=xi+sxi;
	}while(i<kmax);
}// END NIG
if(model == 3/*VGP*/)
{
       lpm1=q-log(-lp1*lm1)*cnup;
   
	bl1->array[0]=q; 
	xi=-M_PI/h;
	xip=xi*xi-lp1*lm1;
	xim=-(lp1+lm1)*xi;
	anp=atan(xim/xip);
	bl1->array[1]=lpm1+cnup*log(xip*xip+xim*xim)/2;
	bl1->array[2*kmax]=cnup*anp-mu*xi;
	sxi=xi/kmax;
	xi=sxi;
	
    i=1;
	do
	{   
		xip=xi*xi-lp1*lm1;
		xim=-(lp1+lm1)*xi;
		anp=atan(xim/xip);
		bl1->array[2*i]=lpm1+cnup*log(xip*xip+xim*xim)/2;	
		bl1->array[2*i+1]=cnup*anp-mu*xi;
		i++;
		xi=xi+sxi;
	}while(i<kmax);
}// END VGP
if(model == 4/*KOU*/)
{
	// nup==sigma
     sg2=sigma*sigma/2;
	cpl=cnup*lp1;     //  !!!!!!
	cml=cnum*lm1;
    
	//lpm1=q+pow(lp1,nup)*cnup+pow(-lm1,num)*cnum;

	bl1->array[0]=q; 
	xi=-M_PI/h;
	xip=xi*xi+lp1*lp1;
	xim=xi*xi+lm1*lm1;
	bl1->array[1]=q+(sg2+cnup/xip+cnum/xim)*xi*xi;
	bl1->array[2*kmax]=xi*(cpl/xip+cml/xim)-mu*xi;
	sxi=xi/kmax;
	xi=sxi;
	
	i=1;
	do
	{   
		xip=xi*xi+lp1*lp1;
		xim=xi*xi+lm1*lm1;
		bl1->array[2*i]=q+(sg2+cnup/xip+cnum/xim)*xi*xi;	
		bl1->array[2*i+1]=xi*(cpl/xip+cml/xim)-mu*xi;
		i++;
		xi=xi+sxi;
	}while(i<kmax);
	
}//   END KOU
  //////////////////////////////////////////////////////
  
  alin1=pnl_vect_create(2*kmax);
  
  alin1->array[0]=q/bl1->array[0];
  mod=bl1->array[2*kmax]*bl1->array[2*kmax]+bl1->array[1]*bl1->array[1];
  alin1->array[1]=q*bl1->array[1]/mod;
  //alin1->array[1]=q/bl1->array[1];
  
  i=1;
  do
    {   mod=bl1->array[2*i+1]*bl1->array[2*i+1]+bl1->array[2*i]*bl1->array[2*i];
      alin1->array[2*i+1]=-q*bl1->array[2*i+1]/mod;
      alin1->array[2*i]=q*bl1->array[2*i]/mod;
      i++;
    }while(i<kmax);
  

  k=0;
  while(k<2*kmax)
    {
      k++;
      al1->array[k]=alin1->array[k-1];
 }
  pnl_vect_free(&alin1);
  pnl_vect_free(&bl1);
  return 1;
} 

/////////////////////////////end procedures from wienerhopf.c//////////////////////////////

static int froot(long int kmax,PnlVect * F, double Fx,long int *k0,long int *k1)
  {
			long int k;
           *k0=0; *k1=kmax;
		   if(Fx<GET(F,0)){*k1=*k0;}
		   while (*k1>*k0+1)
		   {k=int(*k0+*k1)/2;
		   if (Fx>GET(F,k)){*k0=k;}
		   else {*k1=k;}
		   }
		   
  return 1;
}
static int frootm(long int kmax,PnlVect * F, double Fx,long int *k0,long int *k1)
  {
			long int k;
           *k0=0; *k1=kmax;
		   if(Fx>=GET(F,0)){*k1=*k0;}
		   while (*k1>*k0+1)
		   {k=int(*k0+*k1)/2;
		   if (Fx<GET(F,k)){*k0=k;}
		   else {*k1=k;}
		   }
		   
  return 1;
}

int TSL_Mc_lookbackfixed_WHMC(int ifCall, double S0, double K,double
T,double r,double divid, double lm1, double lp1, double num,double nup, double cm,double cp, double h, double er,
int generator,int n_paths,double *ptPrice, double *priceError)
{
    
	double cnup, cnum, lpnu, lmnu, mu, qu, q, om, kx;
	double *payoff, *sum_payoff, *sum_square_payoff, *var_payoff, *x;
    double discount;
	long int i, kmax, Nmax, k, k0, k1, kmin;
    int j,nn, n, Nx;
	double Vn, Sn; //
	double dt; // time variables
	PnlVect *coef;
	PnlVect  *vv0,*vv1,*y,*eom,*Fp,*Fm,*tp1,*tm1, *alin1,*al1;
	PnlRng *rng;  
 	{om=lm1<-2. ? 2. : (-lm1+1.)/2.;  }
 	om=0.;
	cnup=cp*pnl_tgamma(-nup);
	cnum=cm*pnl_tgamma(-num);

	lpnu=exp(nup*log(lp1));
	lmnu=exp(num*log(-lm1));

 	mu= r - divid + cnup*(lpnu-exp(nup*log(lp1+1.0))) + cnum*(lmnu-exp(num*log(-lm1-1.0)));// drift
	qu = (pow(lp1,nup) - pow(lp1+om,nup))*cnup + (pow(-lm1,num)-pow(-lm1-om,num))*cnum-mu*om;
	nn=14;
  	Vn=0;
	Nx=10;
	discount=exp(-r*T);
  
    sum_payoff=0;
    sum_square_payoff=0;

   	rng= pnl_rng_create(generator); 
    /* rng must be initialized. */  
    pnl_rng_sseed (rng, 0); 
 


    j=64;
    kx=er*log(2.0)/h;
    while(j<kx)
    {
      j=j*2;
    }
    kmax=j;
    Nmax=2*kmax;//  the number of the space points for the cdfs
    
 //Memory allocation for space grid
	x=new double[Nx];
	payoff=new double[Nx];
	sum_payoff=new double[Nx];
	sum_square_payoff=new double[Nx];
	var_payoff=new double[Nx];
	coef=pnl_vect_create_from_zero(nn+1);  
 y=pnl_vect_create_from_zero(Nmax+1);//space grid points
  eom=pnl_vect_create_from_zero(Nmax+1);//space grid points
 vv0=pnl_vect_create_from_zero(Nmax+1);//payoff
 vv1=pnl_vect_create_from_zero(Nmax+1);//intermediate function
 Fp=pnl_vect_create_from_zero(kmax+1);// CDF for the Plus WH-factor
 Fm=pnl_vect_create_from_zero(kmax+1);// CDF for the Plus WH-factor
 al1=pnl_vect_create_from_zero(2*kmax+2);
 alin1=pnl_vect_create_from_zero(2*kmax+2);
 tp1=pnl_vect_create_from_zero(2*kmax+2); 
 tm1=pnl_vect_create_from_zero(2*kmax+2);

 lp1+=om; // parameters shift to calculate CDFs for Plus and Minus WH-factors
 lm1+=om; //
 
 ///////////////////cdfs calculations///////////////////////////////
   for(k=1;k<=Nx;k++) x[k]=0.02*k;

    if(ifCall==0)
    for(j=1;j<=kmax;j++)
      {
	y->array[j]=(j-kmax)*h;
	vv0->array[j-1]=0;
	eom->array[j]=exp(-om*y->array[j]);
      }
  else
    for(j=kmax;j<=Nmax;j++)
      {
	y->array[j]=(j-kmax)*h;
	eom->array[j]=exp(-om*y->array[j]);
	vv0->array[j-1]=0;
      }
  
  if(ifCall==0) // IF PUT - lookback fixed strike put 
    {
      j=kmax+1;
      do
	    {   y->array[j]=(j-kmax)*h;
	      eom->array[j]=exp(-om*y->array[j]);
	      vv0->array[j-1]=-eom->array[j];
	      j++;
	    }while(j<=Nmax);
	
  }// END IF PUT
  else // IF CALL  lookback fixed strike call
    {
      j=kmax-1;
      do
	    {   y->array[j]=(j-kmax)*h;
	      eom->array[j]=exp(-om*y->array[j]);
	  	     
		vv0->array[j-1]=-eom->array[j];
	      j--;
	    }while(j>0);
	}
    // END IF CALL

fft_real(vv0, 2*kmax, -1);

  for(i=0;i<kmax;i++)
  {
	  Fm->array[i]=0;
	  Fp->array[i]=0;}

//----------------------------------------------weghts------------------------------------- gs

coef->array[1]=0.00277778;
coef->array[2]=-6.40277778;
coef->array[3]=924.05000000;
coef->array[4]=-34597.92777778;
coef->array[5]=540321.11111111;
coef->array[6]=-4398346.36666666;
coef->array[7]=21087591.77777770;
coef->array[8]=-63944913.04444440;
coef->array[9]=127597579.55000000;
coef->array[10]=-170137188.08333300;
coef->array[11]=150327467.03333300;
coef->array[12]=-84592161.50000000;
coef->array[13]=27478884.76666660;
coef->array[14]=-3925554.96666666;
//----------------------------------------------------------------------------------------- gs


//--------------------------------------------------------------loop in weighted terms------------------------
for(n=1; n<=14;) // for(n=1; n<=nn;)
{	
 /*Time step*/
 dt=T/(log(2.0)*n); //the Gaver-Stehfest algorithm parameter 
 q=qu+1/dt;
 findcoef(1, mu, 0, lm1, lp1, num, nup, cnum, cnup, q, r,T, h, kmax, er, 14, al1);
 //===========================================



 for(i=0;i<2*kmax;i++)
 {alin1->array[i]=al1->array[i+1]; 
 }
 
/*TSL*/ 
findfactor(0.0, 0.0, h,  kmax, lp1, lm1, alin1, tp1, tm1);


  //======= EXCHANGE tp AND tm ARRAYS FOR UP-OUT
  if(ifCall==1) ///////check
    {
      alin1 = tp1;
      tp1 = tm1;
      tm1 = alin1;
    }
  //=========================================
  
  for(i=0;i<Nmax;i++)
    {
      tp1->array[i]=tp1->array[i]/q;
      vv1->array[i]=vv0->array[i];
    }
         
  expectation(kmax,vv1,tp1);
  
  for(i=1;i<Nmax+1;i++)                                            
  {
       vv1->array[i-1]=vv1->array[i-1]/eom->array[i]; 
  }

 
  if(ifCall==0)// IF PUT     
  {for(i=0;i<kmax;i++)
  {Fm->array[i]=Fm->array[i]+coef->array[n]*vv1->array[kmax+i-1];}
   }
  else //IF Call 
  {for(i=0;i<kmax;i++)
  {Fp->array[i]=Fp->array[i]+coef->array[n]*vv1->array[kmax+i-1];}
  }
   n++;
}//=====================================END LOOP in n
    kmin=kmax-1;
    if(ifCall==0)// IF PUT     
	{for(i=kmax-1;i>0;i--){Fm->array[i]=Fm->array[i]*log(2.0)/T+1.;if(Fm->array[i]<1./(double)n_paths){Fm->array[i]=0;kmin=i;}}
    Fm->array[0]=1.;
	}
	else
    for(i=0;i<kmax;i++){Fp->array[i]=Fp->array[i]*log(2.0)/T+1.;if(Fp->array[i]<1./(double)n_paths){Fp->array[i]=0;}}

/*
if(ifCall==0)// IF PUT
{
	
	FILE* mytemp=fopen("Fm.dat", "w");
    if (mytemp == NULL)
    {
        perror("Could not open file");
        return 0;
    }
	
	for(i=0;i<kmax;i++)
	{  
		fprintf(mytemp, "%lf %lf\n",y->array[kmax+i], Fm->array[i]);
	//mytemp<<<<" "<<Fm->array[i]<<endl;
	}

	fclose(mytemp);
}
if(ifCall==1)// IF CALL
{
	
	   
	//mytemp<<y->array[kmax+i]<<" "<<Fp->array[i]<<endl;
	
}

*/

 

/////////////end cdfs calculations////////////////
  
 
////////////Monte Carlo simulation//////////////////
 
           
/*	   if(ifCall==1)//call
        {
         for(i=0;i<n_paths;i++)
           Sn=pnl_rng_uni(rng);
		   froot(kmax, Fp, Sn, &k0, &k1);
		 
		   // printf("factors + - = %d %d %d\n", k0,k1,k1-k0);
		   // printf("factors + - = %f %f %f\n", GET(Fp,k0)-Sn,GET(Fp,k1)-Sn,Sn);
			if (k1==0){Sn=0.0;}
			else{Sn=GET(y,k0)+h*(Sn-GET(Fp,k0))/(GET(Fp,k1)-GET(Fp,k0));}
			//printf("factors + - = %f %f %f\n", GET(y,k0),GET(y,k1),Sn);
		   In=pnl_rng_uni(rng1);
		   froot(kmax,Fm, In, &k0, &k1);
		  // printf("factors + - = %d %d %d\n", k0,k1,k1-k0);
		   //printf("factors + - = %f %f %f\n", GET(Fm,k0)-In,GET(Fm,k1)-In,In);
		   In=GET(y,k0+kmax+1)+h*(In-GET(Fm,k0))/(GET(Fm,k1)-GET(Fm,k0));
		//   printf("factors + - = %f %f %f\n", GET(y,k0+kmax+1),GET(y,k1+kmax+1),In);
		   Vn=Vn+Sn+In;
		   Jn=MIN(J0,MIN(Vn,V0+In));
		   tn=tn+dt;
		   if (Jn<=log_l_S0) {j=n_points+2;}
		  }
		 // printf("factors + - = %f %f %f\n", Sn, Jn, Vn);
          payoff=discount*(S0*exp(Vn)-K)*(S0*exp(Vn)>K)*(Jn>log_l_S0)+exp(-r*tn)*rebate*(Jn<=log_l_S0);
          sum_payoff+=payoff;
          sum_square_payoff+=payoff*payoff;
         }
         var_payoff=(sum_square_payoff-sum_payoff*sum_payoff/n_paths)/(n_paths-1);
         *ptPrice=sum_payoff/n_paths;
         *priceError=1.96*sqrt(var_payoff)/sqrt((double)n_paths);
        }*/
        if(ifCall==0)//put
        { for(k=1;k<=Nx;k++)
			{
				sum_payoff[k]=0;
				sum_square_payoff[k]=0;
			}
         for(i=0;i<n_paths;i++)///////////////for paralleling
		 {
           Sn=pnl_rng_uni(rng);
		   frootm(kmin, Fm, Sn, &k0, &k1);
         if (k1==0){Vn=0.0;}
			else{Vn=-GET(y,kmax+k0)-h*(Sn-GET(Fm,k0))/(GET(Fm,k1)-GET(Fm,k0));} 
		    for(k=1;k<=Nx;k++)
			{
				payoff[k]=MAX(K-S0*exp(Vn+x[k]),0);
				sum_payoff[k]+=payoff[k];
				sum_square_payoff[k]+=payoff[k]*payoff[k];
		    }
         }///////////end for paralleling
		 FILE* mytemp=fopen("MCprices.dat", "w");
			if (mytemp == NULL)
			{
			perror("Could not open file");
			return 0;
			}
				
			 for(k=1;k<=Nx;k++)
		 {
         var_payoff[k]=(sum_square_payoff[k]-sum_payoff[k]*sum_payoff[k]/n_paths)/(n_paths-1);
         *ptPrice=discount*sum_payoff[k]/n_paths;
		 *priceError=1.96*discount*sqrt(var_payoff[k])/sqrt((double)n_paths);
		 fprintf(mytemp, "%lf %.8f %.8f\n", x[k], *ptPrice, *priceError);
		 }
		 fclose(mytemp);
        }
       
   
       pnl_rng_free (&rng); 
  /* Frees the generator */  
   return OK;
}



int main ()  
{  
  
  double  S0, Mu, Divid, R,  r; // general parameters
//  double Sigma, Lambda,  P; // Kou parameters
  double LambdaPlus, LambdaMinus; // Kou and TSL parameters
  double AlphaPlus, AlphaMinus, CPlus, CMinus; //TSL parameters
  double T, Strike;
  double h, er;
  int ifCall;
  int RandGen, M;  
  long int N;
  double ptPrice, priceError;
  double start_time, alltime;

  // Option parameters
  S0=100.0;
  T=0.1;
  Strike=100.0;
  ifCall=0; // 0 - Put, 1 - Call
 // End =Option parameters
  //Kou model parameters
 /* Mu=0.0;
  Divid=0.0;
  R=5.12711;
  Sigma=0.3;
  Lambda=0.0;//=7.0;
  LambdaPlus=50.0;
  LambdaMinus=25.0;
  P=0.6; */           
  // End =Kou model parameters
  //TSL model parameters
  Mu=0.0;
  Divid=0.0;
  R=(exp(0.04)-1)*100;
  LambdaPlus=10.0;
  LambdaMinus=3.0;
  AlphaPlus=1.2;
  AlphaMinus=1.2;
  CPlus=0.2395;
  CMinus=0.2395;

  /*Mu=0.0;
  Divid=0.0;
  R=(exp(0.04)-1)*100;
  LambdaPlus=30.0;
  LambdaMinus=25.0;
  AlphaPlus=1.4;
  AlphaMinus=1.4;
  CPlus=0.146568131101155;
  CMinus=0.146568131101155;*/
  // End =TSL model parameters
  // Method parameters
  RandGen=PNL_RNG_MERSENNE;  // generator
  M = 400;  // number of time points
  N = 100000; // number of sample paths
  h=0.001; // space step for CDF
  er=4.0; // scale range in the space
  // End =Method parameters
  // Results init//
  ptPrice=0.0;
  priceError=0.0;
   // End Results init
  r=log(1.+R/100.);

 start_time = ((double)clock())/(double)CLOCKS_PER_SEC;
  TSL_Mc_lookbackfixed_WHMC(ifCall, S0, Strike, T, r, Divid, -LambdaPlus, LambdaMinus, AlphaPlus, AlphaMinus, 
	             CPlus, CMinus, h, er, RandGen, N, &ptPrice, &priceError);
alltime = ((double)clock())/(double)CLOCKS_PER_SEC - start_time;
	
  printf("Price = %f\n", ptPrice);
  printf("Price Error = %f\n", priceError);
   printf("Elapsed (processor) time in seconds= %f\n", alltime);
 

  

  exit (0);  
}
