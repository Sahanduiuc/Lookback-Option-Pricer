//
//  main.cpp
//  volindex
//
//  Created by Александр Гречко on 21.07.15.
//  Copyright (c) 2015 Alexander Grechko. All rights reserved.
//

#include <vector>
#include <deque>
#include <algorithm>
#include <pnl/pnl_finance.h>
#include <pnl/pnl_cdf.h>


typedef struct option {
    double strike;
    double call_price;
    double put_price;
} option;

double get_d2(double k,double t,double v)
{
    double a=v*sqrt(t);
    double d2=-k/a-a/2;
    return d2;
}


//compute distance between 2 points (x1,y1) and (x2,y2)
double get_distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//compute volatility index, fn - file name with option prices
//option strikes in file must be sorted
//s - asset price
//T - time to maturity
double vol_index(const char* fn, double s, double T)
{
    const int max_options=30;
    option options[max_options];
    
    //Reading option data from file fn
    FILE* f=fopen(fn, "r");
    if (f == NULL)
    {
        perror("Could not open file");
        return 0;
    }
    int n=0;
    int iatm=-1; // index of atm put
    while (!feof(f))
    {
        fscanf(f, "%lf %lf %lf\n",&options[n].strike,&options[n].call_price,&options[n].put_price);
        if (n>=1)
            if ((options[n-1].strike<=s) && (s<options[n].strike))
                iatm=n-1;
        n++;
    }
    printf("N = %d\n",n);
    printf("ATM stike = %g\n",options[iatm].strike);
    fclose(f);
    
    std::deque<double> x,y;
    int e;
    for (int i=iatm;i>=0;i--)
    {
        if (options[i].put_price==0)
            continue;
        double v=pnl_bs_implicit_vol(0, options[i].put_price, s, options[i].strike, T, 0, 0, &e);
        double k=log(options[i].strike/s);
        double d2=get_d2(k,T,v);
        if (!x.empty())
            if (d2<x.front())
                break;
        x.push_front(d2);
        y.push_front(v*v);
    }
    
    for (int i=iatm+1;i<n;i++)
    {
        if (options[i].call_price==0)
            continue;
        double v=pnl_bs_implicit_vol(1, options[i].call_price, s, options[i].strike, T, 0, 0, &e);
        double k=log(options[i].strike/s);
        double d2=get_d2(k,T,v);
        if (!x.empty())
            if (d2>x.back())
                break;
        x.push_back(d2);
        y.push_back(v*v);
    }
    
    std::reverse(x.begin(), x.end());
    std::reverse(y.begin(), y.end());

    int m=(int)x.size();
    if (m<=1)
    {
        printf("Points must be atleast 2\n");
        return 0;
    }
    
    std::vector<double> dy(m); // derivative of y
    dy[0]=0;
    dy[m-1]=0;
    if (m>2)
    {
        double l1=get_distance(x[0], y[0], x[1], y[1]);
        for (int j=1;j<m-1;j++)
        {
            double l2=get_distance(x[j], y[j], x[j+1], y[j+1]);
            dy[j]=-((x[j+1]-x[j])/l2-(x[j]-x[j-1])/l1)/((y[j+1]-y[j])/l2-(y[j]-y[j-1])/l1);
            l1=l2;
        }
    }
    
    //create cubic polynomial approximation
    const std::deque<double>& a=y;
    const std::vector<double>& b=dy;
    std::vector<double> c(m-1);
    std::vector<double> d(m-1);
    for (int j=0;j<m-1;j++)
    {
        double dxj=x[j+1]-x[j];
        double dyj=y[j+1]-y[j];
        c[j]=(3*dyj-dxj*dy[j+1]-2*dxj*dy[j])/(dxj*dxj);
        d[j]=(dyj-dy[j]*dxj-c[j]*dxj*dxj)/(dxj*dxj*dxj);
    }
    
    //integrate approximation function
    double cdf_x1=pnl_cdfnor(x[0]);
    double sum=y[0]*cdf_x1;
    double pdf_x1=pnl_normal_density(x[0]);
    for (int j=0;j<m-1;j++)
    {
        double cdf_x2=pnl_cdfnor(x[j+1]);
        double pdf_x2=pnl_normal_density(x[j+1]);
        double ka=cdf_x2-cdf_x1;
        double kb=-(pdf_x2-pdf_x1)-x[j]*(cdf_x2-cdf_x1);
        double kc=-(x[j+1]*pdf_x2-x[j]*pdf_x1)+2*x[j]*(pdf_x2-pdf_x1)+(1+x[j]*x[j])*(cdf_x2-cdf_x1);
        double kd=(1-x[j+1]*x[j+1])*pdf_x2-(1-x[j]*x[j])*pdf_x1+3*x[j]*(x[j+1]*pdf_x2-x[j]*pdf_x1)-3*(1+x[j]*x[j])*(pdf_x2-pdf_x1)-x[j]*(3+x[j]*x[j])*(cdf_x2-cdf_x1);
        sum+=a[j]*ka+b[j]*kb+c[j]*kc+d[j]*kd;
        cdf_x1=cdf_x2;
        pdf_x1=pdf_x2;
    }
    sum+=y[m-1]*(1-cdf_x1);
    
    double vi=sqrt(sum)*100.0;
    
    return vi;
}

int main(int argc, const char * argv[]) {
    double vol=vol_index("opts-17-08-2015.txt", 88300.0, 0.07);
    printf("Volatility index: %g\n",vol);
    return 0;
}
