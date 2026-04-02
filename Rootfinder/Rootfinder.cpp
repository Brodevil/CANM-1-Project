#include "Rootfinder.h"
#include <functional>
#include <iostream>
#include <cmath>
#include <limits>
using namespace std;


Rootfinder::Rootfinder(double tolerance, bool verbose){
    this->tol=tolerance;
    this->verbose=verbose;
}

double Rootfinder::Bisection_method(function<double(double)> f, double higher, double lower,int num_iter){
    if( f(lower)*f(higher)>0){
        cout<<"ERROR: the range specificed has even number of roots within it.\nTry other methods instead, or enter a valid range.";
        return numeric_limits<double>::quiet_NaN();
    }
    double f_low=f(lower);
    double f_high=f(higher);
    for(int i=0;i<num_iter;i++){
        double mid=(lower+higher)/2;
        double f_mid=f(mid);
        if(abs(f_mid)<tol){
            return mid;
        }
        if(f_mid*f_high<0){
            lower=mid;
            f_low=f_mid;
        }
        else{
            higher=mid;
            f_high=f_mid;
        }
    }

    cout<<"ERROR:Bisection method failed to converge within "<<num_iter<< "iterations. Kindly try other methods instead.\n";
    return numeric_limits<double>::quiet_NaN();

}

double Rootfinder::FDMderivative(function<double(double)> f,double x, double h){
    return ((f(x+h)-f(x-h))/(2*h)); // uses Central differences, better than one sided difference
    
}

double Rootfinder::Newton_Raphson(function<double(double)> f,double x0,int num_iter){
    double xprev=x0;
    double x_next=xprev;

    for(int i=0;i<num_iter;i++){
       xprev=x_next;
       if(abs(f(xprev))<(this->tol))return xprev;
       x_next=xprev-f(xprev)/(FDMderivative(f,xprev));
    }

    cout<<"ERROR:Newton-Raphson failed to converge within "<<num_iter<<"iterations. There is likely a numerical instability arriving within the function. Kindly try other methods\n";
    return numeric_limits<double>::quiet_NaN();
}


double Rootfinder::Secant_Method(function<double(double)> f,double x0,int num_iter){
    double x_n2=x0;
    double x_n1=x0*(1+1e-6)+1e-6;
    double x_n=x_n1-f(x_n1)*(x_n1-x_n2)/(f(x_n1)-f(x_n2));
    for(int i=0;i<num_iter;i++){
       if(abs(f(x_n))<(this->tol))return x_n;
       x_n2=x_n1;
       x_n1=x_n;
       x_n=x_n1-f(x_n1)*(x_n1-x_n2)/(f(x_n1)-f(x_n2));
    }

    cout<<"ERROR:Secant method failed to converge within "<<num_iter<<"iterations. There is likely a numerical instability arriving within the function. Kindly try other methods\n";
    return numeric_limits<double>::quiet_NaN();
}



