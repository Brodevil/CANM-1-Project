#include "Rootfinder.h"
#include <functional>
#include <iostream>
#include <cmath>
using namespace std;


Rootfinder::Rootfinder(double tolerance, bool verbose){
    this->tol=tolerance;
    this->verbose=verbose;
}

double Rootfinder::Bisection_method(std::function<double(double)> f, double higher, double lower,int num_iter){
    if( f(lower)*f(higher)>0){
        cout<<"ERROR: the range specificed has even number of roots within it.\nTry other methods instead, or enter a valid range.";
        return 0.0;
    }
    for(int i=0;i<num_iter;i++){
        if(abs(f(lower))<(this->tol)){
            if(this->verbose)cout<<"Root found in "<<i<<"iterations, with"<<f(lower)<<"residue";
            return f(lower);
        }
        if(abs(f(higher))<(this->tol)){
            if(this->verbose)cout<<"Root found in "<<i<<"iterations, with"<<f(higher)<<"residue";
            return f(higher);
        }
        
        double mid=(lower+higher)/2;
        if(f(mid)*f(higher)<0){
            lower=mid;
        }
        else{
            higher=mid;
        }
    }

    cout<<"ERROR:Bisection method failed to converge within "<<num_iter<< "iterations. Kindly try other methods instead.\n";
    return 0.0;

}

double Rootfinder::FDMderivative(std::function<double(double)> f,double x, double h){
    return ((f(x+h)-f(x-h))/(2*h)); // uses Central differences, better than one sided difference
    
}

double Rootfinder::Newton_Raphson(std::function<double(double)> f,double x0,int num_iter){
    double xprev=x0;
    double x_next=xprev;

    for(int i=0;i<num_iter;i++){
       xprev=x_next;
       if(abs(f(xprev))<(this->tol))return xprev;
       x_next=xprev-f(xprev)/(FDMderivative(f,xprev));
    }

    cout<<"ERROR:Newton-Raphson failed to converge within "<<num_iter<<"iterations. There is likely a numerical instability arriving within the function. Kindly try other methods\n";
    return 0.0;
}


double Rootfinder::Secant_Method(std::function<double(double)> f,double x0,int num_iter){
    double x_n2=x0;
    double x_n1=x0*(1+1e-6);
    double x_n=x_n1-(x_n1-x_n2)/(f(x_n1)-f(x_n2));
    for(int i=0;i<num_iter;i++){
       if(abs(f(x_n))<(this->tol))return x_n;
       x_n2=x_n1;
       x_n1=x_n;
       x_n=x_n1-(x_n1-x_n2)/(f(x_n1)-f(x_n2));
    }

    cout<<"ERROR:Newton-Raphson failed to converge within "<<num_iter<<"iterations. There is likely a numerical instability arriving within the function. Kindly try other methods\n";
    return 0.0;
}



