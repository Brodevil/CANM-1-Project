#include "Rootfinder.h"
#include <functional>
#include <iostream>
#include <cmath>
using namespace std;


Rootfinder::Rootfinder(double tolerance, bool verbose){
    this->tol=tolerance;
    this->verbose=verbose;
}

double Rootfinder::Bisection_method(std::function<double(double)> f, double higher, double lower){
    if( f(lower)*f(higher)>0){
        cout<<"ERROR: the range specificed has even number of roots within it.\nTry other methods instead, or enter a valid range.";
        return -1.0;
    }
    for(int i=0;i<2e3;i++){
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

    cout<<"ERROR:Bisection method failed to converge within 2000 iterations. Kindly try other methods instead.\n";
    return -1.0;

}

double Rootfinder::FDMderivative(std::function<double(double)> f,double x, double h){
    return ((f(x+h)-f(x-h))/(2*h)); // uses Central differences, better than one sided difference
    
}

double Rootfinder::Newton_Raphson(std::function<double(double)> f,double x0){
    double xprev=x0;
    double x_next=xprev;

    for(int i=0;i<200;i++){
       xprev=x_next;
       if(abs(f(xprev))<(this->tol))return xprev;
       x_next=xprev-f(xprev)/(FDMderivative(f,xprev));
    }

    cout<<"ERROR:Newton-Raphson failed to converge within 200 iterations. There is likely a numerical instability arriving within the function. Kindly try other methods\n";
    
}


