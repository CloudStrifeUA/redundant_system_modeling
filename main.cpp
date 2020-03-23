#include <iostream>
#include "node.h"
#include <random>
#include <chrono>
#include <iomanip>


int main() {
    //region Генератори
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::knuth_b g(seed);
    const double lambda = 0.1, nu = 4, eta = 3;
    auto *r = new std::exponential_distribution<double>(lambda);    //
    auto *r1 = new std::exponential_distribution<double>(nu);      //
    auto *r2 = new std::exponential_distribution<double>(eta);      //
    //endregion
    const long double max_time = 1e+05;
    const int k1 = 2,k2 = 2,n = 1,m = 1;
    std::vector<double> vt,indicator;
    double time = 0,theta,ttnr = 0, rho_1 = lambda/nu, rho_2 = lambda/eta;
    Node n_1(n), n_2(m);
    vt.push_back(time);
    (n_1.getMu()+n_1.getM() == k1 && n_2.getMu()+n_2.getM() == k2)?
    indicator.push_back(1):indicator.push_back(0);
    while (time <= max_time){
        double min_1 = n_1.getZ(), min_2 = n_2.getZ();
        theta = std::min(std::min(ttnr,min_1),min_2);

        if(theta == ttnr){
            ttnr = (*r)(g);
            n_1.subtractValue(theta);
            if(n_1.getMu() < n){
                n_1.setMu(n_1.getMu()+1);
                n_1.addValue((*r1)(g));
            }
            else if(n_1.getMu() == n){
                n_1.setM(n_1.getM()+1);
            }
            n_2.subtractValue(theta);
        }

        if(theta == min_1){
            ttnr -= theta;
            n_1.subtractValue(theta);
            if(n_1.getM() > 0){
                n_1.addValue((*r1)(g));
                n_1.setM(n_1.getM()-1);
            }
            else if(n_1.getM() == 0){
                n_1.setMu(n_1.getMu()-1);
            }
            n_2.subtractValue(theta);
            if(n_2.getMu() < m){
                n_2.setMu(n_2.getMu()+1);
                n_2.addValue((*r2)(g));
            }
            else if(n_2.getMu() == m){
                n_2.setM(n_2.getM()+1);
            }
        }

        if(theta == min_2){
            ttnr -= theta;
            n_1.subtractValue(theta);
            n_2.subtractValue(theta);
            if(n_2.getM() > 0){
                n_2.addValue((*r2)(g));
                n_2.setM(n_2.getM()-1);
            }
            else if(n_2.getM() == 0){
                n_2.setMu(n_2.getMu()-1);
            }
        }
        time += theta;
        vt.push_back(time);
        (n_1.getMu()+n_1.getM() == k1 && n_2.getMu()+n_2.getM() == k2)?
        indicator.push_back(1):indicator.push_back(0);
        std::cout << std::setprecision(20) << time << std :: endl;

    }
    delete r;
    delete r1;
    delete r2;
    double res = 0;
    for(int i = 0; i < indicator.size()-1; ++i){
        res += (vt[i+1]-vt[i])*indicator[i];
    }
    double eq = res/time, eq1 = (1-rho_1)*pow(rho_1,k1)*(1-rho_2)*pow(rho_2,k2);
    std::cout << "P(k1,k2) = " << eq << std::endl;
    std::cout << eq1 << std::endl;
    std::cout << fabs(eq - eq1) << std::endl;
    return 0;
}