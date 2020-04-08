#include <iostream>
#include "node.h"
#include <random>
#include <chrono>
#include <iomanip>


int main() {
    //region Генератори
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::knuth_b g(seed);
    const double lambda = 1, nu = 1, eta = 3;
    auto *r = new std::exponential_distribution<double>(lambda);    //
    auto *r1 = new std::exponential_distribution<double>(nu);      //
    // auto *r2 = new std::exponential_distribution<double>(eta);      //
    //endregion
    const unsigned replacement = 1;
    const double t = 1;
    const int n = 1, m = 1, N = 1e+05;
    int res = 0;
    for(int i = 0; i < N; i++)
    {
        Node n1(n), n2(m);
        double time = 0.0;
        n1.setM(replacement);
        for(int i = 0; i < n; i++)
            n1.addValue((*r)(g));
        n1.setMu(n);
        while(time < t)
        {
            double t1 = n1.getZ();
            double t2 = n2.getZ();
            double theta = std::min(t1, t2);
            std::cout << n1.getZ() << std::endl;
            std::cout << n2.getZ() << std::endl;
            std::cout << theta << std::endl;
            if(n1.getMu() != 0)
                n1.subtractValue(theta);
            if(n2.getMu() != 0)
                n2.subtractValue(theta);
            if(theta == t1)
            {
                if(n1.getM() > 0)
                {
                    n1.setM(n1.getM() - 1);
                    n1.addValue((*r)(g));
                }
                else n1.setMu(n1.getMu() - 1);
                if(n2.getMu() < m)
                {
                    n2.addValue((*r1)(g));
                    n2.setMu(n2.getMu() + 1);
                }
                else n2.setM(n2.getM() + 1);
            }  
            else if(theta == t2)
            {
                if(n2.getM() > 0)
                {
                    n2.addValue((*r1)(g));
                    n2.setM(n2.getM() - 1);
                }
                else n2.setMu(n2.getMu() - 1);
                if(n1.getMu() < n)
                {
                    n1.addValue((*r)(g));
                    n1.setMu(n1.getMu() + 1);
                }
                else n1.setM(n1.getM() + 1);
            }
            time += theta;
            if(n1.getM() == 0 && n1.getMu() == 0 && time < t)
            {
                res += 1; 
                break;
            }
        }
    }
    std::cout << 1 - (double)res/N << std::endl;
    return 0;
}