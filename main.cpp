#include <iostream>
#include "node.h"
#include <random>
#include <chrono>
#include <iomanip>

double factor(int n)
{
    if(n == 1 || n == 0)
        return 1;
    return n*factor(n - 1);
}

double poisson_prob(double lambda, int k)
{
    return std::pow(lambda, k)* std::exp(-lambda)/factor(k);
}

double poisson_sum(double lambda, int u)
{
    double res = 0;
    for(int i = 0; i < u; i++)
    {
        res += poisson_prob(lambda, i);
    }
    return res;
}


int main() {
    //region Генератори
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::knuth_b g(seed);
    const double lambda = 1, nu = 1;
    auto *r = new std::exponential_distribution<double>(lambda);    //
    auto *r1 = new std::exponential_distribution<double>(nu);      //
    // auto *r2 = new std::exponential_distribution<double>(eta);      //
    //endregion
    const unsigned replacement = 1;
    const double t = 1;
    const int n = 1, m = 1, N = 1e+03;
    double res = 0;
    for(int i = 0; i < N; i++)
    {
        double time = 0.0, p = 1.0;
        int q = 0;
        while(time < t)
        {
            double dt = t - time;       //t - t_k
            double eta = (*r1)(g);
            std::cout << "dt " << dt << std::endl;
            std::cout << "eta " << eta << std::endl;
            if(q == 0)
            {
                p *= (1 - std::exp(-lambda * dt));
                double tau = -std::log(1 - (std::rand()/(double)RAND_MAX) * (1 - std::exp(-lambda * dt)))/lambda;
                std::cout << "tau " << tau << std::endl;
                if(tau + eta >= dt)
                {
                    p *= 1 - poisson_sum(lambda*(dt - tau), replacement);
                    std::cout << "p " << p << std::endl;
                    break;
                }
                else
                {
                    auto *m = new std::poisson_distribution<int>(lambda*eta);
                    int m_k = (*m)(g);
                    delete m;
                    std::cout << "m_k " << m_k << std::endl;
                    if(m_k >= replacement)
                        break;
                    else q+= m_k;
                    time += eta + tau;
                }
            }
            else
            {
                if(eta >= dt)
                {
                    p *= 1 - poisson_sum(lambda*dt, replacement + 1 - q);
                    std::cout << "p " << p << std::endl;
                    break;
                }
                else
                {
                    auto *m = new std::poisson_distribution<int>(lambda*eta);
                    int m_k = (*m)(g);
                    delete m;
                    std::cout << "m_k " << m_k << std::endl;
                    if(m_k >= replacement + 1 - q)
                        break;
                    else q+= m_k - 1;
                    time += eta;
                }
                
            }
        }
        res += p;
    }
    std::cout << 1 - res/N << std::endl;
    return 0;
}