#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>

#define LAMBDA 1.0
#define NU 10
#define RHO LAMBDA/NU
#define MIN_REPLACEMENT 3
#define MAX_REPLACEMENT 100
#define T 10
#define N 1e+04
struct Logger
{
    std::ostream& stream;
    Logger(std::ostream& os):stream(os){}
};

long double factor(int n)
{
    if(n == 1 || n == 0)
        return 1;
    return n*factor(n - 1);
}

long double poisson_prob(long double lambda, int k)
{
    return powl(lambda, k)* expl(-lambda)/factor(k);
}

long double poisson_sum(long double lambda, int u)
{
    long double res = 0;
    for(int i = u; i < u + 10; i++)
    {
        res += poisson_prob(lambda, i);
    }
    return res;
}

long double q(int n)
{
    return (1 - RHO)*powl(RHO, n - 1)/(1 - powl(RHO, n));
}

long double potf(const double lambda, const double nu, const double t, const uint replacement, const uint n, std::vector<double> &probs)
{
//region Генератори
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::knuth_b g(seed);
    auto *r1 = new std::exponential_distribution<double>(nu);      //
//endregion
    long double res = 0;
    for(int i = 0; i < n; i++)
    {
        long double time = 0.0, p = 1.0;
        int q = 0;
        while(time < t)
        {
            long double dt = t - time;       //t - t_k
            long double eta = (*r1)(g);
            if(q == 0)
            {
                p *= (1 - expl(-lambda * dt));
                long double tau = -logl(1 - (std::rand()/(long double)RAND_MAX) * (1 - expl(-lambda * dt)))/lambda;
                if(tau + eta >= dt)
                {
                    p *= poisson_sum(lambda*(dt - tau), replacement);
                    break;
                }
                else
                {
                    auto *m = new std::poisson_distribution<int>(lambda*eta);
                    int m_k = (*m)(g);
                    delete m;
                    if(m_k >= replacement)
                        break;
                    else q += m_k;
                    time += eta + tau;
                }
            }
            else
            {
                if(eta >= dt)
                {
                    p *= poisson_sum(lambda*dt, replacement + 1 - q);
                    break;
                }
                else
                {
                    auto *m = new std::poisson_distribution<int>(lambda*eta);
                    int m_k = (*m)(g);
                    delete m;
                    if(m_k >= replacement + 1 - q)
                        break;
                    else q+= m_k - 1;
                    time += eta;
                }
                
            }
        }
        probs.push_back(p);
        res += p;
    }
    return res/n;
}

int main() 
{
    std::ofstream f("log.txt", std::ios_base::out|std::ios_base::app);
    Logger logger(f);
    logger.stream << "lambda = " << LAMBDA << ", nu = " << NU << ", t = " << T << ", N = " << N << std::endl;
    logger.stream << std::setw(5) << "n" << " |" << std::setw(15) << "delta(t,n) " << std::setw(6) << "|" << " delta(t,n)/(rho^n) " << "|"
    <<  "delta(t,n)/(n*rho^n)" << "|" << " Iмовірність відмови" << "|" << std::setw(11) << "s^2" << std::setw(10) 
    << "|" << std::endl;
    long double a = 1/(LAMBDA*(1 - RHO));
    std::vector<double> probs;
    for(size_t u = MIN_REPLACEMENT; u <= MAX_REPLACEMENT; u++)
    {
        long double p = potf(LAMBDA, NU, T, u, N, probs);
        long double div = 0.0;
        for(auto i:probs)
        {
            div += (i - p)*(i - p);
        }
        div /= probs.size() - 1;
        long double flQ = q(u);
        long double delta = fabsl(p - flQ > 0.01 ? (1 - expl(-T*q(u)/a)):(T*flQ/a - powl(T*flQ/a, 2)/2 + powl(T*flQ/a, 3)/6 - powl(T*flQ/a, 4)/24));
        logger.stream << std::setprecision(12) << std::scientific;
        logger.stream << std::setw(5) << u << " | " << delta << " | " << delta/(powl(RHO, u))  << " | " << delta/(powl(RHO, u)*u)  << " | " << p  << " | " << div  << " | " << std::endl;
    }
    return 0;
}