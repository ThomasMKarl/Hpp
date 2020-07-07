#include "include/algorithm.h"

template <typename T, class M>
void HB::productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, bool cuda)
{
    energies.resize(0);
    
    if(cuda)
    {/*
        //curand
        //copy
        size_t portion = V/threads;
        size_t over = V%threads;
        if(over != 0) ++portion;
        for(size_t i = 0; i < steps; ++i)
        {
	  model.cuda_MCSweep<<<0,0>>>(map,grid,V,beta,B,J,mE,threads,portion,over,state);
          //copy energy
          //std::cout << i << " " << mE/V << std::endl;
          energies.push_back(mE/V);
	  }*/
        //copy grid
    }
    else
    {
        gsl_rng *rng;
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        long seed = time(NULL);
        gsl_rng_set(rng,seed);   

        for(size_t i = 0; i < steps; ++i)
        {
            model.MCSweep(rng);
            //std::cout << i << " " << mE/V << std::endl;
            energies.push_back(model.calcEnergy()/T(V));
        }
    }    
}

template <typename T, class M>
void HB::productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, bool cuda)
{
    energies.resize(0);
    magn.resize(0);
    
    for(size_t i = 0; i < steps; ++i)
    {
        model.MCSweep();
        //std::cout << i << " " << V << std::endl;
        energies.push_back(model.calcEnergy()/T(V));
	//magn.push_back(mean(mgrid_theta));
    }    
}

template <typename T, class M>
void HB::productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, T T_begin, T T_end, T T_step, bool cuda)
{
    energies.resize(0);
    magn.resize(0);

    T beta;

    for(T t = T_begin; t <= T_end; t += T_step)
    {
        beta = 1.0/t;
        for(size_t i = 0; i < steps; ++i)
        {
            model.MCSweep();
            //std::cout << i << " " << mE/V << std::endl;
	    if(t + T_step > T_end)
	    {
                energies.push_back(model.calcEnergy()/T(V));
		//magn.push_back(mean(mgrid_theta));
	    }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
T mean(std::vector<T> &o)
{
    size_t n = o.size();
    
    T mean = 0.0;
    #pragma omp parallel for reduction(+:mean)
    for(size_t i = 0; i < n; ++i) mean += o[i];
    
    return mean/T(n);
}

template <typename T>
T covFunc(size_t t, std::vector<T> &o) //autocovariance-function
{
    size_t n = o.size() - t;
    
    T a = 0.0;
    T b = 0.0;
    T c = 0.0;

    #pragma omp parallel for reduction(+:a,b,c)
    for(size_t i = 0; i < n; ++i)
    {
        a += o[i]*o[i+t];
        b += o[i];
        c += o[i+t];
    } 
    
    a /= T(n);
    b /= T(n);
    c /= T(n);
 
    return a-(b*c);
}

template <typename T>
T intAuto(std::vector<T> &o) //integrated autocorrelation-time
{
    size_t n = o.size();

    T sum = 0;
    T var = cov_func(0, o);
    
    T cov;
    for(size_t t = 1; t < n; t++)
    {
        cov = cov_func(t, o);
        if(cov > 0)
        {
	    sum += (1-t/double(n)) * cov/var;
        }
        else break;
    }

    return 0.5 + sum;
}

template <typename T>
T error(std::vector<T> &o) //error on expectation for given observable
{
    return sqrt(cov_func(0, o)/T(o.size()));
}

template <typename T>
T blockingError(std::vector<T> &o, size_t block_number)
{
    size_t n = o.size();
    size_t block_size = o.size()/block_number;
    std::vector<T> block(block_size);
    std::vector<T> variances(block_number);

    for(size_t i = 0; i < block_number; ++i)
    {
        if(i == block_number - 1) block_size += n%block_number;
	
        #pragma omp parallel for
        for(size_t j = 0; j < block_size; ++j)
        {
            block[j] = o[j + i*block_number];
        }
	
        variances[i] = cov_func(0, block);
    }
    
    return sqrt(cov_func(0, variances)/T(block_number));
}

template <typename T>
T bootstrapError(std::vector<T> &o, size_t sample_size, size_t sample_number, T tau)
{
    size_t n = o.size();
    std::vector<T>   sample(sample_size);
    std::vector<T> resample(sample_number);
   
    for(size_t i = 0; i < sample_number; ++i)
    {
        #pragma omp parallel
        {
	    gsl_rng *rng;
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            long seed = time(NULL);
            gsl_rng_set(rng,seed);
	    
            #pragma omp for
            for(size_t j = 0; j < sample_size; ++j)
            {
	        sample[j] = o[d2i(n/(2*tau))*gsl_rng_uniform(rng)];
            }
	}
        resample[i] = cov_func(0, sample);
    }
    
    return sqrt(cov_func(0, resample));
}

template <typename T>
T HB::errorProp(std::vector<T> &o, T tau)
{
    size_t n = o.size();
    T mean_x;
    {
        std::vector<T> squares(n);
        #pragma omp parallel for
        for(size_t i = 0; i < n; ++i) squares[i] = o[i]*o[i]; 
        mean_x = mean(squares);
    }
    T mean_y = mean(o);

    std::vector<T> f(n);
    #pragma omp parallel
    {
        T sq;
	
        #pragma omp for
        for(size_t i = 0; i < n; ++i) 
        {
            sq = o[i]*o[i];
	    f[i] = (sq - mean_x) - 2.0*o[i]*(sq - mean_y);
        }
    }

    return sqrt( 2.0*cov_func(0,f) * tau/T(f.size()) );
}

/**************************************************************************//**
* \brief Calculates the statistical error on energies using error propagation
* for a give autocorrelation time tau                                         
******************************************************************************/
template <typename T>
T HB::errorPropEnergy(std::vector<T> &energies, T beta, T tau, size_t V)
{
  return beta*beta*error_prop(energies, tau)/T(V);
}

/***************************************************************************************//**
* \brief Calculates the statistical error on magnetizations (magn) using error propagation 
* for a given autocorrelation time tau and volume V                                                     
*******************************************************************************************/
template <typename T>
T HB::errorPropMagnetization(std::vector<T> &magn, T beta, T tau, size_t V)
{
  return beta*error_prop(magn, tau)/T(V);
}

template <typename T>
size_t HB::d2i(T d) //round double to int
{
    if(d < 0) d *= -1;
    return d<0?d-.5:d+.5;
}

template <typename T>
void HB::removeCorr(std::vector<T> &o)
{
    std::vector<T> help;
    T autot = int_auto(o);
    size_t n = o.size();
    size_t therm = d2i(20*autot);
    size_t corr;

    if(therm >= n)
    {
        std::cerr << n << " values given, but thermalisation needed " << therm << " steps!" << "\n";
        return;
    }   
    //std::cout << therm << std::endl;

    n -= therm;
    help.resize(n);
    #pragma omp parallel for
    for(size_t i = 0; i < n; ++i) help[i] = o[i+therm];
        
    autot = int_auto(o);
    //std::cout << autot << std::endl;
    
    corr = d2i(2.0*autot);
    if(corr == 1) return;
    n /= corr;
    
    o.resize(n);
    #pragma omp parallel for
    for(size_t i = 0; i < n; ++i) o[i] = help[i*corr];
}

/****************************************************************************//**
* \brief Computes the specific heat as secondary quantity depending on energies  
********************************************************************************/
template <typename T>
T HB::specificHeat(std::vector<T> &energies, T beta, size_t V)
{
  return beta*beta*cov_func(0, energies)/T(V);
}

/****************************************************************************//**
* \brief Computes the magneic suszeptibility as secondary quantity depending on 
* magnetizations (magn).                                                        
********************************************************************************/
template <typename T>
T HB::magnSusz(std::vector<T> &magn, T beta, size_t V)
{
  return beta*cov_func(0, magn)/T(V);
}
