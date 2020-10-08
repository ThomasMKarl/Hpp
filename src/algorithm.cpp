#include "algorithm.h"

/*__global__
void HB::setupKernel(curandState *state, size_t seed, size_t n)
{
  size_t id = threadIdx.x + blockIdx.x*blockDim.x;
  curand_init(seed, id, 0, &state[id]);
}*/

template <typename T, class M>
void HB::productionRun(M &model, T beta, size_t V, size_t steps, std::vector<T> &energies, bool cuda)
{
    energies.resize(0);

    T energy = model.calcEnergy();

    T Volume = T(V);

    energies.push_back(energy/Volume);
    
    if(cuda)
    {
      /*size_t threads = 256;
        size_t blocks = (V+threads-1)/threads;

	model->getGrid()->upload();

	thrust::device_vector<T> energies_to_reduce;
	energies_to_reduce.resize(V);

        curandState *devStates;
        CUDA_CALL(cudaMalloc((void **)&devStates, V*sizeof(curandState)));
        setup_kernel<<<blocks,threads>>>(devStates, time(NULL), V)
        for(size_t i = 0; i < steps; ++i)
        {
	    cuda_MCSweep<<<blocks,threads>>>(model->getGrid()->getDeviceTable(),
					     model->getGrid()->getDeviceGrid(),
					     V,
					     beta,
					     model->getGrid()->getCoupling(),
					     thrust::raw_pointer_cast(energies_to_reduce.data()),
					     state);
	    energy += thrust::reduce(thrust::plus, energies_to_reduce.first(), energies_to_reduce.last(), 0.0); 	
	    
            //std::cout << i << " " << energy/V << std::endl;
            energies.push_back(energy/Volume);
	}
	
        model->getGrid()->download();

        cudaFree(devStates);*/
    }
    else
    {
        gsl_rng *rng;
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        long seed = time(NULL);
        gsl_rng_set(rng,seed);   

        for(size_t i = 0; i < steps; ++i)
        {
	    energy += HB::MCSweep(model, beta, rng);
            //std::cout << i << " " << energy/V << std::endl;
            energies.push_back(energy/Volume);
        }
    }    
}

template <typename T, class M>
void HB::productionRun(M &model, T beta, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, bool cuda)
{
    energies.resize(0);
    magn.resize(0);

    T energy = model.calcEnergy();

    T Volume = T(V);

    energies.push_back(energy/Volume);
    magn.push_back(model.calcMagn()/Volume);

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long seed = time(NULL);
    gsl_rng_set(rng,seed);       
    
    for(size_t i = 0; i < steps; ++i)
    {
	energy += HB::MCSweep(model, beta, rng);
        //std::cout << i << " " << energy/V << std::endl;
        energies.push_back(energy/Volume);
	magn.push_back(model.calcMagn()/Volume);
    }    
}

template <typename T, class M>
void HB::productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, T T_begin, T T_end, T T_step, bool cuda)
{
    energies.resize(0);
    magn.resize(0);

    T beta;

    T energy = model.calcEnergy();

    T Volume = T(V);

    energies.push_back(energy/Volume);
    magn.push_back(model.calcMagn()/Volume);

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long seed = time(NULL);
    gsl_rng_set(rng,seed);  

    for(T t = T_begin; t <= T_end; t += T_step)
    {
        beta = 1.0/t;
        for(size_t i = 0; i < steps; ++i)
        {
            energy += HB::MCSweep(model, beta, rng);
            //std::cout << i << " " << energy/V << std::endl;
	    if(t + T_step > T_end)
	    {
                energies.push_back(energy/Volume);
	        magn.push_back(model.calcMagn()/Volume);
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
	    sum += (1-t/T(n)) * cov/var;
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

template <typename T>
T HB::errorPropEnergy(std::vector<T> &energies, T beta, T tau, size_t V)
{
  return beta*beta*error_prop(energies, tau)/T(V);
}

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

template <typename T>
T HB::specificHeat(std::vector<T> &energies, T beta, size_t V)
{
  return beta*beta*cov_func(0, energies)/T(V);
}

template <typename T>
T HB::magnSusz(std::vector<T> &magn, T beta, size_t V)
{
  return beta*cov_func(0, magn)/T(V);
}

template <typename T, class M>
T HB::MCSweep(M &model, T beta, gsl_rng *rng) //metropolis algorithm
{
    T new_E = 0.0;

    dim3 N = model->getGrid()->getGridSize();  
    size_t V = N.x*N.y*N.z;
    
    #pragma omp parallel
    {	
        T dE;
	dim3 idx;
	
        #pragma omp for reduction(+:new_E)   	
        for(size_t j = 0; j < V; ++j)
	{
            idx.x = j% N.x;
	    idx.y = j/ N.x%N.y;
            idx.z = j/(N.x*N.y);
	    
    	    dE = model.calcEnergy(idx);
	  
            if(dE < 0 || exp(-beta*dE) > gsl_rng_uniform(rng))
            {
	        model.flip(idx, rng);
                new_E += dE;
            }
        }
    }
    
    return new_E;
}

/*template <typename T, class S> __global__
HB::cuda_MCSweep(HB::NB *map, S *grid, dim3 gridSize, T beta, T J, T *new_E, curandState *state) //metropolis algorithm
{ 
    T dE;

    dim3 idx;
    
    curandState localState = state[id];
    
    idx.x = id% gridSize.x;
    idx.y = id/ gridSize.x%gridSize.y;
    idx.z = id/(gridSize.x*gridSize.y);
	    
    dE = model.cuda_calcEnergy(idx);
	  
    if(dE < 0 || exp(-beta*dE) > curand_uniform(&localState))
    {
	model.cuda_flip(idx, rng);
        new_E[id] = dE;
    }
    else new_E[id] = 0.0;
    
    state[id] = localState;
    }*/
