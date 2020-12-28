#include "simulation.h"

/*__global__
void HB::setupKernel(curandState *state, size_t seed, size_t n)
{
  size_t id = threadIdx.x + blockIdx.x*blockDim.x;
  curand_init(seed, id, 0, &state[id]);
}*/

template <class M>
float HB::MetropolisSimulation::metropolisSweep(const M &model)
{
    float new_E = 0.0;

    const dim3 N = mGrid->getGridSize();  
    const size_t V = N.x*N.y*N.z;
    	
    float dE;
    dim3 idx;
    
    #pragma omp parallel for reduction(+:new_E) private(dE,idx)   	
    for(size_t j = 0; j < V; ++j)
    {
        std::random_device device;
        std::mt19937 generator(device());
        std::uniform_real_distribution<float> distribution(0,1);
	
        idx.x = j% N.x;
	idx.y = j/ N.x%N.y;
        idx.z = j/(N.x*N.y);
	    
    	dE = model.calcEnergy(mGrid, idx);
	  
        if(dE < 0 || exp(-mBeta*dE) > distribution(generator))
        {
	    model.flip(idx);
            new_E += dE;
        }
    }
    
    return new_E;
}

template <class M>
float HB::MetropolisSimulationQt::metropolisSweep(const M &model)
{
    float new_E = 0.0;

    const dim3 N = mGrid->getGridSize();  
    const size_t V = N.x*N.y*N.z;
    	
    float dE;
    dim3 idx;
    
    #pragma omp parallel for reduction(+:new_E) private(dE,idx)   	
    for(size_t j = 0; j < V; ++j)
    {
        std::random_device device;
        std::mt19937 generator(device());
        std::uniform_real_distribution<float> distribution(0,1);
	
        idx.x = j% N.x;
	idx.y = j/ N.x%N.y;
        idx.z = j/(N.x*N.y);
	    
    	dE = model.calcEnergy(mGrid, idx);
	  
        if(dE < 0 || exp(-mBeta*dE) > distribution(generator))
        {
	    model.flip(idx);
            new_E += dE;
        }
    }
    
    return new_E;
}

/*template <typename T>
__global__
void HB::cuda_MCSweep(T *n, T *gridData, size_t V, T beta, T J, T *energies, curandState *devStates) //metropolis algorithm
{
    size_t id = threadIdx.x + blockIdx.x*blockDim.x;
    if(id < V)
    {
        T new_E = 0.0;	
        T dE;
        dim3 idx;
       	
        idx.x = j% N.x;
	idx.y = j/ N.x%N.y;
        idx.z = j/(N.x*N.y);
	    
    	dE = HB::Device::calcEnergy(gridData, n, idx);
	  
        if(dE < 0 || exp(-beta*dE) > gsl_rng_uniform(rng))
        {
            HB::Device::flip(gridData, idx, rng);
            new_E += dE;
        }

	energies[idx] = new_E;
    }
}*/

/*template <typename T, class M>
void HB::metropolisRun(M &model, HB::DeviceGrid &grid, T beta, size_t steps, std::vector<T> &energies)
{
    energies.resize(0);

    T energy = model.calcEnergy();
    
    dim3 N = grid.getGridSize();
    size_t Volume = N.x*N.y*N.z;

    energies.push_back(energy/Volume);
    
    size_t threads = 256;
    size_t blocks = (V+threads-1)/threads;

    thrust::device_vector<T> energies_to_reduce;
    energies_to_reduce.resize(V);

    curandState *devStates;
    CUDA_CALL(cudaMalloc((void **)&devStates, V*sizeof(curandState)));
    setup_kernel<<<blocks,threads>>>(devStates, time(NULL), V);
    CUDA_CALL(cudaDeviceSynchronize());
	      
    for(size_t i = 0; i < steps; ++i)
    {
        cuda_MCSweep<<<blocks,threads>>>(grid->getTable().data(),
					 grid.data(),
					 V,
					 beta,
					 grid.getCoupling(),
					 thrust::raw_pointer_cast(
					     energies_to_reduce.data()),
					 devStates);
	
	CUDA_CALL(cudaDeviceSynchronize());
	
	energy += thrust::reduce(thrust::plus,
				 energies_to_reduce.first(),
				 energies_to_reduce.last(),
				 0.0); 	
	    
        //std::cout << i << " " << energy/V << std::endl;
        energies.push_back(energy/Volume);
    }

    CUDA_CALL(cudaFree(devStates));
}*/

template <class M>
void HB::MetropolisSimulation::simulate(const M &model)
{
    mEnergies.resize(0);
    mMagnetization.resize(0);

    float energy = model.calcEnergy(mGrid);

    const dim3 N = mGrid->getGridSize();
    const double Volume = double(N.x*N.y*N.z);

    mEnergies.push_back(energy/Volume);
    mMagnetization.push_back(model.calcMagn(mGrid)/Volume);
    
    for(size_t i = 0; i < mSteps; ++i)
    {
      energy += metropolisSweep(model);
      //std::cout << i << " " << energy/V << std::endl;
      mEnergies.push_back(energy/Volume);
      mMagnetization.push_back(model.calcMagn(mGrid)/Volume);
    }    
}

template <class M>
void HB::MetropolisSimulationQt::simulate(const M &model)
{
    QPainter qp(this);
    unsigned PenSize = 10;
    unsigned offset = PenSize ? PenSize : 1;
    QPen WhitePen(Qt::white, PenSize, Qt::SolidLine),
         BlackPen(Qt::black, PenSize, Qt::SolidLine);

    float energy = model.calcEnergy(mGrid);

    metropolisSweep(model);
    
    const dim3 N = mGrid->getGridSize();
    dim3 index;
    for(unsigned int i = 0; i < N.x; i++)
    {
        for(unsigned int j = 0; i < N.y; j++)
        {
            index.x = i; index.y = j; index.z = 0;
            if(mGrid->getSpin(index) == 1) qp.setPen(WhitePen);
            else                           qp.setPen(BlackPen);

	    qp.drawPoint(offset * i, offset * j);
	    //draw white point if spin 1, black if -1
        }
    }   
}

template <class M>
void HB::SimulatedAnnealing::simulate(const M &model)
{
    mEnergies.resize(0);
    mMagnetization.resize(0);

    float energy = model.calcEnergy(mGrid);

    const dim3 N = mGrid->getGridSize();
    const double Volume = double(N.x*N.y*N.z);

    mEnergies.push_back(energy/Volume);
    mMagnetization.push_back(model.calcMagn(mGrid)/Volume);
    
    float beta;
    for(float t = mTemperatureSteps.x; t <= mTemperatureSteps.y; t += mTemperatureSteps.z)
    {
        beta = 1.0f/t;
        for(size_t i = 0; i < mSteps; ++i)
        {
	    energy += metropolisSweep(model);
            //std::cout << i << " " << energy/V << std::endl;
	    if(t + mTemperatureSteps.z > mTemperatureSteps.y)
	    {
                mEnergies.push_back(energy/Volume);
	        mMagnetization.push_back(model.calcMagn(mGrid)/Volume);
	    }
        }
    }
}

void simulate(const HB::Models &models)
{
   for(const auto &model : models)
   {
     model->simulate();
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////


template <typename T>
T mean(const std::vector<T> &o)
{
    const size_t n = o.size();
    
    T mean = 0.0;
    #pragma omp parallel for reduction(+:mean)
    for(size_t i = 0; i < n; ++i) mean += o[i];
    
    return mean/T(n);
}

template <typename T>
T covFunc(const size_t t, const std::vector<T> &o) //autocovariance-function
{
    const size_t n = o.size() - t;
    
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
T intAuto(const std::vector<T> &o) //integrated autocorrelation-time
{
    const size_t n = o.size();

    T sum = 0;
    const T var = cov_func(0, o);
    
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
T error(const std::vector<T> &o) //error on expectation for given observable
{
    return sqrt(cov_func(0, o)/T(o.size()));
}

template <typename T>
T blockingError(const std::vector<T> &o, const size_t block_number)
{
    const size_t n = o.size();
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
T bootstrapError(const std::vector<T> &o, const size_t sample_size, const size_t sample_number, const T tau)
{
    const size_t n = o.size();
    std::vector<T>   sample(sample_size);
    std::vector<T> resample(sample_number);
   
    for(size_t i = 0; i < sample_number; ++i)
    {
        #pragma omp parallel for
        for(size_t j = 0; j < sample_size; ++j)
        {
	    std::random_device device;
            std::mt19937 generator(device());
            std::uniform_real_distribution<float> distribution(0,1);
	    
	    sample[j] = o[d2i(n/(2*tau))*distribution(generator)];
        }
	    
        resample[i] = cov_func(0, sample);
    }
    
    return sqrt(cov_func(0, resample));
}

template <typename T>
T HB::errorProp(const std::vector<T> &o, const T tau)
{
    const size_t n = o.size();
    T mean_x;
    {
        std::vector<T> squares(n);
        #pragma omp parallel for
        for(size_t i = 0; i < n; ++i) squares[i] = o[i]*o[i]; 
        mean_x = mean(squares);
    }
    const T mean_y = mean(o);

    std::vector<T> f(n);
    T square;	
    #pragma omp parallel for private(square)
    for(size_t i = 0; i < n; ++i) 
    {
        square = o[i]*o[i];
	f[i] = (square - mean_x) - 2.0*o[i]*(square - mean_y);
    }

    return sqrt( 2.0*cov_func(0,f) * tau/T(f.size()) );
}

template <typename T>
T HB::errorPropEnergy(const std::vector<T> &energies, const T beta, const T tau, const size_t V)
{
  return beta*beta*error_prop(energies, tau)/T(V);
}

template <typename T>
T HB::errorPropMagnetization(const std::vector<T> &magn, const T beta, const T tau, const size_t V)
{
  return beta*error_prop(magn, tau)/T(V);
}

template <typename T>
size_t HB::d2i(const T d) //round double to int
{
    if(d < 0) d *= -1;
    return d<0?d-.5:d+.5;
}

template <typename T>
void HB::removeCorr(const std::vector<T> &o)
{
    std::vector<T> help;
    T autot = int_auto(o);
    size_t n = o.size();
    const size_t therm = d2i(20*autot);
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
T HB::specificHeat(const std::vector<T> &energies, const T beta, const size_t V)
{
  return beta*beta*cov_func(0, energies)/T(V);
}

template <typename T>
T HB::magnSusz(const std::vector<T> &magn, const T beta, const size_t V)
{
  return beta*cov_func(0, magn)/T(V);
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
