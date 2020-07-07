#include "include/model/ising/ising.h"

double HB::Ising::calcEnergy() //store initial value for energy
{
    int spin_sum = 0;
    int neighbour_sum = 0;
    
    unsigned int n = m_Nx*m_Ny*m_Nz;
    
    #pragma omp parallel for reduction(+:spin_sum,neighbour_sum)
    for(unsigned int i = 0; i < n; ++i)
    {
        spin_sum += mgrid[i];
        neighbour_sum += mgrid[i]*mgrid[mmap[1][i]];
	if(dim > 1)
	{
            neighbour_sum += mgrid[i]*mgrid[mmap[2][i]];
	    if(dim > 2)
	    {
	        neighbour_sum += mgrid[i]*mgrid[mmap[4][i]];
	    }
	}
    }
    
    return -(mJ*neighbour_sum + mB*spin_sum);
}

double HB::Ising::calcEnergy(dim3 index)
{
    int spin_sum = mgrid[random_point];
    int neighbour_sum = 0;
   
    neighbour_sum += mgrid[random_point]*mgrid[mmap[0][random_point]]; //mid point
    neighbour_sum += mgrid[random_point]*mgrid[mmap[1][random_point]];
    if(dim > 1)
    {
        neighbour_sum += mgrid[random_point]*mgrid[mmap[2][random_point]];
        neighbour_sum += mgrid[random_point]*mgrid[mmap[3][random_point]];
	if(dim > 2)
	{
            neighbour_sum += mgrid[random_point]*mgrid[mmap[4][random_point]];
            neighbour_sum += mgrid[random_point]*mgrid[mmap[5][random_point]];
	}
    }

    return 2*mJ*neighbour_sum + mB*spin_sum;
}

double HB::Ising::calcMagn() //store initial value for energy
{
}

void HB::Ising::MC_sweep(gsl_rng *rng) //metropolis algorithm
{
    unsigned int n = m_Nx*m_Ny*m_Nz;
    double new_E = 0.0;
 
    #pragma omp parallel
    {

	
        double dE;
	
        #pragma omp for reduction(+:new_E)   	
        for(unsigned int j = 0; j < n; ++j)
	{
    	    dE = delta_energy(j);
	  
            if(dE < 0 || exp(-m_beta*dE) > gsl_rng_uniform(rng))
            {
	        flip(j, rng);
                new_E += dE;
            }
        }
    }
    
    mE += new_E;
}
