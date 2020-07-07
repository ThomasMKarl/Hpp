#include "xy.h"

/****************************************************//**
* \brief Computes the overall energy of the entire grid 
********************************************************/
double XY::compute_energy() //store initial value for energy
{
    unsigned int n = m_Nx*m_Ny*m_Nz;
    
    int neighbour_sum = 0;
    int spin_sum = 0;
    
    #ifndef _OPENMP
    unsigned int l = 0;
    unsigned int e = n;
    #endif

    #pragma omp parallel
    {	
	double x, x1, y, y1;
	
        #pragma omp for reduction(+:spin_sum,neighbour_sum)
	for(unsigned int i = 0; i < n; ++i)
        {
	    spin_sum += mgrid[i];
	  
            x = cos(mgrid[i]); y = sin(mgrid[i]);
	   
	    x1 = cos(mgrid[mmap[1][i]]); y1 = sin(mgrid[mmap[1][i]]);
            neighbour_sum += x*x1 + y*y1;

	    if(dim > 1)
	    {
	        x1 = cos(mgrid[mmap[2][i]]); y1 = sin(mgrid[mmap[2][i]]);
                neighbour_sum += x*x1 + y*y1;

		if(dim > 2)
		{
	            x1 = cos(mgrid[mmap[4][i]]); y1 = sin(mgrid[mmap[4][i]]);
	            neighbour_sum += x*x1 + y*y1;
		}
	    }
        }
    }
 
    return -(mJ*neighbour_sum + mB*spin_sum);
}

/******************************************************************//**
* \brief Computes the energy difference, when randm_point is flipped.  
**********************************************************************/
double XY::delta_energy(unsigned int random_point)
{ 
    int neighbour_sum = 0;
    double x, x1, y, y1;


    x = cos(mgrid[random_point]); y = sin(mgrid[random_point]);
    
    x1 = cos(mgrid[mmap[0][random_point]]); y1 = sin(mgrid[mmap[0][random_point]]);
    neighbour_sum += x*x1 + y*y1;
    
    x1 = cos(mgrid[mmap[1][random_point]]); y1 = sin(mgrid[mmap[1][random_point]]);
    neighbour_sum += x*x1 + y*y1;

    if(dim > 1)
    {
        x1 = cos(mgrid[mmap[2][random_point]]); y1 = sin(mgrid[mmap[2][random_point]]);
        neighbour_sum += x*x1 + y*y1;
    
        x1 = cos(mgrid[mmap[3][random_point]]); y1 = sin(mgrid[mmap[3][random_point]]);
        neighbour_sum += x*x1 + y*y1;

        if(dim > 2)
	{
            x1 = cos(mgrid[mmap[4][random_point]]); y1 = sin(mgrid[mmap[4][random_point]]);
            neighbour_sum += x*x1 + y*y1;
    
            x1 = cos(mgrid[mmap[5][random_point]]); y1 = sin(mgrid[mmap[5][random_point]]);
            neighbour_sum += x*x1 + y*y1;
	}
    }

    
    return 2*mJ*neighbour_sum - mB*mgrid[random_point];
}

/******************************************************************//**
* \brief Computes the energy difference, when randm_point is flipped.  
**********************************************************************/
void XY::flip(unsigned int point, gsl_rng *rng)
{
    double randn = M_PI*(gsl_rng_uniform(rng) - 1.5);
  
    mgrid[point] += randn;
    if(mgrid[point] < 0) mgrid[point] += 2.0*M_PI;
}

/*****************************************************************//**
* \brief Performs one Monte-Carlo step, updating the entire grid and 
* the energy \param mE of the grid. 
*********************************************************************/
void Model::MC_sweep(gsl_rng *rng) //metropolis algorithm
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
