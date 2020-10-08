#include "hb.h"

/****************************************************//**
* \brief Computes the overall energy of the entire grid 
********************************************************/
double Heisenberg::compute_energy() //store initial value for energy
{   
    unsigned int n = m_Nx*m_Ny*m_Nz;
    
    int neighbour_sum = 0;
    int spin_sum = 0;
	
    #pragma omp parallel
    {	
	double x, x1, y, y1, z, z1;
	
        #pragma omp parallel for reduction(+:spin_sum,neighbour_sum)
	for(unsigned int i = 0; i < n; ++i)
        {
	    spin_sum += mgrid_theta[i]*mgrid_phi[i];
	  
	    x = sin(mgrid_theta[i])*cos(mgrid_phi[i]);
            y = sin(mgrid_theta[i])*sin(mgrid_phi[i]);
            z = cos(mgrid_theta[i]);
	    
	    x1 = sin(mgrid_theta[mmap[1][i]])*cos(mgrid_phi[mmap[1][i]]);
            y1 = sin(mgrid_theta[mmap[1][i]])*sin(mgrid_phi[mmap[1][i]]);
            z1 = cos(mgrid_theta[mmap[1][i]]);
            neighbour_sum += x*x1 + y*y1 + z*z1;

	    if(dim > 1)
	    {
	        x1 = sin(mgrid_theta[mmap[2][i]])*cos(mgrid_phi[mmap[2][i]]);
                y1 = sin(mgrid_theta[mmap[2][i]])*sin(mgrid_phi[mmap[2][i]]);
                z1 = cos(mgrid_theta[mmap[2][i]]);
                neighbour_sum += x*x1 + y*y1 + z*z1;

	        if(dim > 2)
	        {	    
                    x1 = sin(mgrid_theta[mmap[4][i]])*cos(mgrid_phi[mmap[4][i]]);
                    y1 = sin(mgrid_theta[mmap[4][i]])*sin(mgrid_phi[mmap[4][i]]);
                    z1 = cos(mgrid_theta[mmap[4][i]]);
	            neighbour_sum += x*x1 + y*y1 + z*z1;
	        }
	    }
        }
    }
    
    return -(mJ*neighbour_sum + mB*spin_sum);
}

/******************************************************************//**
* \brief Computes the energy difference, when randm_point is flipped.  
**********************************************************************/
double Heisenberg::delta_energy(unsigned int random_point)
{  
    int neighbour_sum = 0;
    double x, x1, y, y1, z, z1;

 
    x  = sin(mgrid_theta[random_point])*cos(mgrid_phi[random_point]);
    y  = sin(mgrid_theta[random_point])*sin(mgrid_phi[random_point]);
    z  = cos(mgrid_theta[random_point]);


    x1 = sin(mgrid_theta[mmap[0][random_point]])*cos(mgrid_phi[mmap[0][random_point]]);
    y1 = sin(mgrid_theta[mmap[0][random_point]])*sin(mgrid_phi[mmap[0][random_point]]);
    z1 = cos(mgrid_theta[mmap[0][random_point]]);
    neighbour_sum += x*x1 + y*y1 + z*z1;
    
    x1 = sin(mgrid_theta[mmap[1][random_point]])*cos(mgrid_phi[mmap[1][random_point]]);
    y1 = sin(mgrid_theta[mmap[1][random_point]])*sin(mgrid_phi[mmap[1][random_point]]);
    z1 = cos(mgrid_theta[mmap[1][random_point]]);
    neighbour_sum += x*x1 + y*y1 + z*z1;

    if(dim > 1)
    {
        x1 = sin(mgrid_theta[mmap[2][random_point]])*cos(mgrid_phi[mmap[2][random_point]]);
        y1 = sin(mgrid_theta[mmap[2][random_point]])*sin(mgrid_phi[mmap[2][random_point]]);
        z1 = cos(mgrid_theta[mmap[2][random_point]]);
        neighbour_sum += x*x1 + y*y1 + z*z1;
    
        x1 = sin(mgrid_theta[mmap[3][random_point]])*cos(mgrid_phi[mmap[3][random_point]]);
        y1 = sin(mgrid_theta[mmap[3][random_point]])*sin(mgrid_phi[mmap[3][random_point]]);
        z1 = cos(mgrid_theta[mmap[3][random_point]]);
        neighbour_sum += x*x1 + y*y1 + z*z1;

	if(dim > 2)
	{
            x1 = sin(mgrid_theta[mmap[4][random_point]])*cos(mgrid_phi[mmap[4][random_point]]);
            y1 = sin(mgrid_theta[mmap[4][random_point]])*sin(mgrid_phi[mmap[4][random_point]]);
            z1 = cos(mgrid_theta[mmap[4][random_point]]);
            neighbour_sum += x*x1 + y*y1 + z*z1;
    
            x1 = sin(mgrid_theta[mmap[5][random_point]])*cos(mgrid_phi[mmap[5][random_point]]);
            y1 = sin(mgrid_theta[mmap[5][random_point]])*sin(mgrid_phi[mmap[5][random_point]]);
            z1 = cos(mgrid_theta[mmap[5][random_point]]);
            neighbour_sum += x*x1 + y*y1 + z*z1;
	}
    }
    

    return 2*mJ*neighbour_sum /*- mB*mgrid[random_point]*/;
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
