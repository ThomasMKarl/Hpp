#include "include/grid/grid.h"

void HB::hotStart() //all spins 1
{
    unsigned int n = mGridSize.x*mGridSize.y*mGridSize.z;
    
    mGrid.resize(n);
    
    #pragma omp parallel
    {
	gsl_rng *rng;
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        long seed = time(NULL);
        gsl_rng_set(rng,seed);

        #pragma omp for
        for(unsigned int i = 0; i < n; ++i)
	{
	    mGrid[i].x = 2.0*M_PI*gsl_rng_uniform(rng);
	    mGrid[i].y = 2.0*M_PI*gsl_rng_uniform(rng);
	    mGrid[i].z = 2.0*M_PI*gsl_rng_uniform(rng);
	}
    }
}

void HB::coldStart() //all spins 1
{
    unsigned int n = m_Nx*m_Ny*m_Nz;
    
    mgrid_phi.resize(n);
    mgrid_theta.resize(n);
    
    #pragma omp parallel for
    for(unsigned int i = 0; i < n; ++i)
    {
	mgrid_phi[i]   = 0.0;
	mgrid_theta[i] = 0.5*M_PI;
    }   
    
    mE = compute_energy(); //store initial value for energy
}

void HB::neighbourTab()
{
    unsigned int n = m_Nx*m_Ny*m_Nz;
    
    mTable.resize(n);
        
    #pragma omp parallel
    {	
        unsigned int i;
	unsigned int j = 0;
	unsigned int k = 0;
	
	#pragma omp for
        for(unsigned int h = 0; h < n; ++h)
        {
            i = h%m_Nx;
	    if(dim > 1)
	    {
                j = h/m_Nx%m_Ny;
	        if(dim > 2)
	        {
                    k = h/(m_Nx*m_Ny);
	        }
	    }
            //std::cout << h << " " << i << " " << j << " " << k << "\n";
	    
            if(i == 0)      mmap[0][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*j + m_Nx - 1;
            else            mmap[0][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*j + i    - 1;

            if(i == m_Nx-1) mmap[1][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*j;
            else            mmap[1][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*j + i    + 1;

	    if(dim > 1)
	    {
                if(j == 0)      mmap[2][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*(m_Ny-1) + i;
                else            mmap[2][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*(j-1)    + i;
	
                if(j == m_Ny-1) mmap[3][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + i;
                else            mmap[3][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*k + m_Nx*(j+1) + i;

	        if(dim > 2)
	        {
                    if(k == 0)      mmap[4][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*(m_Nz-1) + m_Nx*j + i;
                    else            mmap[4][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*(k-1)    + m_Nx*j + i;
	
                    if(k == m_Nz-1) mmap[5][m_Ny*m_Nx*k + m_Nx*j + i] = m_Nx*j + i;
                    else            mmap[5][m_Ny*m_Nx*k + m_Nx*j + i] = m_Ny*m_Nx*(k+1) + m_Nx*j + i;
	        }
	    }

	//std::cout << mmap[0][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[1][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[2][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[3][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[4][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[5][m_Ny*m_Nx*k + m_Nx*j + i]<<"\n\n";
        }
    }
}

void HB::saveGrid(std::string path)
{
}
