#pragma once

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "include/grid/grid.h"

namespace HB
{
  template <class T>
  class Model
  {
   public:
    /****************************************************************************//**
  * \brief Constructs a Model.                                          
  *                                                                               
  * Initializes a neighbour table and reads the inverse \param beta, \param J, \param B, 
  * \param Nx, \param Ny and \param Nz 
  ********************************************************************************/
    Model(double beta, double J, double B, unsigned int Nx, unsigned int Ny, unsigned int Nz)
    {    
      /*m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    if(m_Nx == 0 || m_Ny == 0 || m_Nz == 0)
    {
        std::cerr << "Volume is zero. Exiting...\n";
	exit(1);
    }
    
    dim = 0;
    if(m_Nx != 1) dim++;
    if(m_Ny != 1) dim++;
    if(m_Nz != 1) dim++;
    if(dim == 2 && m_Nx == 1)
    {
        m_Nx = m_Nz;
        m_Nz = 1;
    }
    if(dim == 1 && m_Ny != 1)
    {
        m_Nx = m_Ny;
        m_Ny = 1;
    }
    if(dim == 1 && m_Nz != 1)
    {
        m_Nx = m_Nz;
        m_Nz = 1;
    }
    
    m_beta = beta;
    if(m_beta <= 0.0)
    {
        std::cerr << "Inverse temperature can not be zero or negative. Exiting...\n";
	exit(1);
    }
    mJ = J;
    mB = B;
    
    mGrid->neighbour_tab(); //defines a neighbour table*/
    };

    //functions
    virtual double calcEnergy() = 0;
    virtual double calcMagn() = 0;
    dim3 getBfield() const {return mB;};
    Grid<T>* getGrid() const {return mGrid;};
    T getSpin() const {return mSpin;};
    
   protected:
    virtual double calcEnergy(dim3 index) = 0;
    __device__ virtual double calcEnergy(short int*, double*, double*, unsigned int) = 0;
    virtual void flip(dim3 index, gsl_rng* rng) = 0;
    __device__ virtual void flip(short int*, unsigned int) = 0;
    virtual void MCSweep() = 0;
    __global__ virtual void cuda_MCSweep() = 0;
  
    dim3 mB; ///< The absolute value of an external, homogenous magnetic field
    double mJ; ///< The coupling constant
    T *mSpin;
    Grid<T> *mGrid;
};
}
