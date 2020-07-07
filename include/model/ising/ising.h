#ifndef ISING_H
#define ISING_H

#include "include/model/model.h"

namespace HB::Ising
{
class Ising : public Model
{
 public:
  Ising(double beta, double J, double B, unsigned int Nx, unsigned int Ny, unsigned int Nz, bool hc = true) : Model(beta, J, B, Nx, Ny, Nz)
  { 
    if(hc) hot_start();
    else  cold_start();
  };
  /****************************************************//**
  * \brief Computes the overall energy of the entire grid 
  ********************************************************/
  double calcEnergy() final;
  double calcMagn() final;

 private:
  /******************************************************************//**
  * \brief Computes the energy difference, when randm_point is flipped.  
  **********************************************************************/
  double calcEnergy(dim3 index) final;
  __device__ double calcEnergy(short int*, double*, double*, unsigned int) final;
  
  /******************************************************************//**
  * \brief Computes the energy difference, when randm_point is flipped.  
  **********************************************************************/
  void flip(unsigned int grid_point, gsl_rng *rng) final {rng = rng; mgrid[grid_point] *= -1;};
  __device__ void flip(short int *grid, unsigned int grid_point, state) final {state = state; grid[grid_point] *= -1;};
  /*****************************************************************//**
  * \brief Performs one Monte-Carlo step, updating the entire grid and 
  * the energy \param mE of the grid. 
  *********************************************************************/
  void MCSweep() final;
  __global__ void cuda_MCSweep() final;
};
}
