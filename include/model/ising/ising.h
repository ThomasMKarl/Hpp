#pragma once

#include "model/model.h"

namespace HB::Ising
{
class Ising : public Model<Spin<short int>>
{
 public:
  Ising(double J, float3 B) : Model(J, B) {};
  Ising(double J, float3 B, Grid<Spin<short int>> *grid) : Model(J, B, grid) {};
  Ising(double J, float3 B, Grid<Spin<short int>> *grid, bool hc) : Model(J, B, grid)
  {
    if(hc) grid->hotStart();
    else   grid->coldStart();
  };
  /****************************************************//**
  * \brief Computes the overall energy of the entire grid 
  ********************************************************/
  double calcEnergy() final;
  double calcMagn() final;

 private:
  /******************************************************************//**
  * \brief Computes the energy difference, when random_point is flipped.  
  **********************************************************************/
  double calcEnergy(dim3 index) final;
  //__device__ double calcEnergy(short int*, double*, double*, unsigned int) final;
  
  /******************************************************************//**
  * \brief Computes the energy difference, when a random point is flipped.  
  **********************************************************************/
  void flip(dim3 index, gsl_rng *rng) final
  {
    rng = rng;
    Spin<short int> s;
    s.x = -this->mGrid->getSpin(index).x;
    s.y = 0;
    this->mGrid->setSpin(index, s);
  };
  //__device__ void flip(short int *grid, unsigned int grid_point, curandState_t *state) final {state = state; grid[grid_point] *= -1;};

};
}
