#pragma once

#include "model/model.h"
#include "simulation.h"

namespace HB::Ising
{  
class Ising : public Model
{
 public:  
  Ising(float J, float3 B) : Model(J, B) {};

  /****************************************************//**
  * \brief Computes the overall energy of the entire grid 
  ********************************************************/
  float calcEnergy(Grid<short int> &grid) const final;

  /***********************************************************//**
  * \brief Computes the overall magnetization of the entire grid 
  ***************************************************************/
  float calcMagnetization(Grid<short int> &grid) const final;

  /*******************************************************************//**
  * \brief Computes the energy difference, when random_point is flipped.  
  ***********************************************************************/
  float calcEnergy(Grid<short int> &grid, const dim3 index) const final;
  
  /*********************************************************************//**
  * \brief Computes the energy difference, when a random point is flipped.  
  *************************************************************************/
  void flip(const dim3 index, Grid<short int> &grid) const final;
};
}
