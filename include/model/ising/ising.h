#pragma once

#include "model/model.h"

namespace HB::Ising
{  
class Ising : public Model
{
 public:
  using SimulationStrategy = std::function<void(const Ising&)>;
  
  Ising(float J, float3 B, SimulationStrategy s) : Model(J, B), mS(s) {};

  /****************************************************//**
  * \brief Computes the overall energy of the entire grid 
  ********************************************************/
  float calcEnergy(const Grid<short int> &grid) final;

  /***********************************************************//**
  * \brief Computes the overall magnetization of the entire grid 
  ***************************************************************/
  float calcMagnetization(const Grid<short int> &grid) final;

 private:
  /*******************************************************************//**
  * \brief Computes the energy difference, when random_point is flipped.  
  ***********************************************************************/
  float calcEnergy(const Grid<short int> &grid, const dim3 index) final;
  
  /*********************************************************************//**
  * \brief Computes the energy difference, when a random point is flipped.  
  *************************************************************************/
  void flip(const dim3 index, Grid<short int> &grid) final;
  
  void simulate() const override
  {
    mS(*this);
  }
  
  SimulationStrategy mS;
};
}
