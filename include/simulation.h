#pragma once

#include <QtGui/QPainter>

#include <random>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include "model/ising/ising.h"

namespace HB
{
//__global__ void setupKernel(curandState *state, size_t seed, size_t n);
  
class MetropolisSimulation
{
 public:
  MetropolisSimulation(const Grid<short int> &grid, const float beta, const size_t steps) : mBeta(beta), mSteps(steps)
  {
    mGrid = std::make_shared<Grid<short int>>(grid);
  }
  
  template<class M>
  void operator()(const M& model) const
  {
    simulate(model);
  }
 private:
  float mBeta{1.0f};
  size_t mSteps{1000};
  std::shared_ptr<Grid<short int>> mGrid{};
  std::vector<float> mEnergies{};
  std::vector<float> mMagnetization{};
  
/*************************************************************************************************************//**
* \brief performs a production-run.                                                             
*                                                                                               
* A production run consists of a number of MC sweeps producing correlated energies 
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
*****************************************************************************************************************/
  template<class M>
  void simulate(const M &model);

/********************************************************************//**
* \brief performs one Monte-Carlo step, updating the entire grid 
* associated with a certain model and calculates the energy difference
*
* \param Reference to the model
* \param beta Inverse Temperature
* \param rng Pointer to the state of the GSL RNG
* \return Energy difference
************************************************************************/
  template <class M>
  float metropolisSweep(const M &model);
};
  
class MetropolisSimulationQt
{
 public:
  MetropolisSimulationQt(const Grid<short int> &grid, const float beta) : mBeta(beta)
  {
    mGrid = std::make_shared<Grid<short int>>(grid);
  }
  
  template <class M>
  void operator()(const M &model) const
  {
    simulate(model);
  }
 private:
  float mBeta{1.0f};
  size_t mSteps{1000};
  std::shared_ptr<Grid<short int>> mGrid{};
  
/*************************************************************************************************************//**
* \brief performs a production-run.                                                             
*                                                                                               
* A production run consists of a number of MC sweeps producing correlated energies 
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
*****************************************************************************************************************/
  template <class M>
  void simulate(const M &model);

/********************************************************************//**
* \brief performs one Monte-Carlo step, updating the entire grid 
* associated with a certain model and calculates the energy difference
*
* \param Reference to the model
* \param beta Inverse Temperature
* \param rng Pointer to the state of the GSL RNG
* \return Energy difference
************************************************************************/
  template <class M>
  float metropolisSweep(const M &model);
};
  
class SimulatedAnnealing
{
 public:
  SimulatedAnnealing(const Grid<short int> &grid, const float beta, const size_t steps, const float3 temperatureSteps) : mBeta(beta), mSteps(steps), mTemperatureSteps(temperatureSteps)
  {
    mGrid = std::make_shared<Grid<short int>>(grid);
  }
  
  template<class M>
  void operator()(const M& model) const
  {
    simulate(model);
  }
 private:
  float mBeta{1.0f};
  size_t mSteps{1000};
  float3 mTemperatureSteps{1.0f,0.0f,0.1f};
  std::shared_ptr<Grid<short int>> mGrid{};
  std::vector<float> mEnergies{};
  std::vector<float> mMagnetization{};
  
/*************************************************************************************************************************//**
* \brief Performs a production-run via simulated annealing.                                     
*                                                                                              
* The run starts at Temperature T = T_begin and ends with T = T_end. In steps of T_steps       
* a production run computes the initial state for the next run. In the last   
* one energies and magnetizations are computed. Simulated annealing consists of                
* \f$n = \text{steps}\times \text{floor}\left(\frac{T{\_}\text{begin}-T{\_}\text{end}}{T{\_}\text{step}}\right)\f$ updates producing correlated energies and magnetizations.
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
* \param magnetizations Gets overridden with a standard vector of correlated magnetizations (initial and one for each step)
* \param T_begin Initial temperature value
* \param T_end Minimum temperature value 
* \param T_step Step size from one temperature to another
* \param cuda Set to true foe GPU acceleration
*****************************************************************************************************************************/
  template<class M>
  void simulate(const M &model);
};

///////////////////////////////////////////

/***********************************************************************************//**
* \brief calculates the mean value of a standard vector \f$\mu = \sum_i o_i\f$
*
* \param o Reference to the standard vector
* \param return Mean value
***************************************************************************************/
template <typename T> T mean(const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the covariance function up to a certain index of a standard vector
*
* \param t Index in vector
* \param o Reference to the standard vector
* \return Value of the covariance function
***************************************************************************************/
template <typename T> T covFunc(const size_t t, const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the integrated autocorrelation time of a standard vector
*
* \param o Reference to the standard vector
* \return Integrated autocorrelation time
***************************************************************************************/
template <typename T> T intAuto(const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the standard error (standard deviation) of a standard vector
*
* \param o Reference to the standard vector
* \return Standard error
***************************************************************************************/
template <typename T> T error(const std::vector<T> &o);
  
/***********************************************************************************//**
* \brief calculates the blocking error of a standard vector
*
* \param block_number Number of blocks the vector is divided into
* \param o Reference to the standard vector
* \return Blocking error
***************************************************************************************/
template <typename T> T blockingError(const std::vector<T> &o, const size_t block_number);

/***********************************************************************************//**
* \brief calculates the bootstrap error of a standard vector
*
* \param o Reference to the standard vector
* \param sample_size Size of the bootstrap samples
* \parma sample_number Number of resamplings
* \param tau Integrated autocorrelation time
* \return Bootstrap error
***************************************************************************************/
template <typename T> T bootstrapError(const std::vector<T> &o, const size_t sample_size, const size_t sample_number, const T tau);

/***********************************************************************************//**
* \brief calculates the error based on error propagation of a standard vector
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time
* \return Error propagation
***************************************************************************************/
template <typename T> T errorProp(const std::vector<T> &o, const T tau);

/**************************************************************************//**
* \brief calculates the statistical error on energies using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time 
* \return Error propagation                                    
******************************************************************************/
template <typename T> T errorPropEnergy(const std::vector<T> &energies, const T beta, const T tau, const size_t V);

/***************************************************************************************//**
* \brief calculates the statistical error on magnetizations (magn) using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time   
* \return Error propagation                                                   
*******************************************************************************************/
template <typename T> T errorPropMagnetization(const std::vector<T> &magn, const T beta, const T tau, const size_t V);
  
/***************************************************************************************//**
* \brief rounds a floating point value to an integer
*
* \param d Value to be rounded
* \return Resulting integer
*******************************************************************************************/
template <typename T> size_t d2i(const T d);

/***************************************************************************************//**
* \brief removes correlation from a standard vector
*
* The function calculates the integrated autocorrelation time \f$\tau\f$. The thermalisation 
* of \f$20\tau\f$ is removed. The integrated autocorrelation time is calculated again in order 
* to remove values from the input vector until the remaining ones are uncorrelated. The 
* The function overrides the input vector.
*
* \param o Standard vector of correlated values
*******************************************************************************************/
template <typename T> void removeCorr(const std::vector<T> &o);

/****************************************************************************//**
* \brief computes the specific heat per Volume as secondary quantity for a 
* given vector of energies
*
* \param Reference to a standard vector of energies
* \param beta Inverse temperature
* \param V Volume
* \return Specific heat
********************************************************************************/
template <typename T> T specificHeat(const std::vector<T> &energies, const T beta, const size_t V);

/****************************************************************************//**
* \brief computes the magnetic suszeptibility per Volume as secondary quantity 
* for a given vector of magnetizations
*
* \param Reference to a standard vector of magnetizations
* \param beta Inverse temperature
* \param V Volume
* \return Magnetic suzeptibility
********************************************************************************/
template <typename T> T magnSusz(const std::vector<T> &magn, const T beta, const size_t V);

using Models = std::vector<std::unique_ptr<HB::Model>>;

void simulate(const Models &models);
}
