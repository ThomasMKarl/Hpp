#pragma once

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include "model/model.h"

namespace HB
{
  //__global__ void setupKernel(curandState *state, size_t seed, size_t n);
  
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
* \param cuda Set to true foe GPU acceleration
*****************************************************************************************************************/
template <typename T, class M> void productionRun(M &model, T beta, size_t V, size_t steps, std::vector<T> &energies, bool cuda);

/*************************************************************************************************************************//**
* \brief performs a production-run.                                                             
*                                                                                               
* A production run consists of a number of MC sweeps producing correlated energies 
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
* \param magnetizations Gets overridden with a standard vector of correlated magnetizations (initial and one for each step)
* \param cuda Set to true foe GPU acceleration
*****************************************************************************************************************************/
template <typename T, class M> void productionRun(M &model, T beta, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, bool cuda);

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
template <typename T, class M> void productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, T T_begin, T T_end, T T_step, bool cuda);
  
/***********************************************************************************//**
* \brief calculates the mean value of a standard vector \f$\mu = \sum_i o_i\f$
*
* \param o Reference to the standard vector
* \param return Mean value
***************************************************************************************/
template <typename T> T mean(std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the covariance function up to a certain index of a standard vector
*
* \param t Index in vector
* \param o Reference to the standard vector
* \return Value of the covariance function
***************************************************************************************/
template <typename T> T covFunc(size_t t, std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the integrated autocorrelation time of a standard vector
*
* \param o Reference to the standard vector
* \return Integrated autocorrelation time
***************************************************************************************/
template <typename T> T intAuto(std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the standard error (standard deviation) of a standard vector
*
* \param o Reference to the standard vector
* \return Standard error
***************************************************************************************/
template <typename T> T error(std::vector<T> &o);
  
/***********************************************************************************//**
* \brief calculates the blocking error of a standard vector
*
* \param block_number Number of blocks the vector is divided into
* \param o Reference to the standard vector
* \return Blocking error
***************************************************************************************/
template <typename T> T blockingError(std::vector<T> &o, size_t block_number);

/***********************************************************************************//**
* \brief calculates the bootstrap error of a standard vector
*
* \param o Reference to the standard vector
* \param sample_size Size of the bootstrap samples
* \parma sample_number Number of resamplings
* \param tau Integrated autocorrelation time
* \return Bootstrap error
***************************************************************************************/
template <typename T> T bootstrapError(std::vector<T> &o, size_t sample_size, size_t sample_number, T tau);

/***********************************************************************************//**
* \brief calculates the error based on error propagation of a standard vector
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time
* \return Error propagation
***************************************************************************************/
template <typename T> T errorProp(std::vector<T> &o, T tau);

/**************************************************************************//**
* \brief calculates the statistical error on energies using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time 
* \return Error propagation                                    
******************************************************************************/
template <typename T> T errorPropEnergy(std::vector<T> &energies, T beta, T tau, size_t V);

/***************************************************************************************//**
* \brief calculates the statistical error on magnetizations (magn) using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time   
* \return Error propagation                                                   
*******************************************************************************************/
template <typename T> T errorPropMagnetization(std::vector<T> &magn, T beta, T tau, size_t V);
  
/***************************************************************************************//**
* \brief rounds a floating point value to an integer
*
* \param d Value to be rounded
* \return Resulting integer
*******************************************************************************************/
template <typename T> size_t d2i(T d);

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
template <typename T> void removeCorr(std::vector<T> &o);

/****************************************************************************//**
* \brief computes the specific heat per Volume as secondary quantity for a 
* given vector of energies
*
* \param Reference to a standard vector of energies
* \param beta Inverse temperature
* \param V Volume
* \return Specific heat
********************************************************************************/
template <typename T> T specificHeat(std::vector<T> &energies, T beta, size_t V);

/****************************************************************************//**
* \brief computes the magnetic suszeptibility per Volume as secondary quantity 
* for a given vector of magnetizations
*
* \param Reference to a standard vector of magnetizations
* \param beta Inverse temperature
* \param V Volume
* \return Magnetic suzeptibility
********************************************************************************/
template <typename T> T magnSusz(std::vector<T> &magn, T beta, size_t V);

/********************************************************************//**
* \brief performs one Monte-Carlo step, updating the entire grid 
* associated with a certain model and calculates the energy difference
*
* \param Reference to the model
* \param beta Inverse Temperature
* \param rng Pointer to the state of the GSL RNG
* \return Energy difference
************************************************************************/
template <typename T, class M> T MCSweep(M &model, T beta, gsl_rng *rng); //metropolis algorithm
}
