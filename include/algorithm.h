#pragma once

#include "include/model/model.h"
namespace HB
{
/********************************************************************************************//**
* \brief Performs a production-run.                                                             
*                                                                                               
* A production run consists of steps updates producing steps correlated, unthermalised energies 
************************************************************************************************/
template <typename T, class M> void productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, bool cuda);

/***********************************************************************************//**
* \brief Performs a production-run.                                                    
*                                                                                      
* A production run consists of steps updates producing steps correlated, unthermalised 
* energies and magnetizations                                                          
***************************************************************************************/
template <typename T, class M> void productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, bool cuda);

/*******************************************************************************************//**
* \brief Performs a production-run via simulated annealing.                                     
*                                                                                              
* The run starts at Temperature T = T_begin and ends with T = T_end. In steps of T_steps       
* a production run of steps updates computes the initial state for the next run. In the last   
* one energies and magnetizations are computed. Simulated annealing consists of                
* n = steps*(T_begin-T_end/T_step) updates producing steps correlated, unthermalised energies. 
* For each step the member variable m_beta = 1/T is gets overridden.                           
***********************************************************************************************/
template <typename T, class M> void productionRun(M &model, size_t V, size_t steps, std::vector<T> &energies, std::vector<T> &magn, T T_begin, T T_end, T T_step, bool cuda);

template <typename T> T mean(std::vector<T> &o);
template <typename T> T covFunc(size_t t, std::vector<T> &o);
template <typename T> T intAuto(std::vector<T> &o);
template <typename T> T error(std::vector<T> &o);
template <typename T> T blockingError(std::vector<T> &o, size_t block_number);
template <typename T> T bootstrapError(std::vector<T> &o, size_t sample_size, size_t sample_number, T tau);
template <typename T> T errorProp(std::vector<T> &o, T tau);

/**************************************************************************//**
* \brief Calculates the statistical error on energies using error propagation
* for a give autocorrelation time tau                                         
******************************************************************************/
template <typename T> T errorPropEnergy(std::vector<T> &energies, T beta, T tau, size_t V);

/***************************************************************************************//**
* \brief Calculates the statistical error on magnetizations (magn) using error propagation 
* for a give autocorrelation time tau                                                      
*******************************************************************************************/
template <typename T> T errorPropMagnetization(std::vector<T> &magn, T beta, T tau, size_t V);
template <typename T> size_t d2i(T d);
template <typename T> void removeCorr(std::vector<T> &o);

template <typename T> T specificHeat(std::vector<T> &energies, T beta, size_t V);
template <typename T> T magnSusz(std::vector<T> &magn, T beta, size_t V);
}
  
