#pragma once

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "grid/grid.h"

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

/*! \def CUDA_CALL(x)
    \brief A macro that executes x and returns if it was not a cuda success. 
    An error message is printed in case.
*/
#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return;}} while(0)

/*! \def CUDA_CALL(x)
    \brief A macro that executes x and returns if it was not a curand success. 
    An error message is printed in case.
*/
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return;}} while(0)


namespace HB
{
  //! An abstract class representing the model describing the spins. A model class has to derive from it an override the pure virtual functions.
  class Model
  {
   public:
    Model() = default;
    /****************************************************************************//**
    * \brief constructs a model describing spins on a grid
    *
    * \param J Coupling constant
    * \param B three dimensional magnetic flux density with x component parallel 
    * to y direction of the grid and y component vice versa
    ********************************************************************************/
    Model(float J, float3 B) : mJ(J), mB(B) {};

    virtual ~Model() = default;
    
    /****************************************************************************//**
    * \brief calculates the overall energy of the underlying grid according 
    * to the model
    ********************************************************************************/
    virtual float calcEnergy(Grid<short int> &grid) const = 0;

    /****************************************************************************//**
    * \brief calculates the overall magnetization of the underlying grid according 
    * to the model
    ********************************************************************************/
    virtual float calcMagnetization(Grid<short int> &grid) const = 0;
    
    /****************************************************************************//**
    * \brief returns the magnetic flux density of the model 
    ********************************************************************************/
    float3 getBField() const {return mB;}

    /****************************************************************************//**
    * \brief returns the coupling constant of the model 
    ********************************************************************************/
    float getCouplingConstant() const {return mJ;}
    
    /****************************************************************************//**
    * \brief sets the magnetic flux density of the model 
    *
    * \param B Three dimensional magnetic flux density with x component parallel 
    * to y direction of the grid and y component vice versa
    ********************************************************************************/ 
    void setBField(const float3 B) {mB = B;}

    /****************************************************************************//**
    * \brief sets the coupling constant of the model 
    *
    * \param J Coupling constant
    ********************************************************************************/
    void setCouplingConstant(const float J) {mJ = J;}

    /****************************************************************************//**
    * \brief returns the energy after one spin is flipped
    *
    * \param index Index of the flipped spin
    ********************************************************************************/
    virtual float calcEnergy(Grid<short int> &grid, const dim3 index) const = 0;

    /****************************************************************************//**
    * \brief flips randomly a spin
    *
    * The function should define how a spin has to change 
    * for a Monte-Carlo accept-reject step. 
    *
    * \param index Index of the spin to be flipped
    * \param rng Pointer to the state of GSL RNG
    ********************************************************************************/
    virtual void flip(const dim3 index, Grid<short int> &grid) const = 0;
    
    float3 mB{0.0f,0.0f,0.0f}; /**< absolute value of an external, homogenous magnetic field */
    float mJ{1.0f}; /**< coupling constant */
};

    using Models = thrust::host_vector<std::unique_ptr<Model>>;
}
