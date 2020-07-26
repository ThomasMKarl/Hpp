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
  template <class T>
  class Model
  {
   public:
    /****************************************************************************//**
    * \brief constructs a model describing spins on a grid
    *
    * \param J Coupling constant
    * \param B three dimensional magnetic flux density with x component parallel 
    * to y direction of the grid and y component vice versa
    ********************************************************************************/
    Model(double J, float3 B)
    {
        mJ = J;
	mB = B;
    };

    /****************************************************************************//**
    * \brief constructs a model describing Spins on a grid
    *
    * The grid is NOT referenced. The user has to control an instance of a grid class.
    * The model has exactly one grid, but a grid can be tied to any model.
    * When the underlying grid is destroyed, the model has to find a new one to work 
    * properly.
    *
    * \param J Coupling constant
    * \param B Three dimensional magnetic flux density with x component parallel 
    * to y direction of the grid and y component vice versa
    * \param grid Points to the grid on which the spins of the model live
    ********************************************************************************/
    Model(double J, float3 B, Grid<T> *grid)
    {
        mJ = J;
	mB = B;
	mGrid = grid;
    };
    
    /****************************************************************************//**
    * \brief calculates the overall energy of the underlying grid according 
    * to the model
    ********************************************************************************/
    virtual double calcEnergy() = 0;

    /****************************************************************************//**
    * \brief calculates the overall magnetization of the underlying grid according 
    * to the model
    ********************************************************************************/
    virtual double calcMagn() = 0;
    
    /****************************************************************************//**
    * \brief returns the magnetic flux density of the model 
    ********************************************************************************/
    float3 getBfield() const {return mB;}

    /****************************************************************************//**
    * \brief returns the coupling constant of the model 
    ********************************************************************************/
    double getCoupling() const {return mJ;}
    
    /****************************************************************************//**
    * \brief returns the pointer to the underlying grid 
    ********************************************************************************/
    Grid<T>* getGrid() const {return mGrid;}
    
    /****************************************************************************//**
    * \brief sets the magnetic flux density of the model 
    *
    * \param B Three dimensional magnetic flux density with x component parallel 
    * to y direction of the grid and y component vice versa
    ********************************************************************************/ 
    void setBfield(float3 B) {mB = B;}

    /****************************************************************************//**
    * \brief sets the coupling constant of the model 
    *
    * \param J Coupling constant
    ********************************************************************************/
    void setCoupling(double J) {mJ = J;}

    /****************************************************************************//**
    * \brief sets the pointer to the grid explicitly
    *
    * \param grid Pointer to a grid holding the correct type of spins
    ********************************************************************************/
    void setGrid(Grid<T> *grid) {mGrid = grid;}
    
   protected:
    /****************************************************************************//**
    * \brief returns the energy after one spin is flipped
    *
    * \param index Index of the flipped spin
    ********************************************************************************/
    virtual double calcEnergy(dim3 index) = 0;
    
    //__device__ virtual double calcEnergy(short int*, double*, double*, unsigned int) = 0;

    /****************************************************************************//**
    * \brief flips randomly a spin
    *
    * The function should define how a spin has to change 
    * for a Monte-Carlo accept-reject step. 
    *
    * \param index Index of the spin to be flipped
    * \param rng Pointer to the state of GSL RNG
    ********************************************************************************/
    virtual void flip(dim3 index, gsl_rng* rng) = 0;
    
    //__device__ virtual void flip(short int*, unsigned int) = 0;
  
    float3 mB; /**< absolute value of an external, homogenous magnetic field */
    double mJ; /**< coupling constant */
    Grid<T> *mGrid; /**< points to a grid holding the correct spin type */
};
}
