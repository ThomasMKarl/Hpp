#pragma once

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <hdf5.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace HB
{
  //! A struct containing six neighbour inidices of a spin in three dimensions
  typedef struct {
    size_t up;
    size_t right;
    size_t back;
    size_t down;
    size_t left;
    size_t front;
  }NB;

  //! A class storing spin values represented by two angles
  template<typename T>
  class Spin
  {
   public:
    T x; /**< polar angle */
    T y; /**< azimuthal angle */
  };
  
  //! A class storing spin values as a grid and certain meta information
  template<class T>
  class Grid
  {
   public:
    /*****************************************************************//**
    * \brief constructs a grid, allocates memory an initializes 
    * neighbour tab
    *
    * \param dim Dimension of the spins
    * \param gridSize Grid size
    *********************************************************************/
    Grid(short int dim, dim3 gridSize)
    {
      mDim  = dim;
      mGridSize = gridSize;

      mGrid.resize(gridSize.x*gridSize.y*gridSize.z);
      this->neighbourTab();
    };

    /*****************************************************************//**
    * \brief returns grid size
    *********************************************************************/
    dim3 getGridSize() const
    {
      return mGridSize;
    }

    /*****************************************************************//**
    * \brief sets grid size, resizes grid and computes neighbours
    *
    * \param gridSize Grid size
    *********************************************************************/
    void setGridSize(dim3 gridSize)
    {
      mGridSize = gridSize;

      mGrid.resize(gridSize.x*gridSize.y*gridSize.z);
      this->neighbourTab();
    }

    /*****************************************************************//**
    * \brief returns spin value of a certain grid index
    *
    * \param Grid index
    *********************************************************************/
    T getSpin(dim3 index) const
    {
      size_t idx = index.x*this->mGridSize.x +
	           index.y*this->mGridSize.x*this->mGridSize.y +
                   index.z;
      return this->mGrid[idx];
    }

    /*****************************************************************//**
    * \brief set specific spin to a certain value
    *
    * \param Grid index
    * \param Spin value
    *********************************************************************/
    void setSpin(dim3 index, T s)
    {
      size_t idx = index.x*this->mGridSize.x +
	           index.y*this->mGridSize.x*this->mGridSize.y +
                   index.z;
      this->mGrid[idx] = s;
    }

    /*****************************************************************//**
    * \brief returns neighbours of a certain grid index
    *
    * \param grid index
    *********************************************************************/
    NB getNeighbours(dim3 index) const
    {
      size_t idx = index.x*this->mGridSize.x +
	           index.y*this->mGridSize.x*this->mGridSize.y +
                   index.z;
      return this->mTable[idx];
    }

    /*****************************************************************//**
    * \brief returns spin dimension
    *********************************************************************/
    short int getDim() const {return mDim;}

    /*****************************************************************//**
    * \brief sets spin dimension
    *
    * \param dim Spin dimension
    *********************************************************************/
    void setDim(short int dim) { mDim = dim;}

    /*****************************************************************//**
    * \brief initializes all spins randomly between 0 and 2*pi  
    *********************************************************************/
    void hotStart();
    
    /*****************************************************************//**
    * \brief initializes all spins as up, 
    * i.e. mGrid.x[...] = 0.5*pi and mGrid.y[...] = 0                        
    *********************************************************************/
    void coldStart();
    
    /************************************************************************************************//**
    * \brief saves grid as hdf5 file
    *
    * \param path Path to file 
    ****************************************************************************************************/
    herr_t saveGrid(std::string path);
    
    /************************************************************************************************//**
    * \brief returns a host pointer to the grid data on the device
    ****************************************************************************************************/
    T *getDeviceGrid() const
    {
      return thrust::raw_pointer_cast(d_mGrid.data());
    }
    
    /************************************************************************************************//**
    * \brief returns a host pointer to the table data on the device
    ****************************************************************************************************/
    NB *getDeviceTable() const
    {
      return thrust::raw_pointer_cast(d_mTable.data());
    }
    
    /************************************************************************************************//**
    * \brief copies grid and table to device 
    ****************************************************************************************************/
    void upload()
    {
      d_mGrid = mGrid;
      d_mTable = mTable;
    }
    
    /************************************************************************************************//**
    * \brief copies grid to host
    ****************************************************************************************************/
    void download()
    {
      mGrid = mGrid;
    }
    
   private:
    /*******************************************************************//**
    * \brief generates a table of the neighbouring points of each spin and 
    * stores it
    ***********************************************************************/
    void neighbourTab();
    dim3 mGridSize; /**< three dimensional grid size */
    thrust::host_vector<T> mGrid; /**< vector containg spins in row major format on host */
    thrust::host_vector<NB> mTable; /**< vector containing neighbouring points in row major format on host */
    thrust::device_vector<T> d_mGrid; /**< vector containg spins in row major format on device */
    thrust::device_vector<NB> d_mTable; /**< vector containing neighbouring points in row major format on device */
    short int mDim; /**< dimension of the spins */
  };
}
