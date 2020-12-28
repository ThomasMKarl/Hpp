#pragma once

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <hdf5.h>
#include <random>

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
  typedef struct {
    float phi; /**< polar angle */
    float theta; /**< azimuthal angle */
  }Spin;
  
  //! A class storing spin values as a grid and certain meta information
  template<class T>
  class Grid
  {
   public:
    Grid() = default;

    /*****************************************************************//**
    * \brief constructs a grid and allocates memory
    *********************************************************************/
    Grid(const short int spatialDimensions, const dim3 gridSize) : mDim(spatialDimensions), mGridSize(gridSize)
    {
      mGridData.resize(gridSize.x*gridSize.y*gridSize.z);
    };

    /*******************************************************************//**
    * \brief generates a table of the neighbouring points of each spin and 
    * stores it
    ***********************************************************************/
    void calcNeighbourTable();

    /*****************************************************************//**
    * \brief sets grid size nad resizes grid
    *********************************************************************/
    void setGridSize(const dim3 gridSize)
    {
      mGridSize = gridSize;

      mGridData.resize(gridSize.x*gridSize.y*gridSize.z);
    }
    
    dim3 getGridSize() const
    {
      return mGridSize;
    }

    thrust::host_vector<T> getGridData() const
    {
      return mGridData;
    }
    
    void setGridData(const thrust::host_vector<T> &gridData)
    {
      mGridData = gridData;
    }

    thrust::host_vector<NB> getNeighbourTable() const
    {
      return mNeighbourTable;
    }
    
    void setNeighbourTable(const thrust::host_vector<NB> &neighbourTable)
    {
      mNeighbourTable = neighbourTable;
    }

    T getSpin(const dim3 index) const
    {
      size_t idx = index.x*mGridSize.x +
	           index.y*mGridSize.x*mGridSize.y +
                   index.z;
      return mGridData[idx];
    }

    void setSpin(const dim3 index, const T spin)
    {
      size_t idx = index.x*mGridSize.x +
	           index.y*mGridSize.x*mGridSize.y +
                   index.z;
      mGridData[idx] = spin;
    }
    
    NB getNeighbours(const dim3 index) const
    {
      size_t idx = index.x*mGridSize.x +
	           index.y*mGridSize.x*mGridSize.y +
                   index.z;
      return mNeighbourTable[idx];
    }
    
    short int getDim() const {return mDim;}

    void setDim(const short int dim) {mDim = dim;}

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
    herr_t saveGrid(const std::string path);
    
    /************************************************************************************************//**
    * \brief copies grid and table to device 
    ****************************************************************************************************/
    /*void upload(Grid<T> &deviceGrid)
    {
      deviceGrid.setGridData(mGridData);
      deviceGrid.setNeighbourTable(mNeighbourTable);
      }*/
    
    /************************************************************************************************//**
    * \brief copies grid to host
    ****************************************************************************************************/
    /*void download(Grid<T> &deviceGrid)
    {
      mGridData = deviceGrid.getGridData();
      }*/
    
   private:
    dim3 mGridSize; /**< three dimensional grid size */
    short int mDim; /**< dimension of the spins */
    thrust::host_vector<T> mGridData; /**< vector containg spins in row major format on host */
    thrust::host_vector<NB> mNeighbourTable; /**< vector containing neighbouring points in row major format on device */
  };
}
