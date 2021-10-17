#pragma once

#include <random>
 
#include <hdf5.h>

#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <iostream>
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
  template<typename T>
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

    NB getNeighbours(const dim3 index) const
    {
      size_t idx = index.x +
	           index.y*mGridSize.x +
                   index.z*mGridSize.x*mGridSize.y;
      return mNeighbourTable[idx];
    }

    /*****************************************************************//**
    * \brief sets grid size and resizes grid
    *********************************************************************/
    void setGridSize(const dim3 gridSize)
    {
      mGridSize = gridSize;

      mGridData.resize(gridSize.x*gridSize.y*gridSize.z);
    }    
    dim3 getGridSize() const {return mGridSize;}

    T getSpin(const dim3 index) const
    {
      size_t idx = index.x +
	           index.y*mGridSize.x +
                   index.z*mGridSize.x*mGridSize.y;
      return mGridData[idx];
    }
    void setSpin(const dim3 index, const T spin)
    {
      size_t idx = index.x +
	           index.y*mGridSize.x +
                   index.z*mGridSize.x*mGridSize.y;
      mGridData[idx] = spin;
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
    /*void upload(DeviceGrid &deviceGrid)
    {
      deviceGrid.setDim(mDim);
      deviceGrid.getGridSize(mGridSize);
      deviceGrid.setGridData(mGridData);
    }*/
    
    /************************************************************************************************//**
    * \brief copies grid to host
    ****************************************************************************************************/
    /*void download(DeviceGrid<T> &deviceGrid)
    {
      mDim      = deviceGrid.getDim();
      mGridSize = deviceGrid.setGridSize();
      mGridData = deviceGrid.getGridData();
    }*/
    
   private:
    void setGridData(const thrust::host_vector<T> &gridData) {mGridData = gridData;};
    thrust::host_vector<T> getGridData() const {return mGridData;};
    
    thrust::host_vector<T> mGridData{}; /**< vector containg spins in row major format on host */
    thrust::host_vector<NB> mNeighbourTable{}; /**< vector containing neighbouring points in row major format on device */

    dim3 mGridSize;
    short int mDim{3};
  };

  /*template class Grid<Spin>;
  template class Grid<float>;
  template class Grid<short int>;*/

  using  SintGrid = Grid<short int>;
  using  SpinGrid = Grid<Spin>;
  using FloatGrid = Grid<float>;

  using  SintGrids = thrust::host_vector<std::shared_ptr<Grid<short int>>>;
  using  SpinGrids = thrust::host_vector<std::shared_ptr<Grid<Spin>>>;
  using FloatGrids = thrust::host_vector<std::shared_ptr<Grid<float>>>;
}
