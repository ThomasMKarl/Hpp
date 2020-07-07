#pragma once

#include <vector>
#include <cuda_runtime.h>

namespace HB
{
  typedef struct {
    size_t up;
    size_t right;
    size_t back;
    size_t down;
    size_t left;
    size_t front;
  }NB;
  template<typename T>
  class Spin
  {
   public:
    T x;
    T y;
    T z;
  };
  
  template<class T>
  class Grid
  {
   public:
    Grid();
    ~Grid();
    dim3 getGridSize() const {return mGridSize;};
    void setGridSize(dim3 gridSize) { mGridSize = gridSize;};
    template<T> Spin<T> getGridSize(dim3 index) const
    {
      size_t idx = 0;//!!!!
      return mGrid[idx];
    }
    NB getNeighbours(dim3 index) const
    {
      size_t idx = 0;//!!!!
      return mTable[idx];
    }
    /************************************************************************************************//**
    * \brief Initializes all spins represented by mgrid_phi and mgrid_theta randomly between 0 and 2*pi  
    ****************************************************************************************************/
    void hotStart();
    /*****************************************************************//**
    * \brief Initializes all spins represented by mphi and mtheta as up, 
    * i.e. mgrid_phi[...] = 0.5*pi and mgrid_theta[...] = 0                        
    *********************************************************************/
    void coldStart();
    /*******************************************************************//**
    * \brief Generates a table of the neighbouring points of each spin and 
    * stores it in member variable mmap.                                   
    ***********************************************************************/
    void saveGrid(std::string path);
    
   private:
    void neighbourTab();
    dim3 mGridSize;
    Spin<T> mSpin;
    std::vector<Spin<T>> mGrid;
    std::vector<NB> mTable;
  };
}
