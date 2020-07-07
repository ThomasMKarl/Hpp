#ifndef XY_H
#define XY_H

#include <model.h>


class XY : public Model
{
 public:
  /******************************************************************************************//**
  * \brief Constructs a XY Model.                                                   
  *                                                                                        
  * Initializes a neighbour table and reads the inverse beta, J, B, Nx, Ny and Nz.        
  * Also initializes the grid. If hc reads true, hot_start() is executed, otherwise cold_start().  
  **********************************************************************************************/
  XY(double beta, double J, double B, unsigned int Nx, unsigned int Ny, unsigned int Nz, bool hc = true) : Model(beta, J, B, Nx, Ny, Nz)
  {
    if(hc) hot_start();
    else  cold_start();
  };
  

  //functions
  double compute_energy() final;
  void  hot_start() final;
  void cold_start() final;

 protected:
  double delta_energy(unsigned int) final;
  void flip(unsigned int, gsl_rng*) final;
  
  
 private:
  vector<double> mgrid; ///< A vector containing double values. They represent the polar angulars between 0 and 2*pi of spins of unit length. The position in the grid (x,y,z) is (index%m_Nx, index/m_Nx%m_Ny, index/(m_Nx*m_Ny)) wher / denotes an integer division and % a modulo operation.
};

#endif
