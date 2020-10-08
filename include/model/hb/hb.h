#ifndef HEISENBERG_H
#define Heisenberg_H

#include <model.h>

using std::vector;

class Heisenberg : public Model
{
 public:
  /******************************************************************************************//**
  * \brief Constructs a Heisenberg Model.                                                   
  *                                                                                        
  * Initializes a neighbour table and reads the inverse beta, J, B, Nx, Ny and Nz.        
  * Also initializes the grid. If hc reads true, hot_start() is executed, otherwise cold_start().  
  **********************************************************************************************/
  Heisenberg(double beta, double J, double B, unsigned int Nx, unsigned int Ny, unsigned int Nz, bool hc = true) : Model(beta, J, B, Nx, Ny, Nz)
  { 
    if(hc) hot_start();
    else  cold_start();
  };

  ~Heisenberg() = default;


  //functions
  double compute_energy() final;
  void  hot_start() final;
  void cold_start() final;
  
  
 protected:
  double delta_energy(unsigned int) final;
  void flip(unsigned int, gsl_rng*) final;
  
 private:
  vector<double> mgrid_phi; ///< A vector containing double values. They represent the polar angulars between 0 and 2*pi of spins of unit length. The position in the grid (x,y,z) is (index%m_Nx, index/m_Nx%m_Ny, index/(m_Nx*m_Ny)) wher / denotes an integer division and % a modulo operation.
  vector<double> mgrid_theta; ///< A vector containing double values. They represent the azimutal angulars between 0 and pi of spins of unit length. The position in the grid (x,y,z) is (index%m_Nx, index/m_Nx%m_Ny, index/(m_Nx*m_Ny)) wher / denotes an integer division and % a modulo operation.
};

#endif
