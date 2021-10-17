#ifndef HEISENBERG_H
#define Heisenberg_H

#include <model.h>

using std::vector;

class Heisenberg : public Model
{
 public:
  using SimulationStrategy = std::function<void(const Heisenberg&)>;
  
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
  void simulate() const override {mS(*this);}
  
  SimulationStrategy mS;
};

#endif
