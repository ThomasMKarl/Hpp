#include "model/ising/ising.h"

float HB::Ising::Ising::calcEnergy(const HB::Grid<short int> &grid) //store initial value for energy
{
    int spin_sum = 0;
    int neighbour_sum = 0;
    
    const size_t n = grid.getGridSize().x * grid.getGridSize().y * grid.getGridSize().z;
    
    #pragma omp parallel for reduction(+:spin_sum,neighbour_sum)
    for(size_t i = 0; i < n; ++i)
    {
        dim3 idx;
        idx.x = i% grid.getGridSize().x;
	idx.y = i/ grid.getGridSize().x%grid.getGridSize().y;
        idx.z = i/(grid.getGridSize().x*grid.getGridSize().y);
		 
        short int mgrid = grid.getSpin(idx);
		 
        spin_sum += mgrid;
        neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(idx).right);
	if(grid.getDim() > 1)
	{
	    neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(idx).up);
	    if(grid.getDim() > 2)
	    {
  	        neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(idx).back);
	    }
	}
    }
    
    return -(mJ*neighbour_sum + mB.y*spin_sum);
}

float HB::Ising::Ising::calcEnergy(const HB::Grid<short int> &grid, const dim3 index)
{
    const short int mgrid = grid.getSpin(index);
    int spin_sum = mgrid;
    int neighbour_sum = 0;
   
    neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).left); //mid point
    neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).right);
    if(grid.getDim() > 1)
    {
        neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).down);
        neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).up);
	if(grid.getDim() > 2)
	{
            neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).front);
	    neighbour_sum += mgrid*grid.getSpin(grid.getNeighbours(index).back);
	}
    }

    return 2*mJ*neighbour_sum + mB.y*spin_sum;
}

float HB::Ising::Ising::calcMagnetization(const HB::Grid<short int> &grid) //store initial value for energy
{
    int spin_sum = 0;
    
    const size_t n = grid.getGridSize().x * grid.getGridSize().y * grid.getGridSize().z;
    
    #pragma omp parallel for reduction(+:spin_sum)
    for(size_t i = 0; i < n; ++i)
    {
        dim3 idx;
        idx.x = i% grid.getGridSize().x;
	idx.y = i/ grid.getGridSize().x%grid.getGridSize().y;
        idx.z = i/(grid.getGridSize().x*grid.getGridSize().y);
	
        spin_sum += grid.getSpin(idx);
    }
    
    return spin_sum/float(n);
}

void HB::Ising::Ising::flip(const dim3 index, Grid<short int> &grid)
{
    grid.setSpin(index, -grid.getSpin(index));
}
