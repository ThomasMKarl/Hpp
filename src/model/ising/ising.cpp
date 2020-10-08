#include "model/ising/ising.h"

double HB::Ising::Ising::calcEnergy() //store initial value for energy
{
    int spin_sum = 0;
    int neighbour_sum = 0;
    
    size_t n = this->mGrid->getGridSize().x*this->mGrid->getGridSize().y*this->mGrid->getGridSize().z;
    
    #pragma omp parallel for reduction(+:spin_sum,neighbour_sum)
    for(size_t i = 0; i < n; ++i)
    {
        dim3 idx;
        idx.x = i% this->mGrid->getGridSize().x;
	idx.y = i/ this->mGrid->getGridSize().x%this->mGrid->getGridSize().y;
        idx.z = i/(this->mGrid->getGridSize().x*this->mGrid->getGridSize().y);
		 
        short int mgrid = this->mGrid->getSpin(idx).x;
		 
        spin_sum += mgrid;
        neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(idx).right).x;
	if(this->mGrid->getDim() > 1)
	{
	    neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(idx).up).x;
	    if(this->mGrid->getDim() > 2)
	    {
  	        neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(idx).back).x;
	    }
	}
    }
    
    return -(mJ*neighbour_sum + mB.y*spin_sum);
}

double HB::Ising::Ising::calcEnergy(dim3 index)
{
    short int mgrid = this->mGrid->getSpin(index).x;
    int spin_sum = mgrid;
    int neighbour_sum = 0;
   
    neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).left).x; //mid point
    neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).right).x;
    if(this->mGrid->getDim() > 1)
    {
        neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).down).x;
        neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).up).x;
	if(this->mGrid->getDim() > 2)
	{
            neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).front).x;
	    neighbour_sum += mgrid*this->mGrid->getSpin(this->mGrid->getNeighbours(index).back).x;
	}
    }

    return 2*mJ*neighbour_sum + mB.y*spin_sum;
}

double HB::Ising::Ising::calcMagn() //store initial value for energy
{
    int spin_sum = 0;
    
    size_t n = this->mGrid->getGridSize().x*this->mGrid->getGridSize().y*this->mGrid->getGridSize().z;
    
    #pragma omp parallel for reduction(+:spin_sum)
    for(size_t i = 0; i < n; ++i)
    {
        dim3 idx;
        idx.x = i% this->mGrid->getGridSize().x;
	idx.y = i/ this->mGrid->getGridSize().x%this->mGrid->getGridSize().y;
        idx.z = i/(this->mGrid->getGridSize().x*this->mGrid->getGridSize().y);
        spin_sum += this->mGrid->getSpin(idx).x;
    }
    
    return spin_sum/double(n);
}
