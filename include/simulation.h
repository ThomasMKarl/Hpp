#pragma once

//#include <QtGui/QPainter>

#include <random>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include "model/ising/ising.h"

#include <libxml++/libxml++.h>

namespace HB
{ 
//__global__ void setupKernel(curandState *state, size_t seed, size_t n);

class SimulationSetup
{
 public:
  SimulationSetup() = default;

  void readXML(std::string path)
  {
    /*std::unique_ptr<xmlNode> help;
    std::unique_ptr<xmlChar> J, Bx, By, Bz;
    std::unique_ptr<xmlNode> currentXMLNode, nextXMLNode;
	
    XMLDocument = std::unique_ptr<xmlDoc>(xmlParseFile(path.c_str()));
    if(XMLDocument == NULL)
    {
      std::cerr << "Document not parsed successfully.\n";
      return;
    }

    currentXMLNode = std::unique_ptr<xmlNode>(xmlDocGetRootElement(XMLDocument.get()));
    if(currentXMLNode == NULL)
    {
      std::cerr <<"empty document.\n";
      return;
    }

    if(xmlStrcmp(currentXMLNode->name, (const xmlChar *) "MetropolisSimulation"))
    {
      currentXMLNode = std::unique_ptr<xmlNode>(currentXMLNode->xmlChildrenNode);
      while(currentXMLNode != NULL)
      {
	if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"models")))
        {
	  while (currentXMLNode != NULL)
	  {
	    if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"J")))
	    {
	      nextXMLNode = std::unique_ptr<xmlNode>(currentXMLNode->xmlChildrenNode);
	      J = xmlNodeListGetString(XMLDocument.get(), nextXMLNode.get(), 1);
	    }

	    if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"B")))
	    {
              help = std::unique_ptr<xmlNode>(currentXMLNode->xmlChildrenNode);
	      while (currentXMLNode != NULL)
	      {
	        if((!xmlStrcmp(help->name, (const xmlChar *)"x")))
		{
		  Bx = xmlNodeListGetString(XMLDocument.get(), help->xmlChildrenNode, 1);
 	        }
		if((!xmlStrcmp(help->name, (const xmlChar *)"y")))
		{
		  By = xmlNodeListGetString(XMLDocument.get(), help->xmlChildrenNode, 1);
 	        }
		if((!xmlStrcmp(help->name, (const xmlChar *)"z")))
		{
		  Bz = xmlNodeListGetString(XMLDocument.get(), help->xmlChildrenNode, 1);
 	        }

	        help = std::unique_ptr<xmlNode>(help->next);
	      }
	      
	      currentXMLNode = std::unique_ptr<xmlNode>(currentXMLNode->next);
            }
	  }
	}
	
        if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"sintGrids")))
        {
        }

	if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"floatGrids")))
        {
        }

	if((!xmlStrcmp(currentXMLNode->name, (const xmlChar *)"spinGrids")))
        {
        }
	
	currentXMLNode = std::unique_ptr<xmlNode>(currentXMLNode->next);
      }
    }*/
  }
  
  Models models{};
  SintGrids sintGrids{};
  FloatGrids floatGrids{};
  SpinGrids spinGrids{};
  size_t maxNumberOfSimulations{1};

  //std::unique_ptr<xmlDoc>  XMLDocument{};
};

class MPIsetup
{
 public:
  MPIsetup(dim3 N_) : N(N_)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &num );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ranksPerDimension{} = sqrtf(num);
    if(ranksPerDimension*ranksPerDimension != num)
      std::cout << "Warning: The number of ranks is not quadratic!\n";

    rightBorder.resize(N.y);
    leftBorder.resize(N.y);
    upBorder.resize(N.x);
    downBorder.resize(N.x);

    rankXCoordinate = rank/ranksPerDimension;
    rankYCoordinate = rank%ranksPerDimension;

    leftRank  = rankXCoordinate*ranksPerDimension +
      (rankYCoordinate+ranksPerDimension-1) %ranksPerDimension;
    rightRank = rankXCoordinate*ranksPerDimension +
      (rankYCoordinate+ranksPerDimension+1) %ranksPerDimension;
    upRank    = (rankXCoordinate+ranksPerDimension-1)%ranksPerDimension*
      ranksPerDimension + rankYCoordinate;
    downRank  = (rankXCoordinate+ranksPerDimension+1)%ranksPerDimension*
      ranksPerDimension + rankYCoordinate;
  }

  MPI_Status edgeCommunicationSendRecv(Grid<short int> &grid,
				       MPI_Communciator communicator = MPI_COMM_WORLD)
  {
    for(size_t i = 0; i < N.x; ++i)
    {
      //std::cout << i << " " << m_Ny/size*(m_Nx/size+1)+m_Ny/size+1+i << std::endl;
      downBorder[i]  = mgrid[N.x*(N.y+2)+1+i];
    }
    for(size_t j = 0; j < N.y; ++j)
    {
      //std::cout << j << " " << (j+2)*(m_Nx/size+1)+j << std::endl;
      rightBorder[j] = mgrid[(j+2)*(N.x+1)+j];
    }
	  
    //edge communication
    if(rankXCoordinate == rankYCoordinate)
    {
      MPI_Send(downBorder,  N.x, MPI_INT, downRank,  0, communicator);
      MPI_Recv(upBorder,    N.x, MPI_INT, upRank,    0, communicator, status); 
      //std::cout << rank << ":" << lown << std::endl;

      MPI_Send(rightBorder, N.y, MPI_INT, rightRank, 1, communicator);
      MPI_Recv(leftBorder,  N.y, MPI_INT, leftRank,  1, communicator, status);
      //std::cout << rank << ":" << rn << std::endl;
    }
    else
    {
      MPI_Recv(upBorder,    N.x, MPI_INT, upRank,    0, communicator, status);
      MPI_Send(downBorder,  N.x, MPI_INT, downRank,  0, communicator);	
      //std::cout << rank << ":" << lown << std::endl;
	      
      MPI_Recv(leftBorder,  N.y, MPI_INT, leftRank,  1, communicator, status);
      MPI_Send(rightBorder, N.y, MPI_INT, rightRank, 1, communicator);
      //std::cout << rank << ":" << ln << std::endl;
    }
      
    return status;
  }
  
  MPI_Status edgeCommunicationRecvSend(Grid<short int> &grid,
				       MPI_Communciator communicator = MPI_COMM_WORLD)
  {
    for(size_t i = 0; i < N.x; ++i)
      upBorder[i]   = mgrid[N.x+3+i];

    for(size_t j = 0; j < N.y; ++j)
      leftBorder[j] = mgrid[(N.x+2)*(j+1)+1];

    if(rankXCoordinate == rankYCoordinate)
    {
      MPI_Send(upBorder,    N.x, MPI_INT, upRank,    2, communicator);
      MPI_Recv(downBorder,  N.x, MPI_INT, downRank,  2, communicator, status);
	
      MPI_Send(leftBorder,  N.y, MPI_INT, leftRank,  3, communicator);
      MPI_Recv(rightBorder, N.y, MPI_INT, rightRank, 3, communicator, status);
    }
    else
    {
      MPI_Recv(downBorder,  N.x, MPI_INT, downRank,  2, MPI_COMM_WORLD, status);
      MPI_Send(upBorder,    N.x, MPI_INT, upRank,    2, MPI_COMM_WORLD);

      MPI_Recv(rightBorder, N.y, MPI_INT, rightRank, 3, MPI_COMM_WORLD, status);
      MPI_Send(leftBorder,  N.y, MPI_INT, leftRank,  3, MPI_COMM_WORLD);
    }

    return status;
  }
  
  void updateOverlapUpLeft(   Grid<short int> &grid)
  {
      for(size_t j = 0; j < N.x; ++j)
        mgrid[(N.x+2)*(j+1)]     = upBorder[j];
      for(size_t i = 0; i < N.y; ++i)
        mgrid[i+1]               = leftBorder[i];
  }

  void updateOverlapDownRight(Grid<short int> &grid)
  {
    for(size_t i = 0; i < N.x; ++i)
      mgrid[(N.x+1)*(N.y+2)+1+i] = downBorder[i];

    for(size_t j = 0; j < N.y; ++j)
      mgrid[(j+2)*(N.x+1)+j+1]   = rightBorder[j];
  }

  void computeSum(float partialEnergy,
		  MPI_Communciator communicator = MPI_COMM_WORLD)
  {
    float energy{0.0f};
    MPI_Reduce(&partialEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, communicator);
    return energy;
  }

  
  dim3 N{};
  size_t rank{};
  size_t num{};
  size_t ranksPerDimension{}
  std::vector<short int> rightBorder{}
  std::vector<short int> leftBorder{}
  std::vector<short int> upBorder{}
  std::vector<short int> downBorder{}

  short int rankXCoordinate{};
  short int rankYCoordinate{};

  size_t leftRank{};
  size_t rightRank{};
  size_t upRank{};
  size_t downRank{};
};
  
//////////////////////////////////////////////////////

class MetropolisSimulation
{
 public:
  MetropolisSimulation() = default;
  
  template<class M>
  void operator()(const M& model,
		  Grid<short int> &grid,
		  const float beta,
		  const size_t steps,
		  std::vector<float> &energies,
		  std::vector<float> &magnetization) const
  {
    simulate(model, grid, beta, steps, energies, magnetization);
  }
  
/********************************************************************//**
* \brief performs one Monte-Carlo step, updating the entire grid 
* associated with a certain model and calculates the energy difference
*
* \param Reference to the model
* \param beta Inverse Temperature
* \param rng Pointer to the state of the GSL RNG
* \return Energy difference
************************************************************************/
  template <class M>
  float metropolisSweep(const M &model, Grid<short int> &grid, const float beta) const
  {
    float new_E = 0.0;

    const dim3 N = grid.getGridSize();  
    const size_t V = N.x*N.y*N.z;
    	
    float dE;
    dim3 idx;
    
    #pragma omp parallel
    {
    float dE;
    dim3 idx;

    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<float> distribution(0,1);
    
    #pragma omp for reduction(+:new_E)
    for(size_t j = 0; j < V; ++j)
    {	
        idx.x = j% N.x;
	idx.y = j/ N.x%N.y;
        idx.z = j/(N.x*N.y);
	    
    	dE = model.calcEnergy(grid, idx);
	  
        if(dE < 0 || exp(-beta*dE) > distribution(generator))
        {
	    model.flip(idx, grid);
            new_E += dE;
        }
    }
    }
    
    return new_E;
  }

  float wolffSweep(const Ising &ising, Grid<short int> &grid, const float beta)
  {
    const dim3 N = grid->getGridSize();
    const size_t Volume = double(N.x*N.y*N.z);
    for(unsigned int i = 0; i < N.x; i++)
    {
      for(unsigned int j = 0; j < N.y; j++)
      {
        for(unsigned int k = 0; k < N.y; k++)
        {
	  short int oldSpin = grid.getSpin({N.x,N.y,N.z});
          if(oldSpin == -2) grid.setSpin({N.x,N.y,N.z}, -1);
          if(oldSpin ==  2) grid.setSpin({N.x,N.y,N.z},  1);
        }
      }
    }

    float probability = 1.0f-exp(-2.0f*beta);

    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<float> distribution(0,1);
  
    std::vector<dim3> stack{};
    stack.push_back({d2i(distribution(generator)%N.x),
                     d2i(distribution(generator)%N.y),
                     d2i(distribution(generator)%N.z)});
  
    dim3 seed{};
    short int numberOfNeighbours = 2*grid.getDim();
    for(size_t i = 0; i < stack.size(); ++i)
    {
      seed = stack[i];
      short int currentSpin = grid.getSpin(seed);
      HB::NB neighbours = grid.getNeighbours(seed);
    
      for(short int j = 0; j < numberOfNeighbours; ++j)
      {
        short int neighbourSpin = ???;
        dim3 neighbourIndex = ???;
        if((currentSpin == neighbourSpin)
	   && (std::find(stack.begin(), stack.end(), neighbourSpin) == stack.end())) //if spin is equal to neighbour and not part of stack...
        {
          if(probability > distribution(generator)) //...make probability check
            stack.push_back(neighbourIndex);
        }
      }
      if(currentSpin ==  2 || currentSpin ==  1) grid.setSpin(-2);
      if(currentSpin == -2 || currentSpin == -1) grid.setSpin( 2);
      //multiply with two to make cluster distinguishable from rest of the grid
    }
  }
  
/*************************************************************************************************************//**
* \brief performs a production-run.                                                             
*                                                                                               
* A production run consists of a number of MC sweeps producing correlated energies 
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
*****************************************************************************************************************/
  template <class M>
  void simulate(const M &model, Grid<short int> &grid, const float beta, const size_t steps, std::vector<float> &energies, std::vector<float> &magnetization) const
  {
    energies.resize(0);
    magnetization.resize(0);

    float energy = model.calcEnergy(grid);

    const dim3 N = grid.getGridSize();
    const double Volume = double(N.x*N.y*N.z);

    energies.push_back(energy/Volume);
    magnetization.push_back(model.calcMagnetization(grid)/Volume);
    
    for(size_t i = 0; i < steps; ++i)
    {
      energy += metropolisSweep(model, grid, beta);

      std::cout << i << " " << energy/Volume << std::endl;
      energies.push_back(energy/Volume);
      magnetization.push_back(model.calcMagnetization(grid)/Volume);
    }    
  }
};
  
class MetropolisSimulationQt
{
 public:
  MetropolisSimulationQt(const Grid<short int> &grid, const float beta) : mBeta(beta)
  {
    mGrid = std::make_shared<Grid<short int>>(grid);
  }
  
  template <class M>
  void operator()(const M &model) const
  {
    simulate(model);
  }
 private:
  float mBeta{1.0f};
  size_t mSteps{1000};
  std::shared_ptr<Grid<short int>> mGrid{};
  
/*************************************************************************************************************//**
* \brief performs a production-run.                                                             
*                                                                                               
* A production run consists of a number of MC sweeps producing correlated energies 
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
*****************************************************************************************************************/
  template <class M>
  void simulate(const M &model);

/********************************************************************//**
* \brief performs one Monte-Carlo step, updating the entire grid 
* associated with a certain model and calculates the energy difference
*
* \param Reference to the model
* \param beta Inverse Temperature
* \param rng Pointer to the state of the GSL RNG
* \return Energy difference
************************************************************************/
  template <class M>
  float metropolisSweep(const M &model);
};

class MetropolisSimulationMPI
{
 public:
  MetropolisSimulationMPI() = default;
  
  void operator()(const Ising& ising,
		  Grid<short int> &grid,
		  const float beta,
		  const size_t steps,
		  std::vector<float> &energies,
		  std::vector<float> &magnetization) const
  {
    simulate(ising, grid, beta, steps, energies, magnetization);
  }

  void simulate(const Ising& ising,
     	        Grid<short int> &grid,
	        const float beta,
	        const size_t steps,
	        std::vector<float> &energies,
	        std::vector<float> &magnetization)
  {
    std::unique_ptr<MPI_Status> status{};
    MPIsetup{grid.getGridSize()};

    energies.resize(0);
    magnetization.resize();
    
    float energy{0.0f};
    float partialEnergy = model.calcEnergy(grid);
    const double V = MPIsetup.ranksPerDimension*MPIsetup.ranksPerDimension * N.x*N.y;
     
    for(size_t k = 0; k < steps; ++k)
    {
      status = MPIsetup.edgeCommunicationSendRecv(grid);  
      MPIsetup.updateOverlapUpLeft();
      status = MPIsetup.edgeCommunicationRecvSend(grid);
      MPIsetup.updateOverlapDownRight();
    
      partialEnergy += ising.metropolisSweep(model, grid, beta);
      energy = MPIsetup.computeSum(partialEnergy);

      if(rank == 0)
      {
        energies.push_back(energy/Volume);
        magnetization.push_back(model.calcMagnetization(grid)/Volume);
        //std::cout << E/V << std::endl;
      }
    }
  }
};

class SimulatedAnnealing
{
 public:
  SimulatedAnnealing(const Grid<short int> &grid, const float beta, const size_t steps, const float3 temperatureSteps) : mBeta(beta), mSteps(steps), mTemperatureSteps(temperatureSteps)
  {
    mGrid = std::make_shared<Grid<short int>>(grid);
  }
  
  template<class M>
  void operator()(const M& model) const
  {
    simulate(model);
  }
 private:
  float mBeta{1.0f};
  size_t mSteps{1000};
  float3 mTemperatureSteps{1.0f,0.0f,0.1f};
  std::shared_ptr<Grid<short int>> mGrid{};
  std::vector<float> mEnergies{};
  std::vector<float> mMagnetization{};
  
/*************************************************************************************************************************//**
* \brief Performs a production-run via simulated annealing.                                     
*                                                                                              
* The run starts at Temperature T = T_begin and ends with T = T_end. In steps of T_steps       
* a production run computes the initial state for the next run. In the last   
* one energies and magnetizations are computed. Simulated annealing consists of                
* \f$n = \text{steps}\times \text{floor}\left(\frac{T{\_}\text{begin}-T{\_}\text{end}}{T{\_}\text{step}}\right)\f$ updates producing correlated energies and magnetizations.
*
* \param Reference to the model
* \param beta Inverse temperature
* \param V Volume
* \param steps Number of MC sweeps 
* \param energies Gets overridden with a standard vector of correlated energies (initial and one for each step)
* \param magnetizations Gets overridden with a standard vector of correlated magnetizations (initial and one for each step)
* \param T_begin Initial temperature value
* \param T_end Minimum temperature value 
* \param T_step Step size from one temperature to another
* \param cuda Set to true foe GPU acceleration
*****************************************************************************************************************************/
  template<class M>
  void simulate(const M &model);
};

///////////////////////////////////////////

/***********************************************************************************//**
* \brief calculates the mean value of a standard vector \f$\mu = \sum_i o_i\f$
*
* \param o Reference to the standard vector
* \param return Mean value
***************************************************************************************/
template <typename T> T mean(const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the covariance function up to a certain index of a standard vector
*
* \param t Index in vector
* \param o Reference to the standard vector
* \return Value of the covariance function
***************************************************************************************/
template <typename T> T covFunc(const size_t t, const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the integrated autocorrelation time of a standard vector
*
* \param o Reference to the standard vector
* \return Integrated autocorrelation time
***************************************************************************************/
template <typename T> T intAuto(const std::vector<T> &o);

/***********************************************************************************//**
* \brief calculates the standard error (standard deviation) of a standard vector
*
* \param o Reference to the standard vector
* \return Standard error
***************************************************************************************/
template <typename T> T error(const std::vector<T> &o);
  
/***********************************************************************************//**
* \brief calculates the blocking error of a standard vector
*
* \param block_number Number of blocks the vector is divided into
* \param o Reference to the standard vector
* \return Blocking error
***************************************************************************************/
template <typename T> T blockingError(const std::vector<T> &o, const size_t block_number);

/***********************************************************************************//**
* \brief calculates the bootstrap error of a standard vector
*
* \param o Reference to the standard vector
* \param sample_size Size of the bootstrap samples
* \parma sample_number Number of resamplings
* \param tau Integrated autocorrelation time
* \return Bootstrap error
***************************************************************************************/
template <typename T> T bootstrapError(const std::vector<T> &o, const size_t sample_size, const size_t sample_number, const T tau);

/***********************************************************************************//**
* \brief calculates the error based on error propagation of a standard vector
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time
* \return Error propagation
***************************************************************************************/
template <typename T> T errorProp(const std::vector<T> &o, const T tau);

/**************************************************************************//**
* \brief calculates the statistical error on energies using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time 
* \return Error propagation                                    
******************************************************************************/
template <typename T> T errorPropEnergy(const std::vector<T> &energies, const T beta, const T tau, const size_t V);

/***************************************************************************************//**
* \brief calculates the statistical error on magnetizations (magn) using error propagation
*
* \param o Reference to the standard vector
* \param tau Integrated autocorrelation time   
* \return Error propagation                                                   
*******************************************************************************************/
template <typename T> T errorPropMagnetization(const std::vector<T> &magn, const T beta, const T tau, const size_t V);
  
/***************************************************************************************//**
* \brief rounds a floating point value to an integer
*
* \param d Value to be rounded
* \return Resulting integer
*******************************************************************************************/
template <typename T> size_t d2i(const T d);

/***************************************************************************************//**
* \brief removes correlation from a standard vector
*
* The function calculates the integrated autocorrelation time \f$\tau\f$. The thermalisation 
* of \f$20\tau\f$ is removed. The integrated autocorrelation time is calculated again in order 
* to remove values from the input vector until the remaining ones are uncorrelated. The 
* The function overrides the input vector.
*
* \param o Standard vector of correlated values
*******************************************************************************************/
template <typename T> void removeCorr(const std::vector<T> &o);

/****************************************************************************//**
* \brief computes the specific heat per Volume as secondary quantity for a 
* given vector of energies
*
* \param Reference to a standard vector of energies
* \param beta Inverse temperature
* \param V Volume
* \return Specific heat
********************************************************************************/
template <typename T> T specificHeat(const std::vector<T> &energies, const T beta, const size_t V);

/****************************************************************************//**
* \brief computes the magnetic suszeptibility per Volume as secondary quantity 
* for a given vector of magnetizations
*
* \param Reference to a standard vector of magnetizations
* \param beta Inverse temperature
* \param V Volume
* \return Magnetic suzeptibility
********************************************************************************/
template <typename T> T magnSusz(const std::vector<T> &magn, const T beta, const size_t V);
};
