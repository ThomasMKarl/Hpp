#include <cstdlib>
#include <vector>
#include <iostream>
#include <sys/times.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <simulation.h>
#include <model/ising/ising.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace HB;

int main(int argc, char **argv)
{
    double time_pre;
    #ifdef _OPENMP
      time_pre = omp_get_wtime();
    #endif
    ////////////////////////////////////////////////////////////////////
    
    if(argc < 10)
    {
      std::cout << "Usage: ising <temp> <J> <dim> <Bx> <By> <Bz> <Nx> <Ny> <Nz> <max_step> <hot_start (true/false)> <number of cores>\n";
      return EXIT_FAILURE;
    }
    uint arg = 0;
    float beta   = 1.0/atof(argv[++arg]);
    float J      = atof(argv[++arg]);
    short int dim = atoi(argv[++arg]);
    float3 B;
    B.x  = atof(argv[++arg]);
    B.y  = atof(argv[++arg]);
    B.z  = atof(argv[++arg]);
    dim3 N;
    N.x = atoi(argv[++arg]);
    N.y = atoi(argv[++arg]);
    N.z = atoi(argv[++arg]);
    size_t steps = atof(argv[++arg]);
    bool hot_start = argv[++arg];

    Grid<short int> grid(dim, N);
    grid.calcNeighbourTable();
    grid.hotStart();
    
    //Models models;
    //MetropolisSimulation s{grid,beta,1000};
    //models.push_back(std::make_unique<Ising>(J,B,{grid,beta,1000}));
    //models.push_back(std::make_unique<Heisenberg>(J,B,s));
    
    std::vector<double> energies;
    std::vector<double> magn;
    //simulate(models);
    
    /*removeCorr(energies);
    removeCorr(magn);
    
    float mean_E  =  mean(energies);
    float error_E = error(energies);
    float mean_M  =  mean(magn);
    float error_M = error(magn);
    std::cout << 1.0/beta << ' ' << mean_E << ' ' << error_E << ' '
    << mean_M << ' ' << error_M << '\n';*/

    ////////////////////////////////////////////////////////////////////

    double time_post;
    #ifdef _OPENMP
      time_post = omp_get_wtime();
    #endif
      
    std::cout << "#needed " << time_post - time_pre << " seconds." << std::endl;

    return EXIT_SUCCESS; 
}
