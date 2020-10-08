#include <cstdlib>
#include <vector>
#include <iostream>
#include <sys/times.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm.h>
#include <model/ising/ising.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace HB;

int main(int argc, char **argv)
{
    double time = omp_get_wtime();

    ////////////////////////////////////////////////////////////////////
    
    if(argc < 10)
    {
      std::cout << "Usage: ising <temp> <J> <dim> <Bx> <By> <Bz> <Nx> <Ny> <Nz> <max_step> <hot_start (true/false)> <number of cores>\n";
      return 1;
    }
    uint arg = 0;
    double beta   = 1.0/atof(argv[++arg]);
    double J      = atof(argv[++arg]);
    short int dim = atoi(argv[++arg]);
    float3 B;
    B.x  = atof(argv[++arg]);
    B.y  = atof(argv[++arg]);
    B.z  = atof(argv[++arg]);
    dim3 N;
    N.x = atoi(argv[++arg]);
    N.y = atoi(argv[++arg]);
    N.z = atoi(argv[++arg]);
    double steps = atof(argv[++arg]);
    bool start = argv[++arg];

    Grid<Spin<short int>> grid(dim, N);
    grid.hotStart();
    //Ising::Ising I(J, B, &grid, start);
    std::vector<double> energies; std::vector<double> magn;
    //productionRun(I, N.x*N.y*N.z, steps, energies, false);
    
    //removeCorr(energies);
    //removeCcorr(magn);
    /*
    double mean_E  =  mean(energies);
    double error_E = error(energies);
    double mean_M  =  mean(magn);
    double error_M = error(magn);
    std::cout << 1.0/beta << ' ' << mean_E << ' ' << error_E << ' '
                                 << mean_M << ' ' << error_M << '\n';*/

    ////////////////////////////////////////////////////////////////////

    std::cout << "#needed " << omp_get_wtime() - time << " seconds." << std::endl;

    return 0; 
}
