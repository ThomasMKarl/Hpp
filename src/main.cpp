#include <cstdlib>
#include <vector>
#include <iostream>
#include <sys/times.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ising/ising.h>

int main(int argc, char **argv)
{
    double time = omp_get_wtime();

    ////////////////////////////////////////////////////////////////////
    
    if(argc < 9)
    {
      std::cout << "Usage: ising <temp> <J> <B> <Nx> <Ny> <Nz> <max_step> <hot_start (true/false)> <number of cores>\n";
      return 1;
    }
    double beta = 1.0/atof(argv[1]);
    double J    = atof(argv[2]);
    double mag  = atof(argv[3]);
    unsigned int Nx = atoi(argv[4]);
    unsigned int Ny = atoi(argv[5]);
    unsigned int Nz = atoi(argv[6]);
    double steps = atof(argv[7]);
    bool start = argv[8];

    
    Ising x(beta, J, mag, Nx, Ny, Nz, start);
    std::vector<double> energies; std::vector<double> magn;
    x.production_run(steps, energies, magn);
    
    remove_corr(energies);
    remove_corr(magn);
    double mean_E  =  mean(energies);
    double error_E = error(energies);
    double mean_M  =  mean(magn);
    double error_M = error(magn);
    std::cout << 1.0/beta << ' ' << mean_E << ' ' << error_E << ' '
	                         << mean_M << ' ' << error_M << '\n';

    ////////////////////////////////////////////////////////////////////

    std::cout << "#needed " << omp_get_wtime() - time << " seconds." << std::endl;

    return 0; 
}
