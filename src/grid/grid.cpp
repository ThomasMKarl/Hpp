#include "grid/grid.h"

template<class T>
void HB::Grid<T>::hotStart() //all spins 1
{  
    #pragma omp parallel
    {
	gsl_rng *rng;
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        long seed = time(NULL);
        gsl_rng_set(rng,seed);

        #pragma omp for
        for(size_t i = 0; i < this.mGridSize.x*this.mGridSize.y*this.mGridSize.z; ++i)
	{
	    this.mGrid[i].x   = 2.0*M_PI*gsl_rng_uniform(rng);
	    this.mGrid[i].y = 2.0*M_PI*gsl_rng_uniform(rng);
	}
    }
}

template<class T>
void HB::Grid<T>::coldStart() //all spins 1
{ 
    #pragma omp parallel for
    for(size_t i = 0; i < this.mGridSize.x*this.mGridSize.y*this.mGridSize.z; ++i)
    {
	this.mGrid[i].x   = 0.0;
	this.mGrid[i].y = 0.5*M_PI;
    }   
}

template<class T>
void HB::Grid<T>::neighbourTab()
{
    size_t m_Nx = this.mGridSize.x;
    size_t m_Ny = this.mGridSize.y;
    size_t m_Nz = this.mGridSize.z;
    size_t n = m_Nx*m_Ny*m_Nz;
    this.mTable.resize(n);
        
    #pragma omp parallel
    {	
        size_t i, j, k;
	
	#pragma omp for
        for(size_t h = 0; h < n; ++h)
        {
            i = h%m_Nx;
	    
            if(i == 0)      mTable[m_Ny*m_Nx*k + m_Nx*j + i].left = m_Ny*m_Nx*k + m_Nx*j + m_Nx - 1;
            else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].left = m_Ny*m_Nx*k + m_Nx*j + i    - 1;

            if(i == m_Nx-1) mTable[m_Ny*m_Nx*k + m_Nx*j + i].right = m_Ny*m_Nx*k + m_Nx*j;
            else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].right = m_Ny*m_Nx*k + m_Nx*j + i    + 1;

	    if(mDim > 1)
	    {
	        j = h/m_Nx%m_Ny;
		
                if(j == 0)      mTable[m_Ny*m_Nx*k + m_Nx*j + i].up = m_Ny*m_Nx*k + m_Nx*(m_Ny-1) + i;
                else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].up = m_Ny*m_Nx*k + m_Nx*(j-1)    + i;
	
                if(j == m_Ny-1) mTable[m_Ny*m_Nx*k + m_Nx*j + i].down = m_Ny*m_Nx*k + i;
                else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].down = m_Ny*m_Nx*k + m_Nx*(j+1) + i;

	        if(mDim > 2)
	        {
		    k = h/(m_Nx*m_Ny);
		    
                    if(k == 0)      mTable[m_Ny*m_Nx*k + m_Nx*j + i].front = m_Ny*m_Nx*(m_Nz-1) + m_Nx*j + i;
                    else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].front = m_Ny*m_Nx*(k-1)    + m_Nx*j + i;
	
                    if(k == m_Nz-1) mTable[m_Ny*m_Nx*k + m_Nx*j + i].back = m_Nx*j + i;
                    else            mTable[m_Ny*m_Nx*k + m_Nx*j + i].back = m_Ny*m_Nx*(k+1) + m_Nx*j + i;
	        }
		
		//std::cout << h << " " << i << " " << j << " " << k << "\n";
	    }

	//std::cout << mmap[0][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[1][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[2][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[3][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[4][m_Ny*m_Nx*k + m_Nx*j + i]<<" "<<mmap[5][m_Ny*m_Nx*k + m_Nx*j + i]<<"\n\n";
        }
    }
}

template<class T>
herr_t HB::Grid<T>::saveGrid(std::string path)
{
    herr_t status;

    size_t dims[3];
    dims[0] = this.mGridSize.x;
    dims[1] = this.mGridSize.y;
    dims[2] = this.mGridSize.z;
    auto dataspace_id = H5Screate_simple(3, dims, NULL);

    hid_t file_id = H5Fcreate (path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

    hid_t dataset_id, spinid = H5Tcreate(H5T_COMPOUND, sizeof(T));
    if(sizeof(T) == 8)
    {
        H5Tinsert(spinid, "spin",   HOFFSET(T,x), H5T_NATIVE_SHORT);
	dataset_id = H5Dcreate(file_id, "/dset", H5T_STD_I16BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        H5Tinsert(spinid, "phi",   HOFFSET(T,x), H5T_NATIVE_FLOAT);
        H5Tinsert(spinid, "theta", HOFFSET(T,y), H5T_NATIVE_FLOAT);
	dataset_id = H5Dcreate(file_id, "/dset", H5T_IEEE_F32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	
    }

    size_t i, j, k;
    for(size_t h = 0; h < dims[0]*dims[1]*dims[2]; ++h)
    {
        i = h%dims[0];
        j = h/dims[0]%dims[1];
        k = h/(dims[0]*dims[1]);
	if(H5Dwrite(dataset_id, spinid, H5S_ALL, H5S_ALL, H5P_DEFAULT, this.mGrid[h]) < 0) return status;
    }

    if(H5Dclose(dataset_id)   < 0) return status;
    if(H5Sclose(dataspace_id) < 0) return status;
    if(H5Fclose(file_id)      < 0) return status;

    return status;
}
