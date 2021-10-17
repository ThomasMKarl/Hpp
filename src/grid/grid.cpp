#include "grid/grid.h"

template<class T>
void HB::Grid<T>::calcNeighbourTable()
{
    const size_t m_Nx = mGridSize.x;
    const size_t m_Ny = mGridSize.y;
    const size_t m_Nz = mGridSize.z;
    const size_t n = m_Nx*m_Ny*m_Nz;
    mNeighbourTable.resize(n);

    #pragma omp parallel
    {
    size_t i = 0;
    size_t j = 0;
    size_t k = 0;
      
    #pragma omp for
    for(size_t h = 0; h < n; ++h)
    {
            i = h%m_Nx;
	    
            if(i == 0)      mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].left = m_Ny*m_Nx*k + m_Nx*j + m_Nx - 1;
            else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].left = m_Ny*m_Nx*k + m_Nx*j + i    - 1;

            if(i == m_Nx-1) mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].right = m_Ny*m_Nx*k + m_Nx*j;
            else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].right = m_Ny*m_Nx*k + m_Nx*j + i    + 1;

	    if(mDim > 1)
	    {
	        j = h/m_Nx%m_Ny;
		
                if(j == 0)      mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].up = m_Ny*m_Nx*k + m_Nx*(m_Ny-1) + i;
                else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].up = m_Ny*m_Nx*k + m_Nx*(j-1)    + i;
	
                if(j == m_Ny-1) mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].down = m_Ny*m_Nx*k + i;
                else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].down = m_Ny*m_Nx*k + m_Nx*(j+1) + i;

	        if(mDim > 2)
	        {
		    k = h/(m_Nx*m_Ny);
		    
                    if(k == 0)      mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].front = m_Ny*m_Nx*(m_Nz-1) + m_Nx*j + i;
                    else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].front = m_Ny*m_Nx*(k-1)    + m_Nx*j + i;
	
                    if(k == m_Nz-1) mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].back = m_Nx*j + i;
                    else            mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].back = m_Ny*m_Nx*(k+1) + m_Nx*j + i;
	        }
		
		//std::cout << h << " " << i << " " << j << " " << k << "\n";
	}

	    /*std::cout << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].up    << " "
		  << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].right << " "
		  << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].back  << " "
		  << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].down  << " "
		  << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].left  << " "
		  << mNeighbourTable[m_Ny*m_Nx*k + m_Nx*j + i].front << "\n\n";*/
    }
    }
}

template<>
void HB::Grid<float>::hotStart() //all spins 1
{  
    #pragma omp parallel
    {
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<float> distribution(0,1);

    #pragma omp for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {
	mGridData[i]   = 2.0f*M_PI*distribution(generator);
    }
    }
}

template<>
void HB::Grid<short int>::hotStart() //all spins 1
{  
    #pragma omp parallel
    {
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<float> distribution(0,1);
    
    #pragma omp for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {	  
	if(distribution(generator) >= 0.5) mGridData[i] =  1;
	else                               mGridData[i] = -1;
    }
    }
}

template<>
void HB::Grid<HB::Spin>::hotStart() //all spins 1
{  
    #pragma omp parallel
    {
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_real_distribution<float> distribution(0,1);
    
    #pragma omp for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {	  
	mGridData[i].phi   = 2.0*M_PI*distribution(generator);
	mGridData[i].theta = 2.0*M_PI*distribution(generator);
    }
    }
}

template<>
void HB::Grid<float>::coldStart() //all spins 1
{ 
    #pragma omp parallel for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {
	mGridData[i] = 1.5f*M_PI;
    }   
}

template<>
void HB::Grid<short int>::coldStart() //all spins 1
{ 
    #pragma omp parallel for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {
	mGridData[i] = -1;
    }   
}

template<>
void HB::Grid<HB::Spin>::coldStart() //all spins 1
{ 
    #pragma omp parallel for
    for(size_t i = 0; i < mGridSize.x*mGridSize.y*mGridSize.z; ++i)
    {
	mGridData[i].phi   = 1.5f*M_PI;
	mGridData[i].theta = 0.0f;
    }   
}

template<>
herr_t HB::Grid<HB::Spin>::saveGrid(const std::string path)
{
    herr_t status;

    hsize_t dims[3];
    dims[0] = mGridSize.x;
    dims[1] = mGridSize.y;
    dims[2] = mGridSize.z;
    auto dataspace_id = H5Screate_simple(3, dims, NULL);

    hid_t file_id = H5Fcreate (path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

    hid_t dataset_id, spinid = H5Tcreate(H5T_COMPOUND, sizeof(HB::Spin));

    H5Tinsert(spinid, "phi",   HOFFSET(HB::Spin,phi), H5T_NATIVE_FLOAT);
    H5Tinsert(spinid, "theta", HOFFSET(HB::Spin,theta), H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate(file_id, "/dset", H5T_IEEE_F32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	

    for(size_t h = 0; h < dims[0]*dims[1]*dims[2]; ++h)
    {
	if(H5Dwrite(dataset_id, spinid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mGridData[h]) < 0)
	  return status;
    }

    if(H5Dclose(dataset_id)   < 0) return status;
    if(H5Sclose(dataspace_id) < 0) return status;
    if(H5Fclose(file_id)      < 0) return status;

    return status;
}

template<>
herr_t HB::Grid<short int>::saveGrid(const std::string path)
{
    herr_t status;

    hsize_t dims[3];
    dims[0] = mGridSize.x;
    dims[1] = mGridSize.y;
    dims[2] = mGridSize.z;
    auto dataspace_id = H5Screate_simple(3, dims, NULL);

    hid_t file_id = H5Fcreate (path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

    hid_t dataset_id, spinid = H5Tcreate(H5T_COMPOUND, sizeof(short int));
    
    H5Tinsert(spinid, "spin",   sizeof(short int), H5T_NATIVE_INT);
    dataset_id = H5Dcreate(file_id, "/dset", H5T_IEEE_F32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	

    for(size_t h = 0; h < dims[0]*dims[1]*dims[2]; ++h)
    {
	if(H5Dwrite(dataset_id, spinid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mGridData[h]) < 0) return status;
    }

    if(H5Dclose(dataset_id)   < 0) return status;
    if(H5Sclose(dataspace_id) < 0) return status;
    if(H5Fclose(file_id)      < 0) return status;

    return status;
}

template<>
herr_t HB::Grid<float>::saveGrid(const std::string path)
{
    herr_t status;

    hsize_t dims[3];
    dims[0] = mGridSize.x;
    dims[1] = mGridSize.y;
    dims[2] = mGridSize.z;
    auto dataspace_id = H5Screate_simple(3, dims, NULL);

    hid_t file_id = H5Fcreate (path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

    hid_t dataset_id, spinid = H5Tcreate(H5T_COMPOUND, sizeof(float));
    
    H5Tinsert(spinid, "spin",   sizeof(float), H5T_NATIVE_INT);
    dataset_id = H5Dcreate(file_id, "/dset", H5T_IEEE_F32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	

    for(size_t h = 0; h < dims[0]*dims[1]*dims[2]; ++h)
    {
	if(H5Dwrite(dataset_id, spinid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mGridData[h]) < 0) return status;
    }

    if(H5Dclose(dataset_id)   < 0) return status;
    if(H5Sclose(dataspace_id) < 0) return status;
    if(H5Fclose(file_id)      < 0) return status;

    return status;
}

template class HB::Grid<HB::Spin>;
template class HB::Grid<float>;
template class HB::Grid<short int>;
