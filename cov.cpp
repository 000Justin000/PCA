#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include "omp.h"

 

// F1
// read the vectors from file
std::vector< std::vector<double> > readFile(std::ifstream& fin)
{
    std::vector< std::vector<double> > mydata;
    std::string oneline;

    while( std::getline(fin, oneline) )
    {
        std::istringstream linestream(oneline);
        std::vector<double> current_vector;
        std::string element_str;
        while( linestream >> element_str )
        {
            double element_db = atof(element_str.c_str());
            current_vector.push_back(element_db);
        }

        mydata.push_back(current_vector);
    }

    return mydata;
}



// F2
// calculate the partial sum of vectors
std::vector<double> partialSum(std::vector< std::vector<double> >& mydata)
{
    std::vector<double> mysum( mydata[0].size() );
   
    for (int i = 0; i < mydata.size(); i++)
        for (int j = 0; j < mysum.size(); j++)
            mysum[j] += mydata[i][j];

    return mysum;
}


// F3
//void sleep(unsigned int mseconds)
//{
//    clock_t goal = mseconds + clock();
//    while (goal > clock());
//}



int main(int argc, char** argv)
{
    int iproc;
    int nproc;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    std::ostringstream oss;
    oss << std::setfill(' ') << std::setw(3) << iproc;
    
    std::ifstream fin;
    fin.open( ("pca" + oss.str()).c_str() );

    std::vector< std::vector<double> > mydata = readFile(fin); 

    int vector_dim = mydata[0].size();

    std::vector<double> mysum = partialSum(mydata);

    int my_vector_count = mydata.size();
    int biggest_vector_count = 0;
    int total_vector_count = 0;
    std::vector<int> all_vector_count(nproc);

    MPI_Allreduce(&my_vector_count, &biggest_vector_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&my_vector_count, &total_vector_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgather(&my_vector_count, 1, MPI_INT, &all_vector_count[0], 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<double> totalmean(mysum.size());
    MPI_Allreduce(&mysum[0], &totalmean[0], vector_dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < totalmean.size(); i++)
    {
        totalmean[i] /= total_vector_count;
    }

    if (iproc == 0)
    {
        std::ofstream fmout;
        fmout.open("mean");

        for (int i = 0; i < vector_dim; i++)
            fmout << std::setw(15) << totalmean[i];
    }

    for (int i = my_vector_count; i < biggest_vector_count; i++)
    {
        mydata.push_back(mydata[0]);
    }

    // substract the mean
    for (int i = 0; i < mydata.size(); i++)
    {
        for (int j = 0; j < vector_dim; j++)
            mydata[i][j] -= totalmean[j];
    }

    // determine prev_rank and next_rank in mpi_sendrecv
    int prev_rank = (iproc + nproc - 1) % nproc;
    int next_rank = (iproc + nproc + 1) % nproc;

    // determine how many rows current rank is in charge of
    int block_size = vector_dim / nproc;
    int my_covmat_row = block_size;
    if ( vector_dim % nproc != 0 )
    {
        block_size += 1;
        if ( iproc == nproc - 1 ) 
            my_covmat_row = vector_dim % block_size;
        else
            my_covmat_row = block_size;
    }

    std::cout << "iproc = " << iproc << "   " << "my_covmat_row" << my_covmat_row << std::endl;

    // declare the 2d array for the covmat current rank is in charge of
    std::vector< std::vector<double> > my_covmat;
    for (int i = 0; i < my_covmat_row; i++)
        my_covmat.push_back(std::vector<double> (vector_dim));
    
    for (int i = 0; i < nproc; i++)
    {
        int data_No = (iproc + nproc - i) % nproc;
        
        #pragma omp parallel for
        for (int k = 0; k < my_covmat_row; k++)
            for (int j = 0; j < all_vector_count[data_No]; j++)
                for (int l = 0; l < vector_dim; l++)
                    my_covmat[k][l] += mydata[j][iproc*block_size+k] * mydata[j][l];

        for (int j = 0; j < biggest_vector_count; j++)
        {
            std::vector<double> sendbuf = mydata[j];
            MPI_Sendrecv(&sendbuf[0], vector_dim, MPI_DOUBLE, next_rank, iproc, &mydata[j][0], vector_dim, MPI_DOUBLE, prev_rank, prev_rank, MPI_COMM_WORLD, &status);
        }

        if (iproc == 0)
            std::cout << "i=" << i << "finished" << std::endl;
    }

    std::ofstream fout;
    fout.open( ("covmat" + oss.str()).c_str() );
    for (int i = 0; i < my_covmat_row; i++)
    {
        for (int j = 0; j < vector_dim; j++)
            fout << std::setw(15) << my_covmat[i][j] / (total_vector_count - 1);
        fout << std::endl;
    }

    MPI_Finalize();
    return 0;
}

