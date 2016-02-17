#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "mpi.h"

int main(int argc, char** argv)
{
    int iproc;
    int nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(3) << iproc;

    std::cout << oss.str();

    return 0;

}
