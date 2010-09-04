// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "../matrix.h"
#include "../solvers.h"

#ifdef COMMON_WITH_MUMPS
#include "mpi.h"
#include "dmumps_c.h"

// macro s.t. indices match Fortran documentation
#define ICNTL(I)                icntl[(I)-1]

#define JOB_INIT        -1
#define JOB_END         -2
#define JOB_SOLVE        6
#define USE_COMM_WORLD  -987654

bool CommonSolverMumps::_solve(Matrix *mat, double *res)
{
    // printf("MUMPS solver\n");

    CooMatrix *Acoo = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acoo = mcoo;
    else if (DenseMatrix *mden = dynamic_cast<DenseMatrix *>(mat))
        Acoo = new CooMatrix(mden);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acoo = new CooMatrix(mcsc);
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acoo = new CooMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    int nnz = Acoo->get_nnz();
    int size = Acoo->get_size();

    int argc = 1;
    char *name = "mumps";
    char **argv ;

    DMUMPS_STRUC_C id;

    int mpi_id;
    int ierr;
    argv = &name;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

    // Initialize a MUMPS instance. Use MPI_COMM_WORLD
    id.job = JOB_INIT;
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);

    // Define the problem on the host
    if (mpi_id == 0)
    {
        id.n = size;
        id.nz = nnz;

        id.irn = new int[nnz];
        id.jcn = new int[nnz];
        id.a = new double[nnz];
        Acoo->get_row_col_data(id.irn, id.jcn, id.a);

        for (int i = 0; i < nnz; i++)
        {
            id.irn[i] = id.irn[i]+1;
            id.jcn[i] = id.jcn[i]+1;
        }

        id.rhs = res;
    }
    // No outputs
    id.ICNTL(1) = -1;
    id.ICNTL(2) = -1;
    id.ICNTL(3) = -1;
    id.ICNTL(4) = 0;

    id.ICNTL(20) = 0; // centralized dense RHS
    id.ICNTL(21) = 0; // centralized dense solution

    // Call the MUMPS package.
    id.job = JOB_SOLVE;
    dmumps_c(&id);

    id.job = JOB_END;
    dmumps_c(&id); // Terminate instance

    ierr = MPI_Finalize();

    delete[] id.irn;
    delete[] id.jcn;
    delete[] id.a;

    if (!dynamic_cast<CooMatrix*>(mat))
        delete Acoo;

    return true;
}

bool CommonSolverMumps::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverMumps::solve(Matrix *mat, cplx *res) not implemented.");
}

#else

bool CommonSolverMumps::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverMumps::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverMumps::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverMumps::solve(Matrix *mat, cplx *res) not implemented.");
}

#endif
