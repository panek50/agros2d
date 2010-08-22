// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#include <coord_double.h>
#include <compcol_double.h>
#include <mvvd.h>
#include <mvmd.h>
#include <ilupre_double.h>
#include <icpre_double.h>
#include <diagpre_double.h>
#include <bicg.h>
#include <cg.h>
#include <cgs.h>
#include <bicgstab.h>
#include <cheby.h>
#include <gmres.h>
#include <ir.h>
#include <qmr.h>

bool CommonSolverSparseLib::_solve(Matrix *mat, double *res)
{
    // printf("SparseLib++ solver\n");

    CSCMatrix *Acsc = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acsc = new CSCMatrix(mcoo);
    else if (DenseMatrix *mden = dynamic_cast<DenseMatrix *>(mat))
        Acsc = new CSCMatrix(mden);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acsc = mcsc;
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acsc = new CSCMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    int nnz = Acsc->get_nnz();
    int size = Acsc->get_size();

    CompCol_Mat_double Acc = CompCol_Mat_double(size, size, nnz,
                                                Acsc->get_Ax(), Acsc->get_Ai(), Acsc->get_Ap());

    // rhs
    VECTOR_double rhs(res, size);

    // preconditioner
    CompCol_ILUPreconditioner_double pre(Acc);
    VECTOR_double xv = pre.solve(rhs);

    // method
    int result = -1;
    switch (method)
    {
    case CommonSolverSparseLibSolver_ConjugateGradient:
        {
            result = CG(Acc, xv, rhs, pre, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_ConjugateGradientSquared:
        {
            result = CGS(Acc, xv, rhs, pre, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_BiConjugateGradient:
        {
            result = BiCG(Acc, xv, rhs, pre, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_BiConjugateGradientStabilized:
        {
            result = BiCGSTAB(Acc, xv, rhs, pre, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_Chebyshev:
        {
            result = CHEBY(Acc, xv, rhs, pre, maxiter, tolerance, 0.01, 3.0);
            break;
        }
    case CommonSolverSparseLibSolver_QuasiMinimalResidual:
        {
            result = QMR(Acc, xv, rhs, pre, pre, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_GeneralizedMinimumResidual:
        {
            int restart = 32;
            MV_ColMat_double H(restart+1, restart, 0.0); // storage for upper Hessenberg H

            result = GMRES(Acc, xv, rhs, pre, H, restart, maxiter, tolerance);
            break;
        }
    case CommonSolverSparseLibSolver_RichardsonIterativeRefinement:
        {
            result = IR(Acc, xv, rhs, pre, maxiter, tolerance);
            break;
        }
    default:
        _error("SparseLib++ error. Method is not defined.");
    }

    printf("SparseLib++ solver: maxiter: %i, tol: %e\n", maxiter, tolerance);
    /*
    if (result == 0)
        ; // printf("SparseLib++ solver: maxiter: %i, tol: %e\n", maxiter, tolerance);
    else
        _error("SparseLib++ error.");
    */

    for (int i = 0 ; i < xv.size() ; i++)
        res[i] = xv(i);

    if (!dynamic_cast<CSCMatrix*>(mat))
        delete Acsc;

    return true;
}

bool CommonSolverSparseLib::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSparseLib::solve(Matrix *mat, cplx *res) not implemented.");
}
