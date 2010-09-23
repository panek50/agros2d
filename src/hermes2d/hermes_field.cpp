// This file is part of Agros2D.
//
// Agros2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Agros2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Agros2D.  If not, see <http://www.gnu.org/licenses/>.
//
// hp-FEM group (http://hpfem.org/)
// University of Nevada, Reno (UNR) and University of West Bohemia, Pilsen
// Email: agros2d@googlegroups.com, home page: http://hpfem.org/agros2d/

#include "hermes_field.h"

#include "hermes_general.h"
#include "hermes_electrostatic.h"
#include "hermes_magnetic.h"
#include "hermes_heat.h"
#include "hermes_current.h"
#include "hermes_elasticity.h"
#include "hermes_flow.h"

#include "scene.h"
#include "h2d_reader.h"

// #include <InpMtx.h>

bool isPlanar;
AnalysisType analysisType;
double frequency;
double actualTime;
double timeStep;

HermesField *hermesFieldFactory(PhysicField physicField)
{
    switch (physicField)
    {
    case PhysicField_General:
        return new HermesGeneral();
    case PhysicField_Electrostatic:
        return new HermesElectrostatic();
    case PhysicField_Magnetic:
        return new HermesMagnetic();
    case PhysicField_Heat:
        return new HermesHeat();
    case PhysicField_Current:
        return new HermesCurrent();
    case PhysicField_Elasticity:
        return new HermesElasticity();
    case PhysicField_Flow:
        return new HermesFlow();
    default:
        std::cerr << "Physical field '" + QString::number(physicField).toStdString() + "' is not implemented. hermesObjectFactory()" << endl;
        throw;
        break;
    }
}

Mesh *readMeshFromFile(const QString &fileName)
{
    // save locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    // load the mesh file
    Mesh *mesh = new Mesh();
    H2DReader meshloader;
    meshloader.load(fileName.toStdString().c_str(), mesh);

    // set system locale
    setlocale(LC_NUMERIC, plocale);

    return mesh;
}

void writeMeshFromFile(const QString &fileName, Mesh *mesh)
{
    // save locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    H2DReader meshloader;
    meshloader.save(fileName.toStdString().c_str(), mesh);

    // set system locale
    setlocale(LC_NUMERIC, plocale);
}

SolutionArray *solutionArray(Solution *sln, Space *space = NULL, double adaptiveError = 0.0, double adaptiveSteps = 0.0, double time = 0.0)
{
    SolutionArray *solution = new SolutionArray();
    solution->order = new Orderizer();
    if (space) solution->order->process_solution(space);
    solution->sln = new Solution();
    if (sln) solution->sln->copy(sln);
    solution->adaptiveError = adaptiveError;
    solution->adaptiveSteps = adaptiveSteps;
    solution->time = time;

    return solution;
}

    QList<SolutionArray *> *solveSolutioArray(ProgressItemSolve *progressItemSolve,
                                          void (*cbSpace)(Tuple<Space *>),
                                          void (*cbWeakForm)(WeakForm *, Tuple<Solution *>))
{
    int polynomialOrder = Util::scene()->problemInfo()->polynomialOrder;
    AdaptivityType adaptivityType = Util::scene()->problemInfo()->adaptivityType;
    int adaptivitySteps = Util::scene()->problemInfo()->adaptivitySteps;
    double adaptivityTolerance = Util::scene()->problemInfo()->adaptivityTolerance;
    int numberOfSolution = Util::scene()->problemInfo()->hermes()->numberOfSolution();
    double timeTotal = Util::scene()->problemInfo()->timeTotal.number;
    double initialCondition = Util::scene()->problemInfo()->initialCondition.number;

    timeStep = Util::scene()->problemInfo()->timeStep.number;
    isPlanar = (Util::scene()->problemInfo()->problemType == ProblemType_Planar);
    analysisType = Util::scene()->problemInfo()->analysisType;
    frequency = Util::scene()->problemInfo()->frequency;

    Linearity linearity = Util::scene()->problemInfo()->linearity;
    double linearityNewtonTolerance = Util::scene()->problemInfo()->linearityNewtonTolerance;
    int linearityNewtonMaxSteps = Util::scene()->problemInfo()->linearityNewtonMaxSteps;

    // solution agros array
    QList<SolutionArray *> *solutionArrayList = new QList<SolutionArray *>();

    // load the mesh file
    Mesh *mesh = readMeshFromFile(tempProblemFileName() + ".mesh");
    // refine mesh
    for (int i = 0; i < Util::scene()->problemInfo()->numberOfRefinements; i++)
        mesh->refine_all_elements(0);

    // initialize the shapeset
    H1Shapeset shapeset;

    // create an H1 space
    Tuple<Space *> space;
    // create hermes solution array
    Tuple<Solution *> solution;
    // create reference solution
    Tuple<Solution *> solutionReference;

    for (int i = 0; i < numberOfSolution; i++)
    {
        // space
        space.push_back(new H1Space(mesh, NULL, NULL, 1, &shapeset));
        // set order by element
        for (int j = 0; j < Util::scene()->labels.count(); j++)
            space.at(i)->set_uniform_order(Util::scene()->labels[j]->polynomialOrder > 0 ? Util::scene()->labels[j]->polynomialOrder : polynomialOrder, j);

        // solution agros array
        solution.push_back(new Solution());

        // reference solution
        if ((adaptivityType != AdaptivityType_None))
            solutionReference.push_back(new Solution());
    }

    // callback space
    cbSpace(space);

    int ndof = get_num_dofs(space);
    if (analysisType == AnalysisType_Transient)
    {
        for (int i = 0; i < numberOfSolution; i++)
        {
            // constant initial solution
            solution.at(i)->set_const(mesh, initialCondition);
            solutionArrayList->append(solutionArray(solution.at(i)));
        }
    }

    // initialize the weak formulation
    WeakForm wf(numberOfSolution);
    // callback weakform
    cbWeakForm(&wf, solution);

    // prepare selector
    RefinementSelectors::Selector *selector = NULL;
    switch (adaptivityType)
    {
    case AdaptivityType_H:
        selector = new RefinementSelectors::HOnlySelector();
        break;
    case AdaptivityType_P:
        selector = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_P_ANISO,
                                                                Util::config()->convExp,
                                                                H2DRS_DEFAULT_ORDER);
        break;
    case AdaptivityType_HP:
        selector = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_HP_ANISO,
                                                                Util::config()->convExp,
                                                                H2DRS_DEFAULT_ORDER);
        break;
    }


    // Set up the solver, matrix, and rhs according to the solver selection.
    MatrixSolverType matrix_solver = SOLVER_UMFPACK;

    // assemble the stiffness matrix and solve the system
    double error;

    // set actual time
    actualTime = 0;

    // error marker
    bool isError = false;

    // solution
    int maxAdaptivitySteps = (adaptivityType == AdaptivityType_None) ? 1 : adaptivitySteps;
    int actualAdaptivitySteps = -1;
    for (int i = 0; i<maxAdaptivitySteps; i++)
    {
        // initialize the FE problem
        FeProblem fep(&wf, space, (linearity == Linearity_Linear));

        // initialize matrix, vector and solver
        SparseMatrix *matrix = create_matrix(matrix_solver);
        Vector *rhs = create_vector(matrix_solver);
        Solver *solver = create_solver(matrix_solver, matrix, rhs);

        // assemble stiffness matrix and rhs.
        if (linearity == Linearity_Linear)
        {
            QTime time;
            time.start();
            fep.assemble(matrix, rhs, false);
            qDebug() << "solveSolutioArray: FeProblem::assemble: " << milisecondsToTime(time.elapsed()).toString("mm:ss.zzz");

            if (fep.get_num_dofs() == 0)
            {
                progressItemSolve->emitMessage(QObject::tr("Solver: DOF is zero"), true);
                isError = true;

                delete rhs;
                delete matrix;
                delete solver;

                break;
            }

            // solve the matrix problem.
            time.start();
            if (!solver->solve())
            {
                progressItemSolve->emitMessage(QObject::tr("Matrix solver failed."), true);
                isError = true;

                delete rhs;
                delete matrix;
                delete solver;

                break;
            }
            qDebug() << "solveSolutioArray: CommonSolver::solve: " << milisecondsToTime(time.elapsed()).toString("mm:ss.zzz");
        }
        else
        {
            /*
            // Newton loop
            int it = 0;

            // initial solution for Newton's method
            Solution *initSolution = new Solution();
            initSolution->set_const(mesh, 0.0);

            Vector *coeff = new AVector(ndof);

            // Project to obtain the initial coefficient vector for the Newton's method.
            // The empty solution Tuple means that we do not want the resulting Solution, just the vector
            project_global(space, H2D_H1_NORM, initSolution, Tuple<Solution *>(), coeff);

            while (1)
            {
                dp->assemble(coeff, mat, NULL, rhs, false);

                // multiply the residual vector with -1 since the matrix
                // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
                for (int i = 0; i < dp->get_num_dofs(); i++)
                    rhs->set(i, -rhs->get(i));

                // calculate the l2-norm of residual vector.
                double l2Norm = 0;
                for (int i = 0; i < rhs->get_size(); i++)
                    l2Norm += rhs->get(i)*rhs->get(i);
                l2Norm = sqrt(l2Norm);

                // emit signal
                progressItemSolve->emitMessage(QObject::tr("Newton iteration: %1, L2 norm: %2").
                                               arg(it+1).
                                               arg(l2Norm, 0, 'e', 5), false, 1);

                // if l2 norm of the residual vector is in tolerance, quit.
                if (l2Norm < linearityNewtonTolerance || it > linearityNewtonMaxSteps)
                    break;

                // solve the matrix problem.
                if (!solver->solve(mat, rhs))
                {
                    progressItemSolve->emitMessage(QObject::tr("Matrix solver failed."), true);
                    isError = true;
                    delete solver;
                    break;
                }

                // Add \delta Y^{n+1} to Y^n.
                for (int i = 0; i < ndof; i++)
                    coeff->add(i, rhs->get(i));

                it++;
            }

            // replace rhs
            for (int i = 0; i < ndof; i++)
                rhs->set(i, coeff->get(i));

            delete coeff;
            delete initSolution;
            */
        }

        // convert coefficient vector into a solution.
        vector_to_solutions(solver->get_solution(), space, solution);

        // calculate errors and adapt the solution
        if (adaptivityType != AdaptivityType_None)
        {
            // Construct globally refined reference mesh and setup reference space.
            Tuple<Space *> *spaceRef = construct_refined_spaces(space);

            // initialize the FE problem
            FeProblem fepRef(&wf, *spaceRef, (linearity == Linearity_Linear));

            // assemble ref stiffness matrix and rhs.
            fepRef.assemble(matrix, rhs, false);

            // solve the matrix problem.
            if (!solver->solve())
            {
                progressItemSolve->emitMessage(QObject::tr("Matrix solver for reference solution failed."), true);
                isError = true;
                delete spaceRef;
                break;
            }

            // convert coefficient vector into a solution.
            vector_to_solutions(solver->get_solution(), *spaceRef, solutionReference);

            // project the reference solution on the coarse mesh.
            project_global(space, H2D_H1_NORM, solutionReference, solution);

            // adaptivity
            Adapt hp(space, H2D_H1_NORM);
            hp.set_solutions(solution, solutionReference);
            error = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

            // emit signal
            progressItemSolve->emitMessage(QObject::tr("Relative error: %1 %").
                                           arg(error, 0, 'f', 5), false, 1);
            // add error to the list
            progressItemSolve->addAdaptivityError(error, get_num_dofs(space));

            // delete ref solver
            delete spaceRef;

            if (progressItemSolve->isCanceled())
            {
                isError = true;
                break;
            }

            if (error < adaptivityTolerance || get_num_dofs(space) >= NDOF_STOP)
                break;

            if (i != maxAdaptivitySteps-1) hp.adapt(selector,
                                                    Util::config()->threshold,
                                                    Util::config()->strategy,
                                                    Util::config()->meshRegularity);
            actualAdaptivitySteps = i+1;
        }

        delete rhs;
        delete matrix;
        delete solver;
    }

    // delete selector
    if (selector) delete selector;

    // timesteps
    if (!isError)
    {
        int timesteps = (analysisType == AnalysisType_Transient) ? floor(timeTotal/timeStep) : 1;
        for (int n = 0; n < timesteps; n++)
        {
            // set actual time
            actualTime = (n+1)*timeStep;

            // initialize the FE problem
            FeProblem fep(&wf, space, (linearity == Linearity_Linear));

            // initialize matrix, vector and solver
            SparseMatrix *matrix = create_matrix(matrix_solver);
            Vector *rhs = create_vector(matrix_solver);
            Solver *solver = create_solver(matrix_solver, matrix, rhs);

            if (timesteps > 1)
            {
                // transient - assemble stiffness matrix and rhs.
                QTime time;
                time.start();
                fep.assemble(matrix, rhs, false);
                qDebug() << "solveSolutioArray: FeProblem::assemble"; time.elapsed();

                if (fep.get_num_dofs() == 0)
                {
                    progressItemSolve->emitMessage(QObject::tr("Number of DOFs is zero"), true);
                    isError = true;
                    break;
                }

                // solve the matrix problem.
                if (!solver->solve())
                {
                    progressItemSolve->emitMessage(QObject::tr("Matrix solver failed."), true);
                    isError = true;
                    break;
                }

                // convert coefficient vector into a Solution.
                vector_to_solutions(solver->get_solution(), space, solution);
            }
            else if (n > 0)
            {
                fep.assemble(matrix, rhs, true);
            }

            // output
            for (int i = 0; i < numberOfSolution; i++)
            {
                solutionArrayList->append(solutionArray(solution.at(i), space.at(i), error, actualAdaptivitySteps, (n+1)*timeStep));
            }

            if (analysisType == AnalysisType_Transient)
                progressItemSolve->emitMessage(QObject::tr("Time step: %1/%2").
                                               arg(n+1).
                                               arg(timesteps), false, n+2);

            delete rhs;
            delete matrix;
            delete solver;

            if (progressItemSolve->isCanceled())
            {
                isError = true;
                break;
            }
        }
    }

    // delete mesh
    delete mesh;

    // delete space
    for (int i = 0; i < space.size(); i++)
        delete space.at(i);
    space.clear();

    // delete last solution
    for (int i = 0; i < solution.size(); i++)
        delete solution.at(i);
    solution.clear();

    // delete reference solution
    for (int i = 0; i < solutionReference.size(); i++)
        delete solutionReference.at(i);
    solutionReference.clear();

    if (isError)
    {
        for (int i = 0; i < solutionArrayList->count(); i++)
            delete solutionArrayList->at(i);
        solutionArrayList->clear();
    }

    return solutionArrayList;
}

// *********************************************************************************************************************************************

ViewScalarFilter::ViewScalarFilter(Tuple<MeshFunction *> sln, PhysicFieldVariable physicFieldVariable, PhysicFieldVariableComp physicFieldVariableComp)
    : Filter(sln)
{
    m_physicFieldVariable = physicFieldVariable;
    m_physicFieldVariableComp = physicFieldVariableComp;
}

double ViewScalarFilter::get_pt_value(double x, double y, int item)
{
    error("Not implemented");
}

void ViewScalarFilter::precalculate(int order, int mask)
{
    Quad2D* quad = quads[cur_quad];
    int np = quad->get_num_points(order);
    node = new_node(H2D_FN_DEFAULT, np);

    if (sln[0])
    {
        sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[0]->get_dx_dy_values(dudx1, dudy1);
        value1 = sln[0]->get_fn_values();
    }

    if (num >= 2 && sln[1])
    {
        sln[1]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[1]->get_dx_dy_values(dudx2, dudy2);
        value2 = sln[1]->get_fn_values();
    }

    if (num >= 3 && sln[2])
    {
        sln[2]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[2]->get_dx_dy_values(dudx3, dudy3);
        value3 = sln[2]->get_fn_values();
    }

    update_refmap();

    x = refmap->get_phys_x(order);
    y = refmap->get_phys_y(order);
    Element *e = refmap->get_active_element();

    labelMarker = Util::scene()->labels[e->marker]->marker;

    for (int i = 0; i < np; i++)
    {
        calculateVariable(i);
    }

    replace_cur_node(node);
}
