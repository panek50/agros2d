#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <QObject>
#include <QtPlugin>

#include "util.h"

#include "meshgenerator.h"
#include "field.h"
#include "problem.h"
#include "problem_config.h"
#include "solver.h"
#include "../hermes2d/include/function/exact_solution.h"



class SolverInterface
{

public:
    SolverInterface() {}
    virtual ~SolverInterface() {}
    virtual void solve(FieldInfo*, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh>) = 0;
    virtual Hermes::Hermes2D::ExactSolutionScalar<double>* getSolution() = 0;
};

QT_BEGIN_NAMESPACE
Q_DECLARE_INTERFACE(SolverInterface, "agros2d.SolverInterface/1.0")
QT_END_NAMESPACE

#endif // SOLVER_INTERFACE_H
