#ifndef BEM_INTERFACE_H
#define BEM_INTERFACE_H

#include <QtPlugin>

#include "../util/util.h"
#include "../agros2d-library/hermes2d/solver_interface.h"
#include "bem.h"

class BemInterface: public QObject, SolverInterface
{
    Q_OBJECT
    Q_INTERFACES(SolverInterface)
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
    Q_PLUGIN_METADATA(IID "org.hpfem.agros2d.BemInterface" FILE "")
#endif
public:
    BemInterface() { qDebug() << "BemInterface"; }
    virtual void solve(FieldInfo*, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh>);
    virtual Hermes::Hermes2D::ExactSolutionScalar<double>* getSolution();
private:
    std::tr1::shared_ptr<Bem> m_bem;
};

#endif // BEM_INTERFACE_H
