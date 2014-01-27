#include <QTextStream>


#include "util.h"
#include "util/global.h"

#include "../agros2d-library/scene.h"
#include "../agros2d-library/scenemarker.h"
#include "../agros2d-library/scenebasic.h"
#include "../agros2d-library/scenenode.h"
#include "../agros2d-library/sceneedge.h"
#include "../agros2d-library/scenelabel.h"
#include "../agros2d-library/hermes2d/module.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"
#include "../hermes2d/include/function/exact_solution.h"

#include "bem_interface.h"


Hermes::Hermes2D::ExactSolutionScalar<double>* BemInterface::getSolution()
{            
    return  m_solution;
}

void BemInterface::solve(FieldInfo* field, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh> mesh)
{            
    if(m_bem.isNull())
    {
        m_bem = QSharedPointer<Bem>(new Bem(field, mesh));
    }
    else
    {
        m_bem.clear();
        m_bem = QSharedPointer<Bem>(new Bem(field, mesh));
    }
    m_bem->readMesh();
    m_bem->solve();
    // m_bem->solveComplex();

    m_solution =  new BemSolution<double>(mesh);
    m_solution->setSolver(m_bem);
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(agros2d_bem_interface, BemInterface)
#endif
