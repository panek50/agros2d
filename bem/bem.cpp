#include <QTextStream>

#include "bem.h"
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


#include <cblas.h>
#include <lapacke.h>

EdgeComponent::EdgeComponent(Node firstNode, Node secondNode)
{
    m_firstNode = firstNode;
    m_secondNode = secondNode;
    m_gravity(0) = (firstNode(0) + secondNode(0)) / 2;
    m_gravity(1) = (firstNode(1) + secondNode(1)) / 2;
    m_length = sqrt(pow(firstNode(0) - secondNode(0), 2) + pow(firstNode(1) - secondNode(1), 2));
}

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{
    qDebug() << "BEM solver";
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
    m_solution = new BemSolution<double>(mesh);
    addPhysics();
    qDebug() << toString();
}

// ToDo: probably memory leak
BemSolution<double> * Bem::getSolution()
{
    return m_solution;
}


void Bem::addPhysics()
{

    Hermes::Hermes2D::Element *e;
    int j = 0;


    foreach(SceneEdge* sceneEdge, Agros2D::scene()->edges->items())
    {
        SceneBoundary *boundary = sceneEdge->marker(m_fieldInfo);
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            double value = boundary->value(boundaryType.id()).number();
            bool isEssential = (boundaryType.essential().count() > 0);

            int nElement = 0;
            for_all_active_elements(e, m_mesh)
            {
                nElement++;
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    if(e->vn[i]->bnd == 1)
                    {
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Node firstNode, secondNode;

                            firstNode.id = e->vn[i]->id;
                            firstNode(0) = e->vn[i]->x;
                            firstNode(1) = e->vn[i]->y;

                            secondNode.id = e->vn[e->next_vert(i)]->id;
                            secondNode(0) = e->vn[e->next_vert(i)]->x;
                            secondNode(1) = e->vn[e->next_vert(i)]->y;

                            EdgeComponent component(firstNode, secondNode);

                            component.m_element = e;
                            component.m_value = value;
                            component.m_isEssential = isEssential;
                            component.m_edgeID = j;
                            m_edgeComponents.append(component);
                        }
                    }
                }
            }
            j++;
            // ToDo: improve
            m_nElement = nElement;
            qDebug() << nElement;
        }
    }
}

QString Bem::toString()
{
    QString output = "";

    foreach(EdgeComponent component, m_edgeComponents)
    {
        output += "Edge: ";
        output += QString::number(component.m_edgeID);
        output +=  "\n";
        output += "Boundary conditiom type: ";
        output += QString::number(component.m_isEssential);
        output +=  "\n";
        output += "Boundary conditiom value: ";
        output += QString::number(component.m_value);
        output +=  "\n";
        output += "Edge components: \n";
        output += "First node:";
        output += QString::number(component.firstNode().id);
        output += " ";
        output += QString::number(component.firstNode()(0));
        output += " ";
        output += QString::number(component.firstNode()(1));
        output += "\n";
        output += "Second node:";
        output += QString::number(component.secondNode().id);
        output += " ";
        output += QString::number(component.secondNode()(0));
        output += " ";
        output += QString::number(component.secondNode()(1));
        output += "\n";

    }
    qDebug() << output;
    return output;
}


template<typename Scalar>
BemSolution<Scalar>::BemSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{
    this->mesh = mesh;
}


template<typename Scalar>
Scalar BemSolution<Scalar>::value(double x, double y) const
{
    return x*y;
}

template<typename Scalar>
void BemSolution<Scalar>::derivatives(double x, double y, Scalar &dx, Scalar &dy) const
{
    dx = 0;
    dy = 0;
}

template<typename Scalar>
Hermes::Hermes2D::MeshFunction<Scalar>* BemSolution<Scalar>::clone() const
{
    if(this->sln_type == Hermes::Hermes2D::HERMES_SLN)
        return Hermes::Hermes2D::Solution<Scalar>::clone();
    BemSolution<Scalar>* sln = new BemSolution<Scalar>(this->mesh);
    return sln;
}

void Bem::solve()
{
    int n = m_edgeComponents.count();
    double * matrix_H = new double[n * n];
    double * matrix_G = new double[n * n];


    delete[] matrix_H;
    delete[] matrix_G;
}

template class BemSolution<double>;
