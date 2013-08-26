#include <QTextStream>s

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

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{
    qDebug() << "BEM solver";
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
    m_solution = new BemSolution<double>(mesh);
    addPhysics();
    assemblyMatrices();
}

BemSolution<double> * Bem::getSolution()
{
    return m_solution;
}


void Bem::addPhysics()
{


    int index = 0;
    qDebug() << (m_mesh)->get_max_element_id();
    foreach(SceneEdge* edge, Agros2D::scene()->edges->items())
    {
        SceneBoundary *boundary = edge->marker(m_fieldInfo);
       // qDebug() << boundary->name();
       // qDebug() << boundary->type();
       // qDebug() << boundary->value();
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            Hermes::Hermes2D::Element *e;

            for_all_active_elements(e, m_mesh)
            {
                for (unsigned edge = 0; edge < e->get_nvert(); edge++)
                {
                    bool boundary = false;

                    if (e->en[edge]->marker != -1)
                    {
                        if ((e->en[edge]->bnd == 1) && (!e->en[edge]->elem[1]))
                        {
                            boundary = true;
                            if(atoi(m_fieldInfo->initialMesh()->get_boundary_markers_conversion().get_user_marker(e->en[edge]->marker).marker.c_str()) == index)
                            {
                                qDebug() << e->en[edge]->id;
                                qDebug() << atoi(m_fieldInfo->initialMesh()->get_boundary_markers_conversion().get_user_marker(e->en[edge]->marker).marker.c_str());
                            }
                        }
                    }
                }
            }
        }
        index++;
    }
}

QString Bem::toString()
{
    QString meshString;
    meshString += "\nNodes: \n";
    foreach(Node node, m_nodeList)
    {
        meshString += QString::number(node.m_index);
        meshString += "  ";
        meshString += QString::number(node.m_x);
        meshString += "  ";
        meshString += QString::number(node.m_y);
        meshString += "\n";
    }

    meshString += "\nEdges \n";
    foreach(Edge edge, m_edgeList)
    {
        meshString += QString::number(edge.m_index);
        meshString += "  ";
        meshString += QString::number(edge.m_marker);
        meshString += "  ";
        meshString += QString::number(edge.m_nodes[0]);
        meshString += "  ";
        meshString += QString::number(edge.m_nodes[1]);
        meshString += "  ";
        meshString += QString::number(edge.m_boundary);
        meshString += "\n";
    }

    meshString += "\nElements \n";
    foreach(Element element, m_elementList)
    {
        meshString += QString::number(element.m_index);
        meshString += "  ";
        meshString += QString::number(element.m_marker);
        meshString += "  ";
        meshString += QString::number(element.m_nodes[0]);
        meshString += "  ";
        meshString += QString::number(element.m_nodes[1]);
        meshString += "  ";
        meshString += QString::number(element.m_nodes[2]);
        meshString += "\n";
    }
    return meshString;
}

void Bem::assemblyMatrices()
{
    int n = m_bounderyList.count();
}


template<typename Scalar>
BemSolution<Scalar>::BemSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{
    this->mesh = mesh;
}


template<typename Scalar>
Scalar BemSolution<Scalar>::value(double x, double y) const
{
    return 10;
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

template class BemSolution<double>;
