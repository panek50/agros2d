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


EdgeComponent::EdgeComponent(Node firstNode, Node secondNode)
{
    m_firstNode = firstNode;
    m_secondNode = secondNode;
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
        Edge edge = Edge(j);
        SceneBoundary *boundary = sceneEdge->marker(m_fieldInfo);
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            edge.m_value = boundary->value(boundaryType.id()).number();
            edge.m_isEssential = (boundaryType.essential().count() > 0);

            for_all_active_elements(e, m_mesh)
            {
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    if(e->vn[i]->bnd == 1)
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Node firstNode, secondNode;
                            firstNode.id = e->vn[i]->id;
                            firstNode.x = e->vn[i]->x;
                            firstNode.y = e->vn[i]->y;

                            secondNode.id = e->vn[e->next_vert(i)]->id;
                            secondNode.x = e->vn[e->next_vert(i)]->x;
                            secondNode.y = e->vn[e->next_vert(i)]->y;

                            EdgeComponent component(firstNode, secondNode);
                            edge.m_components.append(component);
                        }
                }
            }
            m_geometry.m_edges.append(edge);
            j++;
        }
    }
}


QString Bem::toString()
{
    QString output = "";

    foreach (Edge edge, m_geometry.m_edges)
    {
        output += "Edge: ";
        output += QString::number(edge.id());
        output +=  "\n";
        output += "Boundary conditiom type: ";
        output += QString::number(edge.m_isEssential);
        output +=  "\n";
        output += "Boundary conditiom value: ";
        output += QString::number(edge.m_value);
        output +=  "\n";
        output += "Edge components: \n";
        foreach(EdgeComponent component, edge.m_components)
        {
            output += "First node:";
            output += QString::number(component.firstNode().id);
            output += " ";
            output += QString::number(component.firstNode().x);
            output += " ";
            output += QString::number(component.firstNode().y);
            output += "\n";
            output += "Second node:";
            output += QString::number(component.secondNode().id);
            output += " ";
            output += QString::number(component.secondNode().x);
            output += " ";
            output += QString::number(component.secondNode().y);
            output += "\n";
        }
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

Bem::solve()
{

}

template class BemSolution<double>;
