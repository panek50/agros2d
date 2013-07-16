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
#include "algebra.h"

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{
    qDebug() << "BEM solver";
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
    readMesh();
    addPhysics();
    assemblyMatrices();
}

void Bem::readMesh()
{
    for (int i = 0; i<Agros2D::scene()->edges->length(); i++)
    {
    }
    //    QTextStream qout(stdout);
    //    QString meshFileName = cacheProblemDir() + "/initial.mesh";
    //    QFile * meshFile = new QFile(meshFileName);
    //    meshFile->open(QIODevice::ReadOnly | QIODevice::Text);
    //    QDomDocument document;
    //    document.setContent(meshFile);

    //    QDomNodeList vertices = document.elementsByTagName("v");
    //    for(int i = 0; i < vertices.count(); i++)
    //    {
    //        Node node;
    //        node.m_index =  vertices.at(i).attributes().namedItem("i").nodeValue().toInt();
    //        node.m_x = vertices.at(i).attributes().namedItem("x").nodeValue().toDouble();
    //        node.m_y = vertices.at(i).attributes().namedItem("y").nodeValue().toDouble();
    //        m_nodeList.append(node);
    //    }

    //    QDomNodeList elements = document.elementsByTagName("domain:t");
    //    for(int i = 0; i < elements.count(); i++)
    //    {
    //        Element element;
    //        element.m_index = elements.at(i).attributes().namedItem("i").nodeValue().toInt();
    //        element.m_marker = elements.at(i).attributes().namedItem("m").nodeValue().toInt();
    //        element.m_nodes[0] = elements.at(i).attributes().namedItem("v1").nodeValue().toDouble();
    //        element.m_nodes[1] = elements.at(i).attributes().namedItem("v2").nodeValue().toDouble();
    //        element.m_nodes[2] = elements.at(i).attributes().namedItem("v3").nodeValue().toDouble();
    //        m_elementList.append(element);
    //    }

    //    QDomNodeList edges = document.elementsByTagName("ed");
    //    for(int i = 0; i < edges.count(); i++)
    //    {
    //        Edge edge;
    //        edge.m_index = edges.at(i).attributes().namedItem("i").nodeValue().toInt();
    //        edge.m_marker = edges.at(i).attributes().namedItem("m").nodeValue().toInt();
    //        edge.m_nodes[0] =  edges.at(i).attributes().namedItem("v1").nodeValue().toDouble();
    //        edge.m_nodes[1] =  edges.at(i).attributes().namedItem("v2").nodeValue().toDouble();
    //        edge.m_boundary = false;
    //        edge.m_boundary_type = BoundaryConditionType_Vector;
    //        edge.m_boundary_value = 0;
    //        m_edgeList.append(edge);
    //    }

    //    QDomNode boundaryEdges = document.elementsByTagName("boundary_edges").at(0);

    //    int n = boundaryEdges.childNodes().count();
    //    for(int i = 0; i < n; i++)
    //    {
    //        int index = boundaryEdges.childNodes().at(i).childNodes().at(0).nodeValue().toInt();
    //        for(int j = 0; j <  m_edgeList.count(); j++)
    //        {
    //            if (m_edgeList.at(j).m_index == index)
    //            {
    //                m_edgeList[j].m_boundary = true;
    //                m_bounderyList.append(m_edgeList[j]);
    //            }
    //        }
    //    }
}

void Bem::addPhysics()
{
    int index = 0;
    foreach(SceneEdge* edge, Agros2D::scene()->edges->items())
    {
        SceneBoundary *boundary = edge->marker(m_fieldInfo);
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            for(int i = 0; i < m_edgeList.count(); i++)
            {
                if (m_edgeList[i].m_marker == index)
                {
                    m_edgeList[i].m_boundary_value = boundary->values().values()[0].number();
                }
            }

            foreach (FormInfo form, boundaryType.essential())
            {
                for(int i = 0; i < m_edgeList.count(); i++)
                {
                    if (m_edgeList[i].m_marker == index)
                    {
                        m_edgeList[i].m_boundary_type = BoundaryConditionType_Essential;
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
    Algebra::Matrix H(n, n);
    Algebra::Matrix G(n, n);
    for(int i = 0; i < n; i++)
    {
        H.setValue(i, i, 0);
        qDebug() <<  H(i,i);
    }
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

template HERMES_API class BemSolution<double>;
