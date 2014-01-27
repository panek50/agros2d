#include <QTextStream>
#include <QTime>
#include <QtTest/QTest>
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
#include <complex>

#include "bem.h"
#include "mesh.h"
#include "quad_tables.h"

# define J std::complex<double>(0,1)

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{
    m_fieldInfo = fieldInfo;
    m_hermesMesh = mesh;
    m_polyOrder = 1;
}

void Bem::readMesh()
{
    Hermes::Hermes2D::Element *e;
    for_all_active_elements(e, m_hermesMesh)
    {
        QList<Node *> elementNodes;
        SceneBoundary *boundary = 0;
        for (unsigned i = 0; i < e->get_nvert(); i++)
        {
            Node node(e->vn[i]->x, e->vn[i]->y);
            // qDebug() << "Boundary:" << e->vn[i]->bnd;
            if(e->vn[i]->bnd == 1)
            {
                for(int j = 0; j < Agros2D::scene()->edges->items().count(); j++)
                {
                    boundary = Agros2D::scene()->edges->items().at(j)->marker(m_fieldInfo);
                    if (boundary && (!boundary->isNone()))
                    {
                        Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
                        double value = boundary->value(boundaryType.id()).number();
                        bool isEssential = (boundaryType.essential().count() > 0);
                        Node * firstNode,  * secondNode;
                        // qDebug() << m_hermesMesh->get_num_vertex_nodes();
                        if(atoi(m_hermesMesh->get_boundary_markers_conversion().get_user_marker(e->en[i]->marker).marker.c_str()) == j)
                        {
                            if(!mesh.m_boundaryNodes.contains(node))
                            {
                                node.globalIndex = mesh.m_boundaryNodes.count();
                                mesh.m_boundaryNodes.append(node);
                                firstNode = &mesh.m_boundaryNodes.last();
                            }
                            else
                            {
                                int index =mesh.m_boundaryNodes.indexOf(node);
                                firstNode = &mesh.m_boundaryNodes[index];
                            }

                            if(!elementNodes.contains(firstNode))
                                elementNodes.append(firstNode);

                            Node nextNode(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y);

                            if(!mesh.m_boundaryNodes.contains(nextNode))
                            {
                                nextNode.globalIndex =mesh.m_boundaryNodes.count();
                                mesh.m_boundaryNodes.append(nextNode);
                                secondNode = &mesh.m_boundaryNodes.last();
                            }
                            else
                            {
                                int index =mesh.m_boundaryNodes.indexOf(nextNode);
                                // qDebug() << "Node is: " <<mesh.m_boundaryNodes[index].toString();
                                secondNode = &mesh.m_boundaryNodes[index];
                            }
                            if(!elementNodes.contains(secondNode))
                                elementNodes.append(secondNode);

                            Segment segment(firstNode, secondNode);
                            segment.m_mesh = & mesh;
                            segment.setElement(e);
                            segment.setValue(value);
                            segment.setEssential(isEssential);
                            segment.setEdgeId(j);

                            /// ToDo: Prepare for higher order - precalculate in integration points
                            ///       Now jacobian does not depend on xi
                            segment.m_jacobian = jacobian(m_polyOrder, 0, segment);
                            segment.m_logJacobian = log(segment.m_jacobian);
                            mesh.m_segments.append(segment);
                        }
                        else
                        {
                            // qDebug() << node.toString();
                            Node * lastNode = 0;
                            if(!mesh.m_boundaryNodes.contains(node))
                            {
                                node.globalIndex =mesh.m_boundaryNodes.count();
                                mesh.m_boundaryNodes.append(node);
                                lastNode = & mesh.m_boundaryNodes.last();
                            }
                            else
                            {
                                int index = mesh.m_boundaryNodes.indexOf(node);
                                lastNode =  & mesh.m_boundaryNodes[index];
                            }

                            if((lastNode) && (!elementNodes.contains(lastNode)))
                                elementNodes.append(lastNode);
                        }
                    }
                }
            }
            else
            {
                if(!mesh.m_innerNodes.contains(node))
                {
                    mesh.m_innerNodes.append(node);
                    elementNodes.append(&mesh.m_innerNodes.last());
                }
                else
                {
                    int index = mesh.m_innerNodes.indexOf(node);
                    elementNodes.append(&mesh.m_innerNodes[index]);
                }
            }
        }
        Element element(elementNodes);
        element.setArea(e->get_area());
        SceneLabel * label =  Agros2D::scene()->labels->at(atoi(m_hermesMesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str()));
        //        element.setValue(label->marker(m_fieldInfo)->value(m_fieldInfo->materialTypeVariables().at(0).id()).number());
        element.setValue(1.);
        mesh.m_elements.append(element);
    }


    for(int i = 0; i < mesh.m_boundaryNodes.count(); i++)
    {
        mesh.m_nodes.append(& mesh.m_boundaryNodes[i]);
    }

    if(m_polyOrder == 1)
    {
        for(int i = 0; i < mesh.m_segments.count(); i++)
        {

            Segment & segment = mesh.m_segments[i];
            if (segment.isEssential())
            {
                segment.firstNode().isEssential = true;
                segment.firstNode().value = segment.value();
                segment.lastNode().isEssential = true;
                segment.lastNode().value = segment.value();
            }
            else
            {
                segment.firstNode().normalDerivation = segment.derivation();
                segment.lastNode().normalDerivation = segment.derivation();
            }

            segment.m_points.clear();
            segment.m_points.append(&segment.firstNode());
            segment.m_points.append(&segment.lastNode());
        }
    }


    if(m_polyOrder == 0)
    {
        mesh.m_nodes.clear();

        for(int i = 0; i < mesh.m_segments.count(); i++)
        {
            Segment & segment = mesh.m_segments[i];

            /// ToDo: Fix memory leak
            Node * node = new Node(segment.gravity());
            if (segment.isEssential())
            {
                node->isEssential = true;
                node->value = segment.value();
            }
            else
            {
                node->normalDerivation = segment.derivation();
            }

            node->globalIndex = i;

            mesh.m_nodes.append(node);
            segment.m_points.clear();
            segment.m_points.append(node);
        }
    }

    //        qDebug() << "Inner nodes:";
    //        foreach(Node node, mesh.m_innerNodes)
    //        {
    //            qDebug() << node.toString();
    //        }

    //        qDebug() << "---------";

    //        qDebug() << "Boundary nodes:";
    //        foreach(Node node,mesh.m_boundaryNodes)
    //        {
    //            qDebug() << node.toString();
    //            qDebug() << node.globalIndex;
    //        }

    //        qDebug() << "Mesh nodes:";
    //        for(int i = 0; i < mesh.m_nodes.count(); i++)
    //        {
    //            qDebug() << mesh.m_nodes[i]->toString();
    //            qDebug() << mesh.m_nodes[i]->globalIndex;
    //        }


    //        qDebug() << "\n Elements:";
    //        foreach(Element element, mesh.m_elements)
    //        {
    //            qDebug() << element.m_nodes.count();
    //            foreach (Node * node, element.m_nodes) {
    //                qDebug() << node->toString();
    //            }
    //            qDebug() << "-------------------------";
    //        }

    //        qDebug() << "Segments:";
    //        foreach(Segment segment, mesh.m_segments)
    //        {
    //            foreach (Node * node, segment.m_points) {
    //                qDebug() << node->toString();
    //            }
    //            qDebug() << "-------------------------";
    //        }
    mesh.m_nElement = mesh.m_elements.count();
    mesh.m_nSegment = mesh.m_segments.count();
}

QString Bem::toString()
{
    QString output = "";

    foreach(Segment component, mesh.m_segments)
    {
        output += "Edge: ";
        output += QString::number(component.id());
        output +=  "\n";
        output += "Boundary conditiom type: ";
        output += QString::number(component.isEssential());
        output +=  "\n";
        output += "Boundary conditiom value: ";
        output += QString::number(component.value());
        output +=  "\n";
        output += "Edge components: \n";
        output += "First Point:";
        output += " ";
        output += QString::number(component.firstNode().x);
        output += " ";
        output += QString::number(component.firstNode().y);
        output += "\n";
        output += "Second Point:";
        output += " ";
        output += QString::number(component.lastNode().x);
        output += " ";
        output += QString::number(component.lastNode().y);
        output += "\n";

    }
    return output;
}


Node Bem::globalCoordinates(double xi, Segment segment)
{
    Node v;
    int n = segment.geometricOrder();
    BemVector Sf(n);

    double Ni[2];
    shapeFunction(n, xi, Ni);

    QList<Node> points;
    if(n == 1);
    points.append(segment.firstNode());
    points.append(segment.lastNode());


    for(int i = 0; i <= n; i++)
    {
        v = v + points[i] * Ni[i];
    }
    return v;
}


void Bem::shapeFunction(int polyOrder, double xi, double * Ni)
{
    if (polyOrder == 0)
    {
        Ni[0] = 1;
        return;
    }

    Ni[0] = 0.5 * (1 - xi);
    Ni[1] = 0.5 * (1 + xi);

    if (polyOrder == 1)
        return;

    Ni[2] = 1.0 - xi * xi;
    Ni[0] = Ni[0] - 0.5 * Ni[2];
    Ni[1] = Ni[1] - 0.5 * Ni[2];
}

BemVector Bem::shapeFunctionDerivative(int polyOrder, double xi)
{
    BemVector Dn = BemVector(polyOrder);

    Dn(0) = -0.5;
    Dn(1) = 0.5;

    if((polyOrder == 1) || (polyOrder == 0))
        return Dn;

    Dn(2) = -2.0 * xi;
    Dn(0) = Dn(0) - 0.5 * Dn(2);
    Dn(1) = Dn(1) - 0.5 * Dn(2);
    return Dn;
}

double Bem::jacobian(int polyOrder, double xi, Segment segment)
{
    double dGamma = 0;

    if((polyOrder == 1) || (polyOrder == 0))
    {
        double x = (segment.lastNode().x - segment.firstNode().x);
        double y = (segment.lastNode().y - segment.firstNode().y);

        dGamma = 0.5 * sqrt(x*x + y*y);
    }

    return dGamma;
}

Node Bem::normalVector(double xi, Segment segment)
{
    int polyOrder = segment.geometricOrder();
    BemVector Dn = shapeFunctionDerivative(polyOrder, xi);
    Node dVec;
    Node n;

    QList<Node> points;
    if(polyOrder == 1);
    {
        points.append(segment.firstNode());
        points.append(segment.lastNode());
    }

    for(int i = 0; i < polyOrder + 1; i++)
    {
        dVec = dVec  + points[i] * Dn(i);
    }
    n.x = dVec.y;
    n.y = - dVec.x;
    n = n * 1 / segment.m_jacobian;
    return n;
}


template<typename Scalar>
BemSolution<Scalar>::BemSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{
    this->mesh = mesh;
}


template<typename Scalar>
Scalar BemSolution<Scalar>::value(double x, double y) const
{
    return m_bem->getValue(x, y);
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
    sln->setSolver(m_bem);
    return sln;
}

double Bem::kernel_length(Node refNode, Segment segment, double xi)
{
    return 1;
}


double Bem::kernel_laplace2D(Node refNode, Segment segment, double xi)
{
    double r = globalCoordinates(xi, segment).distanceOf(refNode);
    return 1.0 / (2 * M_PI) * log(1/r);
}


double Bem::kernel_laplace2D_derivation(Node refNode, Segment segment, double xi)
{
    Node r = globalCoordinates(xi, segment) - refNode;
    Node n = normalVector(xi, segment);
    double rNormSquared = globalCoordinates(xi, segment).distanceOfSquared(refNode);
    return - 1.0 / (2 * M_PI) * (n.x *  r.x  + n.y * r.y) / rNormSquared;
}

std::complex<double> Bem::kernel_helmholtz2D(Node refNode, Segment segment, double xi)
{
    double r = globalCoordinates(xi, segment).distanceOf(refNode);
    return exp(J * r) / (4 * M_PI * r);
}

std::complex<double> Bem::kernel_helmholtz2D_derivation(Node refNode, Segment segment, double xi)
{
    Node r = globalCoordinates(xi, segment) - refNode;
    Node n = normalVector(xi, segment);
    double rNorm = globalCoordinates(xi, segment).distanceOf(refNode);
    return - 1.0 / (2 * M_PI ) * (n.x *  r.x  + n.y * r.y) / (rNorm * rNorm) * (1.0 - J * rNorm) * exp(J * rNorm);
}

void Bem::solveComplex()
{

}

void Bem::solve()
{
    int order = 7;
    int n = mesh.m_segments.count();

    BemMatrix A(n, n);
    BemMatrix C(n, n);
    BemVector diagonal(n);
    BemVector rsv(n);
    BemVector bp(n);

    // Loop over all nodes
    Node node;

    for (int i = 0; i < n; i++)
    {
        node = * mesh.m_nodes[i];

        // Loop over all elements
        for (int j = 0; j < n; j++)
        {
            Segment segment = mesh.m_segments[j];
            // Loop over element nodes
            for(int k = 0; k <= m_polyOrder; k++)
            {
                int index = segment.m_points[k]->globalIndex;
                double dU = 0;
                double dT = 0;

                if(i == index)
                {
                    for(int l = 0; l <= order; l++)
                    {
                        double xi;
                        double jac = segment.m_jacobian;
                        double dxdb = 2;
                        if(m_polyOrder == 0)
                        {
                            xi = gaussLeguerreCoords[order][l];
                        }
                        else
                        {
                            if(k == 0)
                            {
                                xi = 2 * gaussLeguerreCoords[order][l] - 1.0;
                                dxdb = + 2;
                            }
                            else
                            {
                                xi = 1.0 - 2 * gaussLeguerreCoords[order][l];
                                dxdb = + 2;
                            }
                        }
                        double Ni[2];
                        double Sf[2];
                        shapeFunction(m_polyOrder, gaussCoords[order][l], Ni);
                        shapeFunction(m_polyOrder, gaussLeguerreCoords[order][l], Sf);
                        dU  +=   - 1 / ( 2 * M_PI) *  Ni[k] * jac * segment.m_logJacobian * gaussWeights[order][l] + 1 / (2 * M_PI) * dxdb  * Sf[k] * jac * gaussLeguerreWeights[order][l];
                    }

                    if(segment.m_points[k]->isEssential)
                    {
                        A(i, index) +=  -dU;
                    }
                    else
                    {
                        C(i, index) += dU;
                    }

                } else
                {
                    for(int l = 0; l <= order; l++)
                    {
                        double Ni[2];
                        shapeFunction(m_polyOrder, gaussCoords[order][l], Ni);
                        dU += Ni[k] * kernel_laplace2D(node, segment, gaussCoords[order][l]) * segment.m_jacobian * gaussWeights[order][l];
                        dT += Ni[k] * kernel_laplace2D_derivation(node, segment, gaussCoords[order][l]) * segment.m_jacobian * gaussWeights[order][l];
                    }

                    diagonal(i) += dT;
                    if(segment.m_points[k]->isEssential)
                    {
                        A(i, index) -=  dU;
                        C(i, index) -=  dT;
                    }
                    else
                    {
                        A(i, index) +=  dT;
                        C(i, index) +=  dU;
                    }
                }
            }
        }
    }

    for(int i = 0; i < n; i++)
    {
        if(mesh.m_nodes[i]->isEssential)
        {
            C(i, i) =  diagonal(i);
        }
        else
        {
            A(i,i) =  - diagonal(i);
        }
    }


    //                qDebug() << diagonal.toString();
    //                qDebug() << A.toString();
    //                qDebug() << C.toString();


    //    int m = mesh.m_elements.count();
    //    for(int i = 0; i < n; i++ )
    //    {
    //        for(int j = 0; j < m; j++)
    //        {
    //            double R = sqrt((node.x - mesh.m_elements[j].gravity().x) * (node.x - mesh.m_elements[j].gravity().x) + (node.y - mesh.m_elements[j].gravity().y) * (node.y - mesh.m_elements[j].gravity().y));
    //            bp(i) = bp(i) + 1 / (2 * M_PI) * mesh.m_elements[j].value() * log(R) * mesh.m_elements[j].araea();
    //        }
    //    }

    // qDebug() << bp.toString();

    for (int j = 0; j < n; j++)
    {
        int index = mesh.m_nodes[j]->globalIndex;
        if(mesh.m_nodes[j]->isEssential)
        {
            rsv(index) = mesh.m_nodes[j]->value;
        }

        else
        {
            rsv(index) = mesh.m_nodes[j]->normalDerivation;
        }
    }

    //    qDebug() << rsv.toString();

    /// ToDo: Poisson - right side vector


    bp.clear();
    BemVector b(n);
    b = C * rsv + bp;

    BemVector results(n);


    results = A.solve(b);


    for (int i = 0; i < n; i++)
    {
        if(mesh.m_nodes[i]->isEssential)
        {
            mesh.m_nodes[i]->normalDerivation = results(i);

        }
        else
        {
            mesh.m_nodes[i]->value = results(i);
        }
    }

    //        for (int i = 0; i < n; i++)
    //        {
    //            qDebug() << "Index:" << i;
    //            qDebug() << mesh.m_nodes[i]->toString();
    //            qDebug() << mesh.m_nodes[i]->value;
    //            qDebug() << mesh.m_nodes[i]->normalDerivation;
    //            qDebug() << mesh.m_nodes[i]->isEssential;
    //        }


    QTime myTimer;
    myTimer.start();
    domainSolution();
    int nMilliseconds = myTimer.elapsed();
    qDebug() << nMilliseconds;

}

void Bem::domainSolution()
{    
    for(int i = 0; i < mesh.m_innerNodes.count(); i++)
    {
        mesh.m_innerNodes[i].value = potentialInner(mesh.m_innerNodes.at(i).x,  mesh.m_innerNodes.at(i).y);
    }

    for(int i = 0; i < mesh.m_boundaryNodes.count(); i++)
    {
        mesh.m_boundaryNodes[i].value = potentialBoundary(mesh.m_boundaryNodes.at(i).x,  mesh.m_boundaryNodes.at(i).y);
        // qDebug() << potential(mesh.m_nodes.at(i)->x,  mesh.m_nodes.at(i)->y);
    }

    //   qDebug() << "Solution, elements:";
    //    foreach(Element element, mesh.m_elements)
    //    {
    //        for(int i = 0; i < element.m_nodes.count(); i++)
    //        {
    //            element.nodeValues[i] = element.m_nodes.at(i)->value;
    //            //           qDebug() << element.nodeValues[i];
    //        }
    //    }
}

double Bem::getValue(double x, double y)
{
    double result = 10;
    int n = mesh.m_elements.count();
    foreach(Element element, mesh.m_elements)
    {
        Node p(x,y);
        if(element.containsPoint(p))
        {
            result = element.value(p);
        }
    }

    return result;
}

double Bem::potentialInner(double x, double y)
{
    int order = 3;
    int n = mesh.m_nSegment;

    double u = 0;
    Node p(x, y);
    double dU = 0;
    double dT = 0;


    for(int i = 0; i < n; i++)
    {
        Segment segment = mesh.m_segments[i];
        double delta_u = 0;
        for(int k = 0; k <= m_polyOrder; k++)
        {
            for(int l = 0; l <= order; l++)
            {
                double  jac = segment.m_jacobian;
                double Ni[2];
                shapeFunction(m_polyOrder, gaussCoords[order][l], Ni);
                dT += Ni[k] * kernel_laplace2D_derivation(p, segment, gaussCoords[order][l]) * jac * gaussWeights[order][l];
                dU += Ni[k] * kernel_laplace2D(p, segment, gaussCoords[order][l]) * jac * gaussWeights[order][l];
            }
            delta_u =   - dT * segment.m_points[k]->value + dU * segment.m_points[k]->normalDerivation;
            u = u + delta_u;
            dT = 0;
            dU = 0;
        }
    }

    //    int m = mesh.m_elements.count();
    //    for(int i = 0; i < m; i++)
    //    {
    //        double R = sqrt((p.x - mesh.m_elements[i].gravity().x) * (p.x - mesh.m_elements[i].gravity().x) + (p.y - mesh.m_elements[i].gravity().y) * (p.y - mesh.m_elements[i].gravity().y));
    //        u = u - 1 / (2 * M_PI) * mesh.m_elements[i].value() * log(R) * mesh.m_elements[i].araea();
    //    }

    // qDebug("Time elapsed: %f ms", t.elapsed());
    return u;
}

double Bem::potentialBoundary(double x, double y)
{    
    int n = mesh.m_nSegment;
    Node p(x, y);

    for(int i = 0; i < n; i++)
    {
        Segment segment = mesh.m_segments[i];
        if(p.distanceOf(segment.firstNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return segment.m_points[0]->value;

            return segment.firstNode().value;
        }

        if(p.distanceOf(segment.lastNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return segment.m_points[0]->value;

            return segment.lastNode().value;
        }

        if(p.distanceOf(segment.gravity()) < EPS_ZERO)
        {
            return segment.m_points[0]->value;
        }
    }
}




template class BemSolution<double>;
