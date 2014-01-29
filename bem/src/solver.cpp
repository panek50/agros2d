// Bem headers
#include "solver.h"
#include "mesh.h"
#include "quad_tables.h"

// Agros2D headers - common stuff
#include "util.h"
#include "util/global.h"

// Agros2D-library
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

// Hermes libraries
#include "../hermes2d/include/function/exact_solution.h"

// 3d party mathematical libraries
#include <cblas.h>
#include <lapacke.h>

// QT headers
#include <QTextStream>
#include <QTime>
#include <QtTest/QTest>

// Std headers
#include <complex>

// complex unit definition
const std::complex<double> J(0,1);

template <typename Type>
Solver<Type>::Solver(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{
    m_fieldInfo = fieldInfo;
    m_hermesMesh = mesh;
    m_polyOrder = 1;
}

template <typename Type>
void Solver<Type>::readMesh()
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

                        Node *firstNode;
                        Node *secondNode;

                        if(atoi(m_hermesMesh->get_boundary_markers_conversion().get_user_marker(e->en[i]->marker).marker.c_str()) == j)
                        {
                            if(!mesh.m_points.contains(node)) {
                                node.globalIndex = mesh.m_boundaryNodes.count();
                                mesh.m_points.append(node);
                                firstNode =  & mesh.m_points.last();
                                mesh.m_boundaryNodes.append(firstNode);
                            }
                            else
                            {
                                int index = mesh.m_points.indexOf(node);
                                firstNode = & mesh.m_points[index];
                            }

                            if(!elementNodes.contains(firstNode))
                                elementNodes.append(firstNode);

                            Node nextNode(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y);

                            if(!mesh.m_points.contains(nextNode))
                            {
                                nextNode.globalIndex = mesh.m_boundaryNodes.count();
                                mesh.m_points.append(nextNode);
                                mesh.m_boundaryNodes.append(& mesh.m_points.last());
                                secondNode = & mesh.m_points.last();
                            }
                            else
                            {
                                int index =mesh.m_points.indexOf(nextNode);
                                // qDebug() << "Node is: " << mesh.m_boundaryNodes[index].toString();
                                secondNode = & mesh.m_points[index];
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
                            if(!mesh.m_points.contains(node))
                            {
                                node.globalIndex = mesh.m_boundaryNodes.count();
                                mesh.m_points.append(node);
                                mesh.m_boundaryNodes.append(& mesh.m_points.last());
                                lastNode = & mesh.m_points.last();
                            }
                            else
                            {
                                int index = mesh.m_points.indexOf(node);
                                lastNode = & mesh.m_points[index];
                            }

                            if((lastNode) && (!elementNodes.contains(lastNode)))
                                elementNodes.append(lastNode);
                        }
                    }
                }
            }
            else
            {
                if(!mesh.m_points.contains(node))
                {
                    mesh.m_points.append(node);
                    mesh.m_innerNodes.append(&mesh.m_points.last());
                    elementNodes.append(&mesh.m_points.last());
                }
                else
                {
                    int index = mesh.m_points.indexOf(node);
                    elementNodes.append(& mesh.m_points[index]);
                }
            }
        }
        Element element(elementNodes);
        element.setArea(e->get_area());
        // SceneLabel * label =  Agros2D::scene()->labels->at(atoi(m_hermesMesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str()));
        //        element.setValue(label->marker(m_fieldInfo)->value(m_fieldInfo->materialTypeVariables().at(0).id()).number());
        element.setValue(1.);
        mesh.m_elements.append(element);
    }


    if(m_polyOrder == 1)
    {
        for(int i = 0; i < mesh.m_boundaryNodes.count(); i++)
        {
            mesh.m_integrationPoints.append(mesh.m_boundaryNodes[i]);
        }

        for(int i = 0; i < mesh.m_segments.count(); i++)
        {

            Segment & segment = mesh.m_segments[i];
            if (segment.isEssential())
            {
                segment.firstNode().isEssential = true;
                segment.firstNode().real = segment.value();
                segment.lastNode().isEssential = true;
                segment.lastNode().real = segment.value();
            }
            else
            {
                segment.firstNode().normalDerivationReal = segment.derivation();
                segment.lastNode().normalDerivationReal = segment.derivation();
            }

            segment.m_points.clear();
            segment.m_points.append(&segment.firstNode());
            segment.m_points.append(&segment.lastNode());
        }
    }


    if(m_polyOrder == 0)
    {
        mesh.m_integrationPoints.clear();
        for(int i = 0; i < mesh.m_segments.count(); i++)
        {
            Segment & segment = mesh.m_segments[i];

            /// ToDo: Fix memory leak
            Node  node(segment.gravity());
            if (segment.isEssential())
            {
                node.isEssential = true;
                node.real = segment.value();
            }
            else
            {
                node.normalDerivationReal = segment.derivation();
            }

            node.globalIndex = i;
            mesh.m_points.append(node);
            mesh.m_integrationPoints.append(& mesh.m_points.last());

            segment.m_points.clear();
            segment.m_points.append(& mesh.m_points.last());
        }
    }

    //    qDebug() << "Inner nodes:";
    //    foreach(Node * node, mesh.m_innerNodes)
    //    {
    //        qDebug() << node->toString();
    //    }

    //    qDebug() << "---------";

    //    qDebug() << "Boundary nodes:";
    //    foreach(Node * node,mesh.m_boundaryNodes)
    //    {
    //        qDebug() << node->toString();
    //        qDebug() << node->globalIndex;
    //    }

    //    qDebug() << "Mesh nodes:";
    //    for(int i = 0; i < mesh.m_integrationPoints.count(); i++)
    //    {
    //        qDebug() << mesh.m_integrationPoints[i]->toString();
    //        qDebug() << mesh.m_integrationPoints[i]->globalIndex;
    //    }


    //    qDebug() << "\n Elements:";
    //    foreach(Element element, mesh.m_elements)
    //    {
    //        qDebug() << element.m_nodes.count();
    //        foreach (Node * node, element.m_nodes) {
    //            qDebug() << node->toString();
    //        }
    //        qDebug() << "-------------------------";
    //    }

    //    qDebug() << "Segments:";
    //    foreach(Segment segment, mesh.m_segments)
    //    {
    //        foreach (Node * node, segment.m_points) {
    //            qDebug() << node->toString();
    //        }
    //        qDebug() << "-------------------------";
    //    }
    mesh.m_nElement = mesh.m_elements.count();
    mesh.m_nSegment = mesh.m_segments.count();
}

template <typename Type>
QString Solver<Type>::toString()
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

template <typename Type>
Node Solver<Type>::globalCoordinates(double xi, Segment segment)
{
    Node v;
    int n = segment.geometricOrder();

    double Ni[2];
    shapeFunction(n, xi, Ni);

    // Slow - rewrite
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

template <typename Type>
void Solver<Type>::shapeFunction(int polyOrder, double xi, double * Ni)
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


template <typename Type>
BemVector<double> Solver<Type>::shapeFunctionDerivative(int polyOrder, double xi)
{
    BemVector<double> Dn = BemVector<double>(polyOrder);

    Dn(0) = -0.5;
    Dn(1) = 0.5;

    if((polyOrder == 1) || (polyOrder == 0))
        return Dn;

    Dn(2) = -2.0 * xi;
    Dn(0) = Dn(0) - 0.5 * Dn(2);
    Dn(1) = Dn(1) - 0.5 * Dn(2);
    return Dn;
}

template <typename Type>
double Solver<Type>::jacobian(int polyOrder, double xi, Segment segment)
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

template <typename Type>
Node Solver<Type>::normalVector(double xi, Segment segment)
{
    int polyOrder = segment.geometricOrder();
    BemVector<double> Dn = shapeFunctionDerivative(polyOrder, xi);
    Node dVec;
    Node n;

    /// Todo: Rewrite to slow
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

template <typename Scalar>
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

template <typename Type>
double Solver<Type>::kernel_length(Node refNode, Segment segment, double xi)
{
    return 1;
}


template <typename Type>
double Solver<Type>::kernel_laplace2D(Node refNode, Segment segment, double xi)
{
    double r = globalCoordinates(xi, segment).distanceOf(refNode);
    return 1.0 / (2 * M_PI) * log(1/r);
}


template <typename Type>
double Solver<Type>::kernel_laplace2D_derivation(Node refNode, Segment segment, double xi)
{
    Node r = globalCoordinates(xi, segment) - refNode;
    Node n = normalVector(xi, segment);
    double rNormSquared = globalCoordinates(xi, segment).distanceOfSquared(refNode);
    return - 1.0 / (2 * M_PI) * (n.x *  r.x  + n.y * r.y) / rNormSquared;
}

template <typename Type>
std::complex<double> Solver<Type>::kernel_helmholtz2D(Node refNode, Segment segment, double xi)
{
    double r = globalCoordinates(xi, segment).distanceOf(refNode);
    return exp(J * r) / (4 * M_PI * r);
}

template <typename Type>
std::complex<double> Solver<Type>::kernel_helmholtz2D_derivation(Node refNode, Segment segment, double xi)
{
    Node r = globalCoordinates(xi, segment) - refNode;
    Node n = normalVector(xi, segment);
    double rNorm = globalCoordinates(xi, segment).distanceOf(refNode);
    return - 1.0 / (2 * M_PI ) * (n.x *  r.x  + n.y * r.y) / (rNorm * rNorm) * (1.0 - J * rNorm) * exp(J * rNorm);
}


template <>
void Solver<double>::domainSolution()
{
    for(int i = 0; i < mesh.m_innerNodes.count(); i++)
    {
        mesh.m_innerNodes[i]->real = solutionInner(mesh.m_innerNodes.at(i)->x,  mesh.m_innerNodes.at(i)->y);
    }

    for(int i = 0; i < mesh.m_boundaryNodes.count(); i++)
    {
        mesh.m_boundaryNodes[i]->real = solutionBoundary(mesh.m_boundaryNodes.at(i)->x,  mesh.m_boundaryNodes.at(i)->y);
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

template <>
void Solver<std::complex<double> >::domainSolution()
{
    for(int i = 0; i < mesh.m_innerNodes.count(); i++)
    {
        mesh.m_innerNodes[i]->real = solutionInner(mesh.m_innerNodes.at(i)->x,  mesh.m_innerNodes.at(i)->y).real();
        mesh.m_innerNodes[i]->imag = solutionInner(mesh.m_innerNodes.at(i)->x,  mesh.m_innerNodes.at(i)->y).imag();
    }

    for(int i = 0; i < mesh.m_boundaryNodes.count(); i++)
    {
        mesh.m_boundaryNodes[i]->real = solutionBoundary(mesh.m_boundaryNodes.at(i)->x,  mesh.m_boundaryNodes.at(i)->y).imag();
        mesh.m_boundaryNodes[i]->imag = solutionBoundary(mesh.m_boundaryNodes.at(i)->x,  mesh.m_boundaryNodes.at(i)->y).imag();
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

template <typename Type>
void Solver<Type>::solve()
{
    int order = 7;
    int n = mesh.m_nSegment;

    BemMatrix<Type> A(n, n);
    BemMatrix<Type> C(n, n);
    BemVector<Type> diagonal(n);
    BemVector<Type> rsv(n);
    BemVector<Type> bp(n);


    // Loop over all nodes
    for (int i = 0; i < n; i++)
    {
        Node node = *mesh.m_integrationPoints[i];

        // Loop over all segments
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

    //    qDebug() << diagonal.toString();
    //    qDebug() << A.toString();
    //    qDebug() << C.toString();


    for(int i = 0; i < n; i++)
    {
        if(mesh.m_integrationPoints[i]->isEssential)
        {
            C(i, i) =  diagonal(i);
        }
        else
        {
            A(i,i) =  - diagonal(i);
        }
    }


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
        int index = mesh.m_integrationPoints[j]->globalIndex;
        if(mesh.m_integrationPoints[j]->isEssential)
        {
            rsv(index) = mesh.m_integrationPoints[j]->real;
        }

        else
        {
            rsv(index) = mesh.m_integrationPoints[j]->normalDerivationReal;
        }
    }

    //    qDebug() << rsv.toString();

    /// ToDo: Poisson - right side vector


    bp.clear();
    BemVector<Type> b(n);
    b = C * rsv + bp;
    BemVector<Type> results(n);
    results = A.solve(b);




    //        for (int i = 0; i < n; i++)
    //        {
    //            qDebug() << "Index:" << i;
    //            qDebug() << mesh.m_nodes[i]->toString();
    //            qDebug() << mesh.m_nodes[i]->value;
    //            qDebug() << mesh.m_nodes[i]->normalDerivation;
    //            qDebug() << mesh.m_nodes[i]->isEssential;
    //        }

    fillResults(results);
    QTime myTimer;
    myTimer.start();
    domainSolution();
    int nMilliseconds = myTimer.elapsed();
    qDebug() << nMilliseconds;

}


template <>
void Solver<double>::fillResults(BemVector<double> & results) const
{
    for(int i = 0; i < mesh.m_nSegment; i++)
    {
        if(mesh.m_integrationPoints[i]->isEssential)
        {
            mesh.m_integrationPoints[i]->normalDerivationReal = results(i);

        }
        else
        {
            mesh.m_integrationPoints[i]->real = results(i);
        }
    }
}

template <>
void Solver<std::complex<double> >::fillResults(BemVector<std::complex<double> > & results) const
{
    for(int i = 0; i < mesh.m_nSegment; i++)
    {
        if(mesh.m_integrationPoints[i]->isEssential)
        {
            mesh.m_integrationPoints[i]->normalDerivationReal = results(i).real();
            mesh.m_integrationPoints[i]->normalDerivationImag = results(i).imag();

        }
        else
        {
            mesh.m_integrationPoints[i]->real = results(i).real();
            mesh.m_integrationPoints[i]->imag = results(i).imag();
        }
    }
}

template <>
double Solver<double>::getValue(double x, double y)
{
    double result;
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

template <>
double Solver<std::complex<double> >::getValue(double x, double y)
{
    double result;
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

template <typename Type>
Type Solver<Type>::solutionInner(double x, double y)
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
            delta_u =   - dT * segment.m_points[k]->real + dU * segment.m_points[k]->normalDerivationReal;
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

template <typename Type>
Type Solver<Type>::solutionBoundary(double x, double y)
{    
    int n = mesh.m_nSegment;
    Node p(x, y);

    for(int i = 0; i < n; i++)
    {
        Segment segment = mesh.m_segments[i];
        if(p.distanceOf(segment.firstNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return segment.m_points[0]->real;

            return segment.firstNode().real;
        }

        if(p.distanceOf(segment.lastNode()) < EPS_ZERO)
        {
            if(m_polyOrder == 0)
                return segment.m_points[0]->real;

            return segment.lastNode().real;
        }

        if(p.distanceOf(segment.gravity()) < EPS_ZERO)
        {
            return segment.m_points[0]->real;
        }
    }
}




template class BemSolution<double>;
template class Solver<double>;
template class Solver<std::complex<double> >;
