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

#include "bem.h"
#include "mesh.h"

Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{    
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
}


void Bem::addPhysics()
{
    QList<Segment> edgeComponents;
    Hermes::Hermes2D::Element *e;
    int j = 0;
    for(int j = 0; j < Agros2D::scene()->edges->items().count(); j++)
    {        
        SceneBoundary *boundary = Agros2D::scene()->edges->items().at(j)->marker(m_fieldInfo);
        if (boundary && (!boundary->isNone()))
        {
            Module::BoundaryType boundaryType = m_fieldInfo->boundaryType(boundary->type());
            double value= boundary->value(boundaryType.id()).number();
            bool isEssential = (boundaryType.essential().count() > 0);

            int nElement = 0;
            for_all_active_elements(e, m_mesh)
            {
                nElement++;
                QList<Node *> points;
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    Node * point = new Node(e->vn[i]->x, e->vn[i]->y);
                    points.append(point);

                    if(e->vn[i]->bnd == 1)
                    {
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Node * firstPoint = new Node(e->vn[i]->x, e->vn[i]->y);
                            Node * secondPoint = new Node(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y);

                            int index = mesh.m_nodes.indexOf(firstPoint);
                            if(index == -1)
                            {                                
                                firstPoint->id = mesh.m_nodes.count();
                                mesh.m_nodes.append(firstPoint);
                            }
                            else
                            {                                                             
                                delete firstPoint;
                                firstPoint = mesh.m_nodes[index];
                            }
                            index = mesh.m_nodes.indexOf(secondPoint);
                            if(index == -1)
                            {                                
                                secondPoint->id = mesh.m_nodes.count();
                                mesh.m_nodes.append(secondPoint);
                            }
                            else
                            {                                
                                delete secondPoint;
                                secondPoint = mesh.m_nodes[index];
                            }

                            Segment segment(firstPoint, secondPoint);
                            segment.m_mesh = & mesh;
                            segment.setElement(e);
                            segment.setValue(value);
                            segment.setEssential(isEssential);
                            segment.setEdgeId(j);
                            edgeComponents.append(segment);
                        }
                    }
                }
                Element element(points);
                element.setArea(e->get_area());

                mesh.m_elements.append(element);
                // ToDo: read material properties

            }            
            mesh.m_nElement = nElement;
        }
    }


    mesh.m_segments.append(edgeComponents.at(0));
    for(int j = 0; j < edgeComponents.count(); j++)
    {
        for(int i = 1; i < edgeComponents.count(); i++)
        {
            if(edgeComponents[i].firstNode() == mesh.m_segments[j].secondNode())
            {
                mesh.m_segments.append(edgeComponents.at(i));
                break;
            }
        }
    }    

    foreach(Segment segment, mesh.m_segments)
    {

        if (segment.isEssential())
        {
            segment.firstNode().isEssential = true;
            segment.firstNode().value = segment.value();
            segment.secondNode().isEssential = true;
            segment.secondNode().value = segment.value();
        }
        else
        {
            segment.firstNode().normalDerivation = segment.derivation();
            segment.secondNode().normalDerivation = segment.derivation();
        }
    }

    foreach(Node * node, mesh.m_nodes)
    {
        qDebug() << node->x << node->y;
        qDebug() << node->isEssential;
        qDebug() << node->value;
        qDebug() << node->normalDerivation;
    }
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
        output += QString::number(component.secondNode().x);
        output += " ";
        output += QString::number(component.secondNode().y);
        output += "\n";

    }
    return output;
}


Node Bem::globalCoordinates(double xi, Segment segment)
{
    Node v;
    int n = segment.geometricOrder();
    BemVector Sf(n);


    Sf = shapeFunction(n, xi);

    QList<Node> points;
    if(n == 1);
    points.append(segment.firstNode());
    points.append(segment.secondNode());


    for(int i = 0; i <= n; i++)
    {
        v = v + points[i] * Sf(i);
    }
    return v;
}


BemVector Bem::shapeFunction(int polyOrder, double xi)
{     
    BemVector Ni = BemVector(polyOrder + 1);
    if (polyOrder == 0)
    {
        Ni(0) = 1;
        return Ni;
    }

    Ni(0) = 0.5 * (1 - xi);
    Ni(1) = 0.5 * (1 + xi);

    if (polyOrder == 1)
        return Ni;

    Ni(2) = 1.0 - xi * xi;
    Ni(0) = Ni(0) - 0.5 * Ni(2);
    Ni(2) = Ni(1) - 0.5 * Ni(2);
    return Ni;
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
        double x = (segment.secondNode().x - segment.firstNode().x);
        double y = (segment.secondNode().y - segment.firstNode().y);

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
        points.append(segment.secondNode());
    }

    for(int i = 0; i < polyOrder + 1; i++)
    {
        dVec = dVec  + points[i] * Dn(i);
    }
    n.x = dVec.y;
    n.y = - dVec.x;
    n = n * 1 / jacobian(polyOrder, xi, segment);
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
    double rNorm = globalCoordinates(xi, segment).distanceOf(refNode);
    return - 1.0 / (2 * M_PI) * (n.x *  r.x  + n.y * r.y) / (rNorm * rNorm);
}

void Bem::solve()
{        
    // Gauss coordinates and weights
    QList<QList<double> > weights;
    QList<QList<double> > coords;


    QList<double> wi;
    QList<double> coord;

    wi.append(2.0);
    coord.append(0);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.577350269);
    coord.append(-0.577350269);
    wi.append(1.0);
    wi.append(1.0);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.774596669);
    coord.append(0);
    coord.append(-0.774596669);
    wi.append(0.5555555555);
    wi.append(0.8888888888);
    wi.append(0.5555555555);
    weights.append(wi);
    coords.append(coord);


    wi.clear();
    coord.clear();
    coord.append(0.861136311);
    coord.append(0.339981043);
    coord.append(-0.339981043);
    coord.append(-0.861136311);
    wi.append(0.347854845);
    wi.append(0.652145154);
    wi.append(0.652145154);
    wi.append(0.347854845);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.906179845);
    coord.append(0.538469310);
    coord.append(0);
    coord.append(-0.538469310);
    coord.append(-0.906179845);

    wi.append(0.236926885);
    wi.append(0.478628670);
    wi.append(0.568888888);
    wi.append(0.478628670);
    wi.append(0.236926885);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.932469514);
    coord.append(0.661209386);
    coord.append(0.238619186);
    coord.append(-0.238619186);
    coord.append(-0.661209386);
    coord.append(-0.932469514);

    wi.append(0.171324492);
    wi.append(0.360761573);
    wi.append(0.467913934);
    wi.append(0.467913934);
    wi.append(0.360761573);
    wi.append(0.171324492);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.949107912);
    coord.append(0.741531185);
    coord.append(0.405845151);
    coord.append(0);
    coord.append(-0.405845151);
    coord.append(-0.741531185);
    coord.append(-0.949107912);
    wi.append(0.129484966);
    wi.append(0.279705391);
    wi.append(0.381830050);
    wi.append(0.417959183);
    wi.append(0.381830050);
    wi.append(0.279705391);
    wi.append(0.129484966);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.960289856);
    coord.append(0.796666477);
    coord.append(0.525532409);
    coord.append(0.183434642);
    coord.append(-0.183434642);
    coord.append(-0.525532409);
    coord.append(-0.796666477);
    coord.append(-0.960289856);
    wi.append(0.101228536);
    wi.append(0.222381034);
    wi.append(0.313706645);
    wi.append(0.362683783);
    wi.append(0.362683783);
    wi.append(0.313706645);
    wi.append(0.222381034);
    wi.append(0.101228536);
    weights.append(wi);
    coords.append(coord);

    QList<QList<double> > glWeights;
    QList<QList<double> > glCoords;


    QList<double> glWi;
    QList<double> glCoord;
    QList<double> rads;

    glWi.append(1);
    glCoord.append(0.5);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.112008806);
    glCoord.append(0.602276908);

    glWi.append(0.718539319);
    glWi.append(0.281460680);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.063890793);
    glCoord.append(0.368997063);
    glCoord.append(0.766880303);

    glWi.append(0.513404552);
    glWi.append(0.391980041);
    glWi.append(0.0946154065);
    weights.append(glWi);
    coords.append(glCoord);


    glWi.clear();
    glCoord.clear();
    glCoord.append(0.0414484801);
    glCoord.append(0.245274914);
    glCoord.append(0.556165453);
    glCoord.append(0.848982394);

    glWi.append(0.383464068);
    glWi.append(0.386875317);
    glWi.append(0.190435126);
    glWi.append(0.0392254871);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.0291344721);
    glCoord.append(0.173977213);
    glCoord.append(0.411702520);
    glCoord.append(0.677314174);
    glCoord.append(0.894771361);

    glWi.append(0.297893471);
    glWi.append(0.349776226);
    glWi.append(0.234488290);
    glWi.append(0.0989304595);
    glWi.append(0.0189115521);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.021634005);
    glCoord.append(0.129583391);
    glCoord.append(0.314020449);
    glCoord.append(0.538657217);
    glCoord.append(0.756915337);
    glCoord.append(0.922668851);

    glWi.append(0.238763662);
    glWi.append(0.308286573);
    glWi.append(0.245317426);
    glWi.append(0.142008756);
    glWi.append(0.055454622);
    glWi.append(0.0101689586);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.016719355);
    glCoord.append(0.100185677);
    glCoord.append(0.246294246);
    glCoord.append(0.433463493);
    glCoord.append(0.632350988);
    glCoord.append(0.811118626);
    glCoord.append(0.940848166);
    glWi.append(0.196169389);
    glWi.append(0.270302644);
    glWi.append(0.239681873);
    glWi.append(0.165775774);
    glWi.append(0.088943227);
    glWi.append(0.0331943043);
    glWi.append(0.0059327870);
    weights.append(glWi);
    coords.append(glCoord);

    glWi.clear();
    glCoord.clear();
    glCoord.append(0.0133202441);
    glCoord.append(0.0797504290);
    glCoord.append(0.197871029);
    glCoord.append(0.354153994);
    glCoord.append(0.529458575);
    glCoord.append(0.701814529);
    glCoord.append(0.849379320);
    glCoord.append(0.953326450);
    glWi.append(0.164416604);
    glWi.append(0.237525610);
    glWi.append(0.226841984);
    glWi.append(0.175754079);
    glWi.append(0.112924030);
    glWi.append(0.0578722107);
    glWi.append(0.0209790737);
    glWi.append(0.00368640710);
    glWeights.append(glWi);
    glCoords.append(glCoord);


    int polyOrder = 1;
    int order = 5;
    int n = mesh.m_segments.count();

    // double (Bem::*kernel)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D;
    // double (Bem::*kernel_derivation)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D_derivation;

    BemMatrix dU(n, n);
    BemMatrix dT(n, n);
    dU.clear();
    dT.clear();
    // Loop over all nodes
    Node node;
    for (int i = 0; i < n; i++)
    {
        if(polyOrder == 0)
            node = mesh.m_segments[i].gravity();
        if(polyOrder == 1)
            node = mesh.m_segments[i].firstNode();

        // Loop over all elements
        for (int j = 0; j < n; j++)
        {
            Segment segment = mesh.m_segments[j];
            // Loop over element nodes
            for(int k = 0; k <= polyOrder; k++)
            {
                {
                    // Analytic integration
                    for(int l = 0; l <= order; l++)
                    {
                        double  jac =jacobian(polyOrder, coords[order][l],segment);
                        if(i == (j + k))
                        {
                            dU(i, j + k) +=  1 / (2 * M_PI) * (- 2 * log(segment.length()) * shapeFunction(polyOrder, coords[order][l])(k) *  jac * weights[order][l] +
                                                               2 * shapeFunction(polyOrder, coords[order][l])(k) *  jac * weights[order][l]);
                            dT(i, j + k) = 0.5;
                        }
                        else
                        {
                            dT(i, j + k) += shapeFunction(polyOrder, coords[order][l])(k) * kernel_laplace2D_derivation(node, segment, coords[order][l]) * jac * weights[order][l];
                            dU(i, j + k) += shapeFunction(polyOrder, coords[order][l])(k) * kernel_laplace2D(node, segment, coords[order][l]) * jac * weights[order][l];
                        }
                    }
                }
            }
        }
    }
    BemVector bp(n);
//    for(int i = 0; i < n; i++ )
//    {
//        for(int j = 0; j < mesh.m_elements.count(); j++)
//        {
//            double R = sqrt((mesh.m_segments[i].gravity().x - mesh.m_elements[j].gravity().x) * (mesh.m_segments[i].gravity().x - mesh.m_elements[j].gravity().x) +
//                            (mesh.m_segments[i].gravity().y - mesh.m_elements[j].gravity().y) * (mesh.m_segments[i].gravity().y - mesh.m_elements[j].gravity().y));
//            bp(i) += 1/(2 * M_PI) * mesh.m_elements[j].f() * log(R) * mesh.m_elements[j].araea();
//        }
//    }

    // Rearranging matrices
    BemMatrix A(n, n);
    BemMatrix C(n,n);
    BemVector rsv(n);

    for (int i = 0; i < n; i++)
    {
        if(mesh.m_nodes[i]->isEssential)
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) =  - dU(j, i);
                C(j, i) =  - dT(j, i);
            }
            rsv(i) = mesh.m_nodes[i]->value;
        }

        else
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) =    dT(j, i);
                C(j, i) =    - dU(j, i);
            }
            rsv(i) = mesh.m_nodes[i]->normalDerivation;
        }
    }

    // ToDo: Poisson - right side vector


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
    domainSolution();

    foreach(Node * node, mesh.m_nodes)
    {
        qDebug() << node->x << node->y;
        qDebug() << node->isEssential;
        qDebug() << node->value;
        qDebug() << node->normalDerivation;
    }
}

void Bem::domainSolution()
{

    for(int i = 0; i < mesh.m_elements.count(); i++)
    {
        Element & element = mesh.m_elements[i];
        for(int i = 0; i < 3; i++)
        {
            element.nodeValues[i] = potential(element.m_nodes.at(i)->x, element.m_nodes.at(i)->y);
        }
    }
}

double Bem::getValue(double x, double y)
{
    double result = 0;
    for(int i = 0; i < mesh.m_elements.count(); i++)
    {
        Node p(x,y);
        if(mesh.m_elements[i].containsPoint(p))
        {
            result = mesh.m_elements[i].value(p);
        }
    }

    return result;
}

double Bem::potential(double x, double y)
{
    // QTime t;
    // t.start();
    // double (Bem::*kernel)(Node, Segment, double) = & Bem::kernel_laplace2D;
    // double (Bem::*kernel_derivation)(Node, Segment, double) = & Bem::kernel_laplace2D_derivation;

    double u = 0;
    Node p(x, y);
    int n = mesh.m_segments.count();
    for(int i = 0; i < n; i++)
    {
        Segment edge = mesh.m_segments[i];
        if(p.distanceOf(edge.firstNode()) < EPS_ZERO)
        {
            // qDebug() << edge.value();
            // qDebug() << edge.firstPoint().toString();
            return edge.firstNode().value;
        }

        if(p.distanceOf(edge.secondNode()) < EPS_ZERO)
        {
            // qDebug() << mesh.m_edgeComponents[(i + 1) % n].value();
            // qDebug() << edge.secondNode().toString();
            return edge.secondNode().value;
        }

        if(edge.distanceOf(p) < 100 * EPS_ZERO)
            return edge.value();

        Node a = edge.firstNode();
        Node b = edge.secondNode();
        double H = 0;
        double G = 0;
        Node integ = integral(p, a, b);
        G = integ.y;
        H = integ.x;
        for(int k = 0; k < 1; k++)
        {
            G = integ.y;
            H = integ.x;
            double delta_u;
            delta_u = - H * edge.value() +  G * edge.derivation();

            u = u + delta_u;
        }
    }

    //    for(int i = 0; i < m; i++)
    //    {
    //        double R = sqrt((p.x - m_elements[i].gravity().x) * (p.x - m_elements[i].gravity().x) - (p.y - m_elements[i].gravity().y) * (p.y - m_elements[i].gravity().y));
    //        u = u - 1 / (2 * M_PI) * m_elements[i].f() * log(R) * m_elements[i].araea();
    //    }

    // qDebug("Time elapsed: %f ms", t.elapsed());
    return u;
}


Node Bem::integral(Node v, Node pa, Node pb)
{
    // center
    Node center = (pa + pb) / 2;

    // shift of center of the edge on the position [0, 0]
    Node at = pa - center;
    Node bt = pb - center;

    // shift of reference point
    Node vt = v - center;

    // rotation
    double phi = - atan2(bt.y, bt.x);
    vt.rotate(phi);
    at.rotate(phi);
    bt.rotate(phi);

    double H;
    double G;

    if(v.distanceOf(center) < EPS_ZERO)
    {
        // singular point
        H = 0;
        double m = abs(bt.x);
        G = -1 / ( M_PI)  * (m * log(m) - m);
        return Node(H, G);
    } else
    {

        double x = 0;
        double y = 0;
        double a = (vt.x - at.x);
        double b = (vt.x - bt.x);

        if(abs(vt.y) < EPS_ZERO)
        {

            if((abs(vt.x - at.x) < EPS_ZERO) && (abs(vt.x - bt.x) < EPS_ZERO))
            {
                assert(0);
            }

            if(abs(b) < EPS_ZERO)
            {
                x = (a > 0) ? M_PI_2 : - M_PI_2;
                y = M_PI_4;
                G =  - ( a * (-2 + log(a * a)))/(4 * M_PI);
            } else
                if(abs(a) < EPS_ZERO)
                {
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    x = M_PI_4;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI);
                }
                else
                {
                    x = (a > 0) ? M_PI_2 : - M_PI_2;
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI) - (a * (-2 + log(a * a)))/(4 * M_PI);
                }
        }
        else
        {
            x = atan( a / vt.y);
            y = atan( b / vt.y);
            G = (2 * vt.y * y + b * (-2 + log(vt.y * vt.y + b * b)))/(4 * M_PI) -
                    (2 * vt.y * x + a * (-2 + log(vt.y * vt.y + a * a)))/(4 * M_PI);
        }
        H = (y - x)  / (2 * M_PI);
    }

    return Node(H, G);
}

template class BemSolution<double>;
