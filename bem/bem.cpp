#include <QTextStream>
#include <QTime>

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

Element::Element(Point a, Point b, Point c)
{
    m_points.append(a);
    m_points.append(b);
    m_points.append(c);
    m_gravity = (a + b + c) / 3;
}

Element::Element(QList<Point> points)
{
    m_points = points;
    int i = 0;
    foreach (Point point, m_points) {
        m_gravity = m_gravity + point;
        i++;
    }
    m_gravity = m_gravity / i;

    m_area = 0;
    m_f = 0;
}

EdgeComponent::EdgeComponent(Point firstPoint, Point lastPoint)
{
    m_nodes.append(firstPoint);
    m_nodes.append(lastPoint);
    m_gravity.x = (firstPoint.x + lastPoint.x) / 2;
    m_gravity.y = (firstPoint.y + lastPoint.y) / 2;
    m_length = sqrt(pow(firstPoint.x - lastPoint.x, 2) + pow(firstPoint.y - lastPoint.y, 2));
    m_derivation = 0;
    m_value = 0;
}

bool EdgeComponent::isLyingPoint(Point point)
{

    double dx = m_nodes.last().x - m_nodes.first().x;
    double dy = m_nodes.last().y - m_nodes.first().y;

    Point sp = m_nodes.first();

    double t = ((point.x - sp.x)*dx + (point.y - sp.y)*dy);

    if (t < 0.0)
        t = 0.0;
    else if (t > (dx * dx + dy * dy))
        t = 1.0;
    else
        t /= (dx * dx + dy * dy);

    Point p(sp.x + t*dx, sp.y + t*dy);

    Point dv = point - p;
    return (dv.magnitudeSquared() < EPS_ZERO);
}


Bem::Bem(FieldInfo * const fieldInfo, MeshSharedPtr mesh)
{    
    m_fieldInfo = fieldInfo;
    m_mesh = mesh;
}


void Bem::addPhysics()
{

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
                QList<Point> points;
                for (unsigned i = 0; i < e->get_nvert(); i++)
                {
                    Point point(e->vn[i]->x, e->vn[i]->y);
                    points.append(point);
                    if(e->vn[i]->bnd == 1)
                    {
                        if(atoi(m_mesh->get_boundary_markers_conversion().get_user_marker(m_mesh->get_base_edge_node(e, i)->marker).marker.c_str()) == j)
                        {
                            Point firstPoint(e->vn[i]->x, e->vn[i]->y);
                            Point secondPoint(e->vn[e->next_vert(i)]->x, e->vn[e->next_vert(i)]->y);
                            EdgeComponent component(firstPoint, secondPoint);
                            component.m_element = e;
                            component.m_value = value;
                            component.m_isEssential = isEssential;
                            component.m_edgeID = j;
                            m_edgeComponents.append(component);
                        }
                    }
                }
                Element element(points);
                element.setArea(e->get_area());

                m_elements.append(element);
                // ToDo: read material properties

            }
            // ToDo: improve
            m_nElement = nElement;
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
        output += "First Point:";
        output += " ";
        output += QString::number(component.firstPoint().x);
        output += " ";
        output += QString::number(component.firstPoint().y);
        output += "\n";
        output += "Second Point:";
        output += " ";
        output += QString::number(component.lastPoint().x);
        output += " ";
        output += QString::number(component.lastPoint().y);
        output += "\n";

    }
    return output;
}


Point Bem::globalCoordinates(int polyOrder, double xi, EdgeComponent segment)
{
    Point v;
    BemVector Sf(polyOrder + 1);
    int n = segment.m_nodes.count();
    Sf = shapeFunction(n - 1, xi);

    for(int i = 0; i < n; i++)
    {
        v = v + segment.m_nodes.at(i) * Sf(i);
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
    if(polyOrder == 1)
        return Dn;
    Dn(2) = -2.0 * xi;
    Dn(0) = Dn(0) - 0.5 * Dn(2);
    Dn(1) = Dn(1) - 0.5 * Dn(2);
    return Dn;
}

double Bem::jacobian(int polyOrder, double xi, EdgeComponent segment)
{
    double dGamma = 0;
    if(polyOrder == 1)
    {
        double x = (segment.lastPoint().x - segment.firstPoint().x);
        double y = (segment.lastPoint().y - segment.firstPoint().y);

        dGamma = 0.5 * sqrt(x*x + y*y);
    }

    return dGamma;
}

Point Bem::normalVector(int polyOrder, double xi, EdgeComponent segment)
{
    BemVector Dn = shapeFunctionDerivative(polyOrder, xi);
    Point dVec;
    Point n;
    for(int i = 0; i < polyOrder + 1; i++)
    {
        dVec = dVec  + segment.m_nodes.at(i) * Dn(i);
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
    return m_bem->potential(x, y);
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

double Bem::kernel_length(Point refNode, EdgeComponent segment, double xi)
{
    return 1;
}


double Bem::kernel_laplace2D(Point refNode, EdgeComponent segment, double xi)
{
    double r = globalCoordinates(2, xi, segment).distanceOf(refNode);
    return 1.0 / (2 * M_PI) * log(1/r);
}


double Bem::kernel_laplace2D_derivation(Point refNode, EdgeComponent segment, double xi)
{
    Point r = globalCoordinates(1, xi, segment) - refNode;
    Point n = normalVector(1, xi, segment);
    double rNorm = r.magnitude();
    return - 1.0 / (2 * M_PI) * (n.x *  r.x  + n.y * r.y) / (rNorm * rNorm);
}

void Bem::solve()
{        
    int n = m_edgeComponents.count();
    double (Bem::*kernel)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D;
    double (Bem::*kernel_derivation)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D_derivation;

    BemMatrix H(n, n);
    BemMatrix G(n, n);
    for (int i = 0; i < n; i++)
    {
        Point v(m_edgeComponents.at(i).gravity());
        for (int j = 0; j < n; j++)
        {
            double length = 0;
            double normDeriv = 0;
            Point a = m_edgeComponents[j].firstPoint();
            Point b = m_edgeComponents[j].lastPoint();
            // Point integ = integral(v, a, b);
            if( i == j)
            {
                length = gaussLaguerre(3, v, m_edgeComponents[j], kernel);
                normDeriv = gaussLaguerre(3, v, m_edgeComponents[j], kernel_derivation);
            }
            else
            {
                length = quad(3, v, m_edgeComponents[j], kernel);
                normDeriv = quad(3, v, m_edgeComponents[j], kernel_derivation);
            }
            // qDebug() << normDeriv << "  " << integ.x;

            H(i, j) = normDeriv;
            G(i, j) = length;

            if (i == j)
                H(i, j) = H(i, j) + 0.5;
        }
    }

    BemVector bp(n);
    for(int i = 0; i < n; i++ )
    {
        for(int j = 0; j < m_elements.count(); j++)
        {
            double R = sqrt((m_edgeComponents[i].gravity().x - m_elements[j].gravity().x) * (m_edgeComponents[i].gravity().x - m_elements[j].gravity().x) +
                            (m_edgeComponents[i].gravity().y - m_elements[j].gravity().y) * (m_edgeComponents[i].gravity().y - m_elements[j].gravity().y));
            bp(i) += 1/(2 * M_PI) * m_elements[j].f() * log(R) * m_elements[j].araea();
        }
    }

    // Rearranging matrices
    BemMatrix A(n, n);
    BemMatrix C(n,n);
    BemVector rsv(n);

    for (int i = 0; i < n; i++)
    {
        if(m_edgeComponents[i].m_isEssential)
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) = - G(j, i);
                C(j, i) = - H(j, i);
            }
            rsv(i) = m_edgeComponents[i].m_value;
        }

        else
        {
            for(int j = 0; j < n; j++)
            {
                A(j, i) =   H(j, i);
                C(j, i) = - G(j, i);
            }
            rsv(i) = m_edgeComponents[i].m_derivation;
        }
    }

    // ToDo: Poisson - right side vector


    BemVector b(n);
    b = C * rsv + bp;


    BemVector results(n);
    results = A.solve(b);



    for (int i = 0; i < n; i++)
    {
        if(m_edgeComponents[i].m_isEssential)
        {
            m_edgeComponents[i].m_derivation = results(i);
        }
        else
        {
            m_edgeComponents[i].m_value = results(i);
        }
    }
}

double Bem::potential(double x, double y)
{    
    // QTime t;
    // t.start();

    double (Bem::*kernel)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D;
    double (Bem::*kernel_derivation)(Point, EdgeComponent, double) = &Bem::kernel_laplace2D_derivation;

    double u = 0;
    Point p(x, y);
    int n = m_edgeComponents.count();
    for(int i = 0; i < n; i++)
    {

        EdgeComponent edge = m_edgeComponents[i];
        Point a = edge.firstPoint();
        Point b = edge.lastPoint();
        double length = 0;
        double normDeriv = 0;

        Point integ = integral(p, a, b);
        length = quad(7, p, m_edgeComponents[i], kernel);
        normDeriv = quad(7, p, m_edgeComponents[i], kernel_derivation);

        if(abs(integ.x - normDeriv) > 1e-3)
            if(!edge.isLyingPoint(p));
                qDebug() << p.toString();

        double G = length;
        double H = normDeriv;
        double delta_u;
        delta_u = - H * edge.m_value +  G * edge.m_derivation;

        u = u + delta_u;
    }

    //    for(int i = 0; i < m; i++)
    //    {
    //        double R = sqrt((p.x - m_elements[i].gravity().x) * (p.x - m_elements[i].gravity().x) - (p.y - m_elements[i].gravity().y) * (p.y - m_elements[i].gravity().y));
    //        u = u - 1 / (2 * M_PI) * m_elements[i].f() * log(R) * m_elements[i].araea();
    //    }

    // qDebug("Time elapsed: %f ms", t.elapsed());
    return u;
}



double Bem::quad(int order, Point refNode, EdgeComponent segment, double (Bem::*kernel)(Point, EdgeComponent, double))
{
    // Gauss coordinates and weights
    QList<QList<double> > weights;
    QList<QList<double> > coords;


    QList<double> wi;
    QList<double> coord;
    QList<double> rads;

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

    double integral = 0;


    for(int i = 0; i <= order; i++)
    {
        double  jac =jacobian(1, 0,segment);
        integral += (this->*kernel)(refNode, segment, coords[order][i]) * jac * weights[order][i];
    }

    return integral;
}

double Bem::gaussLaguerre(int order, Point refNode, EdgeComponent segment, double (Bem::*kernel)(Point, EdgeComponent, double))
{
    // Gauss-Laguerre coordinates and weights
    QList<QList<double> > weights;
    QList<QList<double> > coords;


    QList<double> wi;
    QList<double> coord;
    QList<double> rads;

    wi.append(1);
    coord.append(0.5);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.112008806);
    coord.append(0.602276908);

    wi.append(0.718539319);
    wi.append(0.281460680);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.063890793);
    coord.append(0.368997063);
    coord.append(0.766880303);

    wi.append(0.513404552);
    wi.append(0.391980041);
    wi.append(0.0946154065);
    weights.append(wi);
    coords.append(coord);


    wi.clear();
    coord.clear();
    coord.append(0.0414484801);
    coord.append(0.245274914);
    coord.append(0.556165453);
    coord.append(0.848982394);

    wi.append(0.383464068);
    wi.append(0.386875317);
    wi.append(0.190435126);
    wi.append(0.0392254871);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.0291344721);
    coord.append(0.173977213);
    coord.append(0.411702520);
    coord.append(0.677314174);
    coord.append(0.894771361);

    wi.append(0.297893471);
    wi.append(0.349776226);
    wi.append(0.234488290);
    wi.append(0.0989304595);
    wi.append(0.0189115521);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.021634005);
    coord.append(0.129583391);
    coord.append(0.314020449);
    coord.append(0.538657217);
    coord.append(0.756915337);
    coord.append(0.922668851);

    wi.append(0.238763662);
    wi.append(0.308286573);
    wi.append(0.245317426);
    wi.append(0.142008756);
    wi.append(0.055454622);
    wi.append(0.0101689586);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.016719355);
    coord.append(0.100185677);
    coord.append(0.246294246);
    coord.append(0.433463493);
    coord.append(0.632350988);
    coord.append(0.811118626);
    coord.append(0.940848166);
    wi.append(0.196169389);
    wi.append(0.270302644);
    wi.append(0.239681873);
    wi.append(0.165775774);
    wi.append(0.088943227);
    wi.append(0.0331943043);
    wi.append(0.0059327870);
    weights.append(wi);
    coords.append(coord);

    wi.clear();
    coord.clear();
    coord.append(0.0133202441);
    coord.append(0.0797504290);
    coord.append(0.197871029);
    coord.append(0.354153994);
    coord.append(0.529458575);
    coord.append(0.701814529);
    coord.append(0.849379320);
    coord.append(0.953326450);
    wi.append(0.164416604);
    wi.append(0.237525610);
    wi.append(0.226841984);
    wi.append(0.175754079);
    wi.append(0.112924030);
    wi.append(0.0578722107);
    wi.append(0.0209790737);
    wi.append(0.00368640710);
    weights.append(wi);
    coords.append(coord);

    double integral = 0;


    for(int i = 0; i <= order; i++)
    {
        double  jac =jacobian(1, 0, segment);
            integral += (this->*kernel)(refNode, segment, globalCoordinates(order, coords[order][i], segment).distanceOf(refNode)) * jac * weights[order][i];
    }
    return integral;
}

Point Bem::integral(Point v, Point pa, Point pb)
{    
    // center
    Point center = (pa + pb) / 2;

    // shift of center of the edge on the position [0, 0]
    Point at = pa - center;
    Point bt = pb - center;

    // shift of reference point
    Point vt = v - center;

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
        return Point(H, G);
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

    return Point(H, G);
}

template class BemSolution<double>;
