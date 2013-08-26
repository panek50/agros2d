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

EdgeComponent::EdgeComponent(Point firstPoint, Point secondPoint)
{
    m_firstPoint = firstPoint;
    m_secondPoint = secondPoint;
    m_gravity.x = (firstPoint.x + secondPoint.x) / 2;
    m_gravity.y = (firstPoint.y + secondPoint.y) / 2;
    m_length = sqrt(pow(firstPoint.x - secondPoint.x, 2) + pow(firstPoint.y - secondPoint.y, 2));
    m_derivation = 0;
    m_value = 0;
}

bool EdgeComponent::isLyingPoint(Point point)
{

    double dx = m_secondPoint.x - m_firstPoint.x;
    double dy = m_secondPoint.y - m_firstPoint.y;

    Point sp = m_firstPoint;

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
        output += QString::number(component.secondPoint().x);
        output += " ";
        output += QString::number(component.secondPoint().y);
        output += "\n";

    }
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

void Bem::solve()
{        
    int n = m_edgeComponents.count();
    BemMatrix H(n, n);
    BemMatrix G(n, n);
    for (int i = 0; i < n; i++)
    {
        Point v(m_edgeComponents.at(i).gravity());
        for (int j = 0; j < n; j++)
        {
            Point a = m_edgeComponents[j].firstPoint();
            Point b = m_edgeComponents[j].secondPoint();
            Point integ = integral(v, a, b);

            H(i, j) = integ.x;
            G(i, j) = integ.y;

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


    BemVector b = C * rsv + bp;
    BemVector results = A.solve(b);

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
    double u = 0;
    Point p(x, y);
    int n = m_edgeComponents.count();
    for(int i = 0; i < n; i++)
    {

        EdgeComponent edge = m_edgeComponents[i];
        Point a = edge.firstPoint();
        Point b = edge.secondPoint();

        if ((p.distanceOf(edge.gravity())) < EPS_ZERO)
        {
            return edge.m_value;
        }

        Point integ = integral(p, a, b);

        double H = integ.x;
        double G = integ.y;
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
