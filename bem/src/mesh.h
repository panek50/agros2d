#ifndef MESH_H
#define MESH_H

#include <QList>
#include <QString>
#include <QObject>

#include "../agros2d-library/meshgenerator.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"
#include "../hermes2d/include/function/exact_solution.h"
#include "../agros2d-library/hermes2d/solver.h"
#include "bem_matrix.h"


struct Node;
class Segment;
class Element;

/*!
 * \brief Class describing mesh of the problem.
 *
 *
 */

class Mesh
{
public:
    void read(MeshSharedPtr mesh);       /*! read the Hermes mesh and build own structures */

    QList<Segment> m_segments;          /*! Surface of the body is divided into segments */
    QList<Element> m_elements;          /*! Volume of the body is divided into elements */
    QList<Node>    m_points;            /*! All used points */
    QList<Node *>  m_integrationPoints; /*! Integration points */
    QList<Node *>  m_innerNodes;        /*! Domain nodes */
    QList<Node *>  m_boundaryNodes;     /*! Boundary nodes */

    int m_nElement;                     /*! Number of elements */
    int m_nSegment;                     /*! Number of segments */
};

struct  Node
{    
    int globalIndex;
    double x, y;
    double real;
    double imag;
    double normalDerivationReal;
    double normalDerivationImag;
    bool isEssential;

    Node() { this->x = 0; this->y = 0; this->real = 0; this->normalDerivationReal = 0; isEssential = false;}
    Node(double x, double y) { this->x = x; this->y = y; this->real = 0; this->normalDerivationReal = 0; isEssential = false;}

    inline Node operator+(const Node & vec) const { return Node(x + vec.x, y + vec.y); }
    inline Node operator-(const Node & vec) const { return Node(x - vec.x, y - vec.y); }
    inline Node operator*(double num) const { return Node(x * num, y * num); }
    inline Node operator/(double num) const { return Node(x / num, y / num); }
    inline double operator&(const Node & vec) const { return x * vec.x + y * vec.y; } // dot product
    inline double operator%(const Node & vec) const { return x * vec.y - y * vec.x; } // cross product
    bool operator!=(const Node &vec) const;
    bool operator==(const Node &vec) const;

    inline double magnitude() const { return sqrt(x * x + y * y); }
    inline double magnitudeSquared() const { return (x * x + y * y); }
    inline double angle() const { return atan2(y, x); }
    inline double distanceOf( Node p) { return ((* this) - p).magnitude(); }
    inline double distanceOfSquared( Node p) { return ((* this) - p).magnitudeSquared(); }
    Node rotate(double angle);

    Node normalizeNode() const
    {
        double m = magnitude();

        double mx = x / m;
        double my = y / m;

        return Node(mx, my);
    }

    QString toString() const
    {
        return QString("x = %1, y = %2, magnitude = %3").
                arg(x).
                arg(y).
                arg(magnitude());
    }
};


class Element
{
public:
    Element (QList<Node *> nodes);
    int id() const { return m_id; }
    Node gravity() { return m_gravity; }
    double araea() { return m_area;}
    void setArea(double area) { m_area = area; }
    void setValue(double value) { m_value = value; }
    double value() { return m_value; }
    double value(Node p);
    double f() { return m_f; }
    double nodeValues[3];
    bool containsPoint(Node x);    
    QList<Node *> m_nodes;
    Mesh * m_mesh;

private:    

    int m_id;
    Node m_gravity;
    double m_area;
    double m_f;
    double m_value;
    double m_derivation;
};


class Segment
{

public:
    Segment(Node *firstNode, Node *lastNode);
    double elementArea() {return m_element->area;}
    int id() { return m_id; }
    Node & firstNode()  { return * m_nodes.first(); }
    Node & lastNode() { return * m_nodes.last(); }
    Node gravity() const { return m_gravity;}
    bool isLyingPoint (Node node);
    double distanceOf(Node point);

    void setElement(Hermes::Hermes2D::Element * element) { m_element = element; }
    Hermes::Hermes2D::Element * element() { return m_element; }

    void setValue(const double value) { m_value = value; }
    double value() { return m_value; }

    void setDerivation(const double derivation) { m_derivation = derivation; }
    double derivation() { return m_derivation; }

    void setEdgeId (int id) { m_id = id; }
    int edgeId() { return m_id; }

    void setEssential(bool isEssential) { m_isEssential = isEssential; }
    bool isEssential() { return m_isEssential; }

    void setApproximationOrder(int order) { m_approximationOrder = order; }
    int approximationOrder() { return m_approximationOrder; }

    void setGeometricOrder(int order) { m_geometricOrder = order; }
    int geometricOrder() { return m_geometricOrder; }

    double length() const { return m_length; }

    double parametricCoordinate(Node node);

    Mesh * m_mesh;
    QList<Node *> m_nodes;
    QList<Node *> m_points;
    double m_jacobian;
    double m_logJacobian;


private:
    int m_id;
    int m_order;
    int m_edgeID;
    int m_approximationOrder;
    int m_geometricOrder;
    bool m_isEssential;

    double m_value;
    double m_derivation;
    double m_length;
    double jacobian;

    QString m_type;
    QString m_boundaryType;    

    Node m_gravity;    

    Hermes::Hermes2D::Element * m_element;    
};


#endif // MESH_H
