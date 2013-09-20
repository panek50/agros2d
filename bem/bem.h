#ifndef BEM_H
#define BEM_H

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

class Element
{
public:
    Element (QList<Point> nodes);
    Element (Point a, Point b, Point c);
    int id() const { return m_id; }
    Point gravity() { return m_gravity; }
    double araea() { return m_area;}
    void setArea(double area) { m_area = area; }
    double f() { return m_f; }

private:
    QList<Point> m_points;
    int m_id;
    Point m_gravity;
    double m_area;
    double m_f;
};


class EdgeComponent
{

public:
    EdgeComponent(Point firstPoint, Point lastPoint);
    double elementArea() {return m_element->area;}
    int id() { return m_id; }

    // ToDo: make it private
    Hermes::Hermes2D::Element * m_element;
    Point & firstPoint()  { return m_nodes.first(); }
    Point & lastPoint() { return m_nodes.last(); }
    Point gravity() const { return m_gravity;}
    bool isLyingPoint (Point node);
    QString boundaryType;
    double m_value;
    double m_derivation;
    bool m_isEssential;
    int m_edgeID;
    QString m_type;
    double m_length;
    QList<Point> m_nodes;

private:
    int m_id;    
    Point m_gravity;
};


class Bem
{

public:
    Bem(FieldInfo* field, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh> mesh);
    ~Bem() {}
    void readMesh();
    void addPhysics();
    void assemblyMatrices();
    void solve();
    double quad(int oder, Point refNode, EdgeComponent segment,double (Bem::*kernel)(Point, EdgeComponent, double));
    double gaussLaguerre(int oder, Point refNode, EdgeComponent segment,double (Bem::*kernel)(Point, EdgeComponent, double));

    double kernel_length(Point refNode, EdgeComponent segment, double xi);
    double kernel_laplace2D_derivation(Point refNode, EdgeComponent segment, double xi);
    double kernel_laplace2D(Point refNode, EdgeComponent segment, double xi);

    QString toString();
    Point integral(Point v, Point a, Point b);
    double potential(double x, double y);
    Point globalCoordinates(int polyOrder, double xi, EdgeComponent segment);
    BemVector shapeFunction(int polyOrder, double xi);
    BemVector shapeFunctionDerivative(int polyOrder, double xi);
    Point normalVector(int polyOrder, double xi, EdgeComponent segment);
    double jacobian(int polyOrder, double xi, EdgeComponent segment);

    // Todo: make it private
    QList<EdgeComponent> m_edgeComponents;
    QList<Element> m_elements;
    int m_nElement;


private:
    MeshSharedPtr m_mesh;
    FieldInfo * m_fieldInfo;
};

template <typename Scalar>
class BemSolution: public Hermes::Hermes2D::ExactSolutionScalar<Scalar>
{
public:
    BemSolution<Scalar>(MeshSharedPtr mesh);
    virtual ~BemSolution() {}
    virtual Scalar value (double x, double y) const;
    virtual void derivatives (double x, double y, Scalar& dx, Scalar& dy) const;
    virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const { return Hermes::Ord(20); }
    virtual Hermes::Hermes2D::MeshFunction<Scalar>* clone() const;    
    /// Saves the exact solution to an XML file.
    void save(const char* filename) const { }    
    void setSolver(Bem * bem) { m_bem = bem; }    

protected:    
    // virtual void precalculate(int order, int mask) { qDebug() << "OK";}
    double constant;
    Bem * m_bem;
};




#endif // BEM_H
