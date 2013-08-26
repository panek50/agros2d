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


enum BoundaryConditionType
{
    BoundaryConditionType_None,
    BoundaryConditionType_Essential,
    BoundaryConditionType_Vector
};

struct Node
{
    int m_index;
    double m_x;
    double m_y;
};

struct Element
{
    int m_index;
    int m_nodes[3];
    int m_marker;
};

struct Edge
{
    int m_index;
    int m_nodes[2];
    int m_marker;
    bool m_boundary;
    BoundaryConditionType m_boundary_type;
    double m_boundary_value;
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

protected:

    // virtual void precalculate(int order, int mask) { qDebug() << "OK";}
    double constant;
};


class Bem
{

public:
    Bem(FieldInfo* field, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh> mesh);
    void readMesh();
    void addPhysics();
    void assemblyMatrices();
    void solve();
    QString toString();
    BemSolution<double> * getSolution();

private:
    BemSolution<double> * m_solution;
    MeshSharedPtr m_mesh;
    QList<Node> m_nodeList;
    QList<Edge> m_edgeList;
    QList<Edge> m_bounderyList;
    QList<Element> m_elementList;
    FieldInfo * m_fieldInfo;
};

#endif // BEM_H
