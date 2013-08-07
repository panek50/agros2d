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


class EdgeComponent
{

public:
    EdgeComponent(Node firstNode, Node secondNode);
    double elementArea() {return m_element->area;}
    int id() { return m_id; }

    // ToDo: make it private
    Hermes::Hermes2D::Element * m_element;
    Node firstNode() { return m_firstNode; }
    Node secondNode() { return m_secondNode; }    
    Node gravity() const { return m_gravity;}
    QString boundaryType;
    double m_value;
    double m_derivation;
    bool m_isEssential;
    int m_edgeID;
    QString m_type;
    double m_length;

private:
    int m_id;
    Node m_firstNode;
    Node m_secondNode;
    Node m_gravity;
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
    Node integral(Node v, Node a, Node b);
    QList<EdgeComponent> m_edgeComponents;
    int m_nElement;

private:
    BemSolution<double> * m_solution;
    MeshSharedPtr m_mesh;
    FieldInfo * m_fieldInfo;    

};

#endif // BEM_H
