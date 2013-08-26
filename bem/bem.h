#ifndef BEM_H
#define BEM_H
#define SOLVE_BEM

#include <QList>
#include <QString>

#include "util.h"
#include "../agros2d-library/meshgenerator.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"

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


class Bem
{

public:
    Bem(FieldInfo * const fieldInfo);
    void readMesh();
    void addPhysics();
    void assemblyMatrices();
    void solve();
    QString toString();

private:
    QList<Node> m_nodeList;
    QList<Edge> m_edgeList;
    QList<Edge> m_bounderyList;
    QList<Element> m_elementList;
    FieldInfo * m_fieldInfo;
};

#endif // BEM_H
