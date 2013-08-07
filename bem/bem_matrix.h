#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>

#include <exception>
#include<QObject>
#include<QString>

class BemException : public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Matrices are not aligned.";
    }
};

class BemMatrix : public QObject
{
    Q_OBJECT
public:
    BemMatrix() { assert(0); }
    BemMatrix(int m_row, int m_column);
    BemMatrix(const BemMatrix  &m);
    BemMatrix(double array[], int row, int column);
    ~BemMatrix();

    BemMatrix & operator=(const BemMatrix & m);
    bool operator== (const BemMatrix & m) const;
    BemMatrix operator+(const BemMatrix & m) const;
    BemMatrix operator-(const BemMatrix & m);
    BemMatrix operator*(BemMatrix & m);
    BemMatrix operator*(const double & x) const;
    BemMatrix operator/(double x) {return operator *(1/x); }
    BemMatrix solve(const BemMatrix & m);
    friend BemMatrix operator*(double x, BemMatrix m);
    friend BemMatrix operator/(double x, BemMatrix m);

    int row() const { return m_row; }
    void setRow(int row) { m_row = row; }
    int column() const { return m_column; }
    void setColumn(int column) { m_column = column; }

    double & operator() (unsigned row, unsigned col);
    void set(int m_row, int m_column, double value);

    QString toString();
    int countColumn();
    int countRow();
    void eye();
    void clear();
    void gaussElim();
    void luFactorisation();

protected:
    void createArray(int m_row, int m_column);
    double *m_array;
    int m_row;
    int m_column;

};

enum BemVectorType
{
    VectorTypeColumn = 0,
    VectorTypeRow = 1
};

class BemVector : public BemMatrix
{
public:    
    BemVector(int n) : BemMatrix(n, 1) { m_type = VectorTypeColumn; }
    BemVector(BemMatrix m) : BemMatrix(m) { if(column() != 1) throw BemException(); }
    double & operator() (int i) { return BemMatrix::operator()(0, i); }
    double length() const;
    BemVector & transpose();

private:
    BemVectorType m_type;
};


class Node : public BemVector
{
public:
    Node() : BemVector(2) {}
    Node(BemMatrix m) : BemVector(m) {}
    Node(double x, double y) : BemVector(2) { m_array[0] = x; m_array[1] = y; }
    Node & operator=(BemVector & v);
    Node & rotate(double phi);
    double distanceOf(const Node & node);
    int id() { return m_id; }
    void set_id(int id) { m_id = id; }

private:
    int m_id;
};


Q_DECLARE_METATYPE(BemMatrix);

#endif // MATRIX_H


