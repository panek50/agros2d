#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>

#include<QObject>
#include<QString>

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
    BemMatrix operator*(const BemMatrix & m);
    BemMatrix operator*(const double & x) const;
    BemMatrix operator/(double x) {return operator *(1/x); }
    friend BemMatrix operator*(double x, BemMatrix m);
    friend BemMatrix operator/(double x, BemMatrix m);

    int row() const { return m_row; }
    int column() const { return m_column; }
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
    double *m_array;

private:
    void createArray(int m_row, int m_column);
    int m_row;
    int m_column;

};

class BemVector : public BemMatrix
{
public:
    BemVector(BemMatrix m);
    BemVector(const BemVector & v);
    BemVector(int n);
    BemVector & operator=(const BemVector & m);
    double length();
    void set(int n, double value) {  BemMatrix::set(0, n, value); }
    BemVector operator +( const BemVector & v) { return static_cast<BemVector>(BemMatrix::operator +(v)); }
    BemVector operator -( const BemVector & v) { return static_cast<BemVector>(BemMatrix::operator -(v)); }
    BemVector operator *( const double & x) { return static_cast<BemVector>(BemMatrix::operator *(x)); }
    BemVector operator /( const double & x) { return static_cast<BemVector>(BemMatrix::operator /(x)); }
    double & operator()(int n) { return BemMatrix::operator() (0, n); }
    double operator*(BemVector v);

    // BemVector operator*(BemMatrix  m);
    friend BemVector operator*(BemMatrix m, BemVector v);
};


class  Node: public BemVector
{
public:
    Node(BemVector v);
    Node() : BemVector(2) {}
    Node(double x, double y);
    Node & operator=(const Node & m);
    Node operator +( const Node & v) { return static_cast<Node>(BemVector::operator +(v)); }
    Node operator -( const Node & v) { return static_cast<Node>(BemVector::operator -(v)); }
    Node operator *( const double & x) { return static_cast<Node>(BemVector::operator *(x)); }
    Node operator /( const double & x) { return static_cast<Node>(BemVector::operator /(x)); }
    void shift(double x, double y);
    void rotate(double angle);
    int id;
};

Q_DECLARE_METATYPE(BemMatrix);

#endif // MATRIX_H


