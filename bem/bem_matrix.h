#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>

#include <exception>
#include<QObject>
#include<QString>

class BemMatrix;
class BemVector;


enum BemVectorType
{
    VectorTypeColumn = 0,
    VectorTypeRow = 1
};

class BemException : public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Matrices are not aligned.";
    }
};


class BemMatrix
{
public:
    BemMatrix() { assert(0); }
    BemMatrix(int m_row, int m_column);
    BemMatrix(const BemMatrix  &m);
    // BemMatrix(BemVector & v);
    BemMatrix(double array[], int row, int column);
    ~BemMatrix();

    BemMatrix & operator=(const BemMatrix & m);
    bool operator== (const BemMatrix & m) const;
    BemMatrix operator+(const BemMatrix & m) const;
    BemMatrix operator-(const BemMatrix & m);
    BemMatrix operator*(BemMatrix & m);
    BemMatrix operator*(const double & x) const;
    BemVector operator*(BemVector & v);
    BemMatrix operator/(double x) {return operator *(1/x); }
    BemVector solve(BemVector & v);
    friend BemMatrix operator*(double x, BemMatrix m);
    // friend BemVector operator*(BemVector v, BemMatrix m);
    friend BemMatrix operator/(double x, BemMatrix m);

    int row() const { return m_row; }
    void setRow(int row) { m_row = row; }
    int column() const { return m_column; }
    void setColumn(int column) { m_column = column; }

    double & operator() (unsigned row, unsigned col) const;
    void set(int m_row, int m_column, double value);

    QString toString();
    int countColumn();
    int countRow();
    void eye();
    void clear();
    void gaussElim();
    void luFactorisation();
    double *m_array;

protected:
    void createArray(int m_row, int m_column);    
    int m_row;
    int m_column;
};

class BemVector
{
public:
    BemVector() { assert(0); }
    BemVector(int n);
    BemVector(const BemVector & v);
    ~BemVector() { delete[] m_array; }

    double & operator() (int i);
    void operator+=(BemVector & v);
    BemVector operator +(const BemVector & m) const;
    BemVector operator * (const double & x) const;
    BemVector operator * (BemMatrix & m) const;
    BemVector & operator=(const BemVector & v);

    QString toString();
    int n() const { return m_n; }
    double length() const;
    friend BemVector BemMatrix::solve(BemVector & v);
    void clear();

private:
    double * m_array;
    int m_n;
};

Q_DECLARE_METATYPE(BemMatrix);

#endif // MATRIX_H


