#ifndef MATRIX_H
#define MATRIX_H


#include <assert.h>
#include <exception>
#include <complex>

#include <QObject>
#include <QString>

#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <lapacke.h>
#include <cblas.h>

class BemMatrix;
class BemVector;
class BemComplexVector;
class BemComplexMatrix;


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


class BemComplexMatrix
{
public:
    BemComplexMatrix() { assert(0); }
    BemComplexMatrix(int m_row, int m_column);
    BemComplexMatrix(const BemComplexMatrix  &m);
    // BemMatrix(BemVector & v);
    BemComplexMatrix(double array[], int row, int column);
    ~BemComplexMatrix();

    BemComplexMatrix & operator=(const BemComplexMatrix & m);
    bool operator== (const BemComplexMatrix & m) const;
    BemComplexMatrix operator+(const BemComplexMatrix & m) const;
    BemComplexMatrix operator-(const BemComplexMatrix & m);
    BemComplexMatrix operator*(BemComplexMatrix & m);
    BemComplexMatrix operator*(const double & x) const;
    BemComplexVector operator*(BemComplexVector & v);
    BemComplexMatrix operator/(double x) {return operator *(1/x); }
    BemComplexVector solve(BemComplexVector &v);
    friend BemComplexMatrix operator*(double x, BemComplexMatrix m);
    // friend BemVector operator*(BemVector v, BemMatrix m);
    friend BemComplexMatrix operator/(double x, BemComplexMatrix m);

    int row() const { return m_row; }
    void setRow(int row) { m_row = row; }
    int column() const { return m_column; }
    void setColumn(int column) { m_column = column; }


    std::complex<double> & operator() (unsigned row, unsigned col);

    void set(int m_row, int m_column, std::complex<double> number);
    QString toString();
    int countColumn();
    int countRow();
    void eye();
    void clear();
    void gaussElim();
    void luFactorisation();
    std::complex<double> *m_array;

protected:
    void createArray(int m_row, int m_column);
    int m_row;
    int m_column;
};

class BemComplexVector
{
public:
    BemComplexVector() { assert(0); }
    BemComplexVector(int n);
    BemComplexVector(const BemComplexVector & v);
    ~BemComplexVector() { delete[] m_array; }

    std::complex<double> & operator() (int i);
    void operator+=(BemComplexVector & v);
    BemComplexVector operator +(const BemComplexVector & v) const;
    BemComplexVector operator * (const double & x) const;
    BemComplexVector operator * (BemComplexMatrix & m) const;
    BemComplexVector & operator=(const BemComplexVector & v);

    QString toString();
    int n() const { return m_n; }
    double length() const;
    friend BemComplexVector BemComplexMatrix::solve(BemComplexVector & v);
    void clear();
   std::complex<double> * m_array;

private:    
    int m_n;
};

Q_DECLARE_METATYPE(BemMatrix)
Q_DECLARE_METATYPE(BemVector)
Q_DECLARE_METATYPE(BemComplexMatrix)
Q_DECLARE_METATYPE(BemComplexVector)
#endif // MATRIX_H


