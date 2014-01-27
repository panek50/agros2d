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

template <typename Type>
class BemMatrix;

template <typename Type>
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

template <typename Type>
class BemMatrix
{
public:
    BemMatrix() { assert(0); }
    BemMatrix(int m_row, int m_column);
    BemMatrix(const BemMatrix  &m);
    // BemMatrix(BemVector & v);
    BemMatrix(Type array[], int row, int column);
    ~BemMatrix();

    BemMatrix & operator=(const BemMatrix & m);
    bool operator== (const BemMatrix & m) const;
    BemMatrix<Type> operator+(const BemMatrix & m) const;
    BemMatrix<Type> operator-(const BemMatrix & m);
    BemMatrix<Type> operator*(BemMatrix & m);
    BemMatrix<Type> operator*(const Type & x) const;
    BemVector<Type> operator*(BemVector<Type> & v);
    // BemMatrix<Type> operator/(Type x) {return operator * (1 / x); }
    BemVector<Type> solve(BemVector<Type> & v);
    friend BemMatrix<Type> operator*(Type x, BemMatrix<Type> m);
    // friend BemVector operator*(BemVector v, BemMatrix m);
    friend BemMatrix<Type> operator/(Type x, BemMatrix<Type> m);

    int row() const { return m_row; }
    void setRow(int row) { m_row = row; }
    int column() const { return m_column; }
    void setColumn(int column) { m_column = column; }

    Type & operator() (unsigned row, unsigned col) const;
    void set(int m_row, int m_column, Type value);

    QString toString();
    int countColumn();
    int countRow();
    void eye();
    void clear();
    void gaussElim();
    void luFactorisation();
    Type *m_array;

protected:
    void createArray(int m_row, int m_column);    
    int m_row;
    int m_column;
};

template <typename Type>
class BemVector
{
public:
    BemVector() { assert(0); }
    BemVector(int n);
    BemVector(const BemVector & v);
    ~BemVector() { delete[] m_array; }

    Type & operator() (int i);
    void operator+=(BemVector & v);
    BemVector operator +(const BemVector<Type> & m) const;
    BemVector operator * (const Type & x) const;
    BemVector operator * (BemMatrix<Type> & m) const;
    BemVector<Type> & operator=(const BemVector<Type> & v);

    QString toString();
    int n() const { return m_n; }
    // double length() const;
    // friend BemVector<Type> BemMatrix<Type>::solve(BemVector<Type> & v);
    void clear();
    Type * m_array;

private:  
    int m_n;
};

Q_DECLARE_METATYPE(BemMatrix<double>)
Q_DECLARE_METATYPE(BemVector<double>)
#endif // MATRIX_H


