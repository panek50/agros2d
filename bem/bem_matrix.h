#ifndef MATRIX_H
#define MATRIX_H

#include<QObject>
#include<QString>
using namespace std;

class Matrix : public QObject
{
    Q_OBJECT
public:    
    Matrix () {}
    Matrix(int m_row, int m_column);
    Matrix(const Matrix  &m);
    Matrix(double array[], int row, int column);
    ~Matrix();

    Matrix & operator=(const  Matrix & m);
    bool operator== (const Matrix & m) const;
    Matrix operator+(const Matrix & m);
    Matrix operator-(const Matrix & m);
    Matrix operator*(const Matrix & m);

    int row() { return m_row; }
    int colum() { return m_column; }   
    double & operator() (unsigned row, unsigned col);
    void set(int m_row, int m_column, double value);

    QString toString();
    int countColumn();
    int countRow();
    void eye();
    void clear();
    void gaussElim();
    void luFactorisation();

private:
    void createArray(int m_row, int m_column);
    int m_row;
    int m_column;
    double *m_array;
};

class Vector : public Matrix
{
public:
    Vector(int n);
    void set(int n, double value) {  Matrix::set(0, n, value); }
    double operator()(int n) { return Matrix::operator() (0, n); }
};

Q_DECLARE_METATYPE(Matrix);

#endif // MATRIX_H
