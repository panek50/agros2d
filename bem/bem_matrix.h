#ifndef MATRIX_H
#define MATRIX_H

#include<string>
using namespace std;

class BemMatrix
{
public:
    BemMatrix(int m_row, int m_column);
    BemMatrix(const BemMatrix  &m);
    BemMatrix(double ** m_array,int m_row, int m_column);
    ~BemMatrix();

    BemMatrix & operator=(const  BemMatrix & m);
    BemMatrix operator+(const BemMatrix & m);
    BemMatrix operator-(const BemMatrix & m);
    BemMatrix operator*(const BemMatrix & m);

    int row() { return m_row; }
    int colum() { return m_column; }
    double operator()(int m_row, int m_column);
    void set(int m_row, int m_column, double value);

    string toString();
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
    double **m_array;
};

class Vector : public BemMatrix
{
public:
    Vector(int n);
    void set(int n, double value) {  BemMatrix::set(0, n, value); }
    double operator()(int n) { return BemMatrix::operator() (0, n); }
};


#endif // MATRIX_H
