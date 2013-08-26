#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <QString>
namespace Algebra
{
class Matrix
{
public:
    Matrix();
    Matrix(int rows, int columns);
    Matrix(const Matrix & matrix);
    Matrix(double** array, int rows, int columns);

    Matrix operator +(const Matrix &);
    Matrix operator -(const Matrix &);
    Matrix operator *(const Matrix &);
    Matrix operator *(double x);

    double operator() (int row, int column) { return m_array[row][column];}
    void setValue(int row, int column, double value) {m_array[row][column] = value;}

    friend Matrix operator *(double x, const Matrix & matrix);

    Matrix operator /(const Matrix &);
    void switchRow(int first, int second);
    void switchColumn(int first, int second);
    void operator =(const Matrix &);


    int rows() const {return m_rows;}
    int columns() const {return m_columns;}
    QString toString();
    ~Matrix();

private:
    double ** m_array;
    int * m_ri;
    int * m_ci;

    int m_rows;
    int m_columns;
    void allocate();
    void release();
    Matrix sum(const Matrix &, char sign);
};

}
#endif // ALGEBRA_H
