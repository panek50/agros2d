#include <QString>
#include <QDebug>

#include <assert.h>
#include <math.h>
#include <complex>

#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <lapacke.h>
#include <cblas.h>

#include "bem_matrix.h"

BemMatrix::BemMatrix(int row, int column)
{    
    this->createArray(row,column);
    this->clear();
}

BemMatrix::BemMatrix(double array[], int row, int column)
{
    int i,j;
    this->createArray(row, column);
    for(i = 0; i  <  this->m_row; i++)
    {
        for(j = 0; j < this->m_column; j++)
        {
            this->m_array[j + i * m_column] = array[j + i * m_column];
        }
    }
}

BemMatrix::BemMatrix(const BemMatrix & m)
{    
    this->createArray(m.m_row,m.m_column);
    for(int i = 0; i  <  this->m_row; i++)
    {
        for(int j = 0; j < this->m_column; j++)
        {
            this->m_array[j + i * m_column] = m.m_array[j + i * m_column];
        }
    }
}


BemMatrix::~BemMatrix()
{
    delete[] m_array;
}


BemMatrix & BemMatrix::operator=(const BemMatrix & m)
{
    if(this == &m)
        return *this;
    
    delete[] m_array;
    
    this->m_column = m.m_column;
    this->m_row = m.m_row;
    
    m_array = new double[m_row * m_column];

    for(int i = 0; i < m_row; i++)
    {
        for(int j = 0; j < m_column; j++)
        {
            m_array[j + i * m_column] = m.m_array[j + i * m_column];
        }
    }
    return *this;
}

bool BemMatrix::operator==(const BemMatrix & m) const
{
    for(int i = 0; i < m_row; i++)
    {
        for(int j = 0; j < m_column; j++)
        {
            if( m_array[j + i * m_column] != m.m_array[j + i * m_column])
                return false;
        }
    }
    return true;
}

BemMatrix  BemMatrix::operator+(const BemMatrix & m) const
{
    int i,j;
    BemMatrix mtemp(m.m_row, m.m_column);
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j < m_column; j++)
        {
            mtemp.m_array[j + i * m_column] = m.m_array[j + i * m_column] + m_array[j + i * m_column];
        }
    }
    return mtemp;
}


BemMatrix BemMatrix::operator-(const BemMatrix & m)
{
    int i,j;
    BemMatrix  mtemp(m.m_row,m.m_column);
    
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j <  m_column; j++)
        {
            mtemp.m_array[j + i * m_column]=this->m_array[j + i * m_column]-m.m_array[j + i * m_column];
        }
    }
    return mtemp;
}



BemMatrix BemMatrix::operator *(BemMatrix & m)
{
    int i,j,k;
    if (m_column != m.row())
        throw BemException();
    {
        BemMatrix mtemp(m_row, m.column());
        mtemp.clear();

        for(i = 0; i  <  m_row; i++)
        {
            for(j = 0; j  <  m.column(); j++)
            {
                for(k = 0; k  <  m_column; k++)
                {
                    mtemp(i, j) += (* this)(i, k) * m(k, j);
                }
            }
        }
        return mtemp;
    }
}

BemVector BemMatrix::operator *(BemVector & v)
{        
    BemVector vtemp(v.n());
    vtemp.clear();

    for(int i = 0; i  <  m_row; i++)
    {
        for(int k = 0; k  <  m_column; k++)
        {
            vtemp(i) += m_array[i  * m_column + k] * v(k);
        }
    }
    return vtemp;
}

BemMatrix BemMatrix::operator *(const double & x) const
{
    BemMatrix mtemp(m_row, m_column);
    mtemp.clear();
    for(int i = 0; i  <  m_row; i++)
    {
        for(int j = 0; j  <  m_column; j++)
        {
            mtemp.m_array[j + i * m_column] = m_array[j + i * m_column] * x;
        }
    }
    return mtemp;
}

BemMatrix operator *(double x, BemMatrix m)
{
    BemMatrix mtemp(m.row(), m.column());
    for(int i = 0; i  <  m.row(); i++)
    {
        for(int j = 0; j  <  m.column(); j++)
        {
            mtemp(i, j) = m(i,j) * x;
        }
    }
    return mtemp;
}

double & BemMatrix::operator() (unsigned row, unsigned column) const
{
    return m_array[m_column * row + column];
}

void BemMatrix::set(int row, int column, double value)
{
    m_array[m_column * row + column] = value;
}

BemVector BemMatrix::solve(BemVector &v)
{
    // Todo: copy of original matrix is performed, decide if this is usefull or not.
    BemVector result(v);
    BemMatrix temp(*this);
    int * pivots = new int[m_row];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, m_row, 1, temp.m_array, m_row, pivots, result.m_array, 1);
    delete pivots;
    return result;
}

QString BemMatrix::toString()
{
    int i,j;
    QString out = "";

    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j  <  m_column; j++)
        {
            out.append(QString("%1 ").arg(m_array[j + i * m_column]));
        }
        out += "\n";
    }

    return out;
}


void BemMatrix::createArray(int row, int column)
{
    this->m_column = column;
    this->m_row = row;
    m_array = new double [row * column];
}

void BemMatrix::clear()
{
    int i,j;
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j  <  m_column; j++)
        {
            m_array[j + i * m_column] = 0;
        }
    }
}

void BemMatrix::eye()
{
    int i,j;
    for(i = 0; i < this->m_row; i++)
    {
        for(j=0; j < this->m_column; j++)
        {
            this->m_array[j + i * m_column]=0;
            if(i==j)
                this->m_array[j + i * m_column]=1;
        }
    }
}


BemVector::BemVector(int n)
{
    this->m_n = n;
    m_array = new double[n];
    this->clear();
}

BemVector::BemVector(const BemVector & v)
{
    m_n = v.n();
    m_array = new double[m_n];
    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = v.m_array[i];
    }
}

double BemVector::length() const
{
    double length = 0;
    for(int i =  0; i < m_n; i++)
    {
        length += m_array[i] * m_array[i];
    }
    return sqrt(length);
}

void BemVector::operator+=(BemVector & v)
{
    if(v.n() != m_n)
        throw;
    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = m_array[i] + v(i);
    }
}

BemVector BemVector::operator* (const double & x) const
{
    BemVector v(m_n);
    for(int i = 0; i < m_n; i++)
    {
        v(i) = m_array[i] * x;
    }
    return v;
}


BemVector BemVector::operator* (BemMatrix & m) const
{
    {
        BemVector vtemp(m_n);
        vtemp.clear();

        for(int i = 0; i  <  m_n; i++)
        {
            for(int j = 0; j  <  m.row(); j++)
            {
                vtemp(i) += m_array[j] * m(j, i);
            }
        }
        return vtemp;
    }
}

BemVector BemVector::operator+(const BemVector & m) const
{
    BemVector vtemp(m_n);

    for(int i = 0; i  <  m_n; i++)
    {
        vtemp.m_array[i] = m.m_array[i] + m_array[i];

    }
    return vtemp;
}

BemComplexVector & BemComplexVector::operator=(const BemComplexVector & v)
{
    if(this == &v)
        return *this;

    delete[] m_array;

    m_n = v.m_n;
    m_array = new std::complex<double> [m_n];

    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = v.m_array[i];
    }
    return *this;
}


double & BemVector::operator () (int i)
{
    return m_array[i];
}


void BemVector::clear()
{    
    for(int i = 0; i  <  m_n; i++)
    {
        m_array[i] = 0;
    }
}

QString BemVector::toString()
{
    QString out = "";
    for(int i = 0; i < m_n; i++)
    {
        out += QString::number(m_array[i]);
        out += " ";
    }
    return out;
}


BemComplexMatrix::BemComplexMatrix(int row, int column)
{
    this->createArray(row, column);
    this->clear();
}

BemComplexMatrix::BemComplexMatrix(const BemComplexMatrix & m)
{
    this->createArray(m.m_row,m.m_column);
    for(int i = 0; i  <  this->m_row; i++)
    {
        for(int j = 0; j < this->m_column; j++)
        {
            this->m_array[j + i * m_column] = m.m_array[j + i * m_column];
        }
    }
}


void BemComplexMatrix::createArray(int row, int column)
{
    m_column = column;
    m_row = row;
    m_array = new std::complex<double>[row * column];
}

void BemComplexMatrix::clear()
{

    for(int i = 0; i  <  m_row; i++)
    {
        for(int j = 0; j  < m_column; j++)
        {
          m_array[(j + m_column * i)].real() = 0;
          m_array[(j + m_column * i)].imag() = 0;
        }
    }
}

QString BemComplexMatrix::toString()
{

    QString out = "";
    for(int i = 0; i  <  m_row; i++)
    {
        for(int j = 0; j  < m_column; j++)
        {            
            out.append(QString("%1 + ").arg((m_array[(j + i * m_column)].real())));
            out.append(QString("%1j ").arg((m_array[(j + i * m_column)].imag())));
        }
        out += "\n";
    }

    return out;
}

std::complex<double> & BemComplexMatrix::operator() (unsigned row, unsigned column)
{
    return (m_array[(m_column * row + column)]);
}

void BemComplexMatrix::set(int row, int column, std::complex<double> number)
{    
    m_array[column + row * m_column] = number;
}

BemComplexMatrix::~BemComplexMatrix()
{
    delete[] m_array;
}


BemComplexMatrix & BemComplexMatrix::operator=(const BemComplexMatrix & m)
{
    if(this == &m)
        return *this;

    delete[] m_array;

    this->m_column = m.m_column;
    this->m_row = m.m_row;

    m_array = new std::complex<double>[m_row * m_column];

    for(int i = 0; i < m_row; i++)
    {
        for(int j = 0; j < m_column; j++)
        {
            m_array[j + i * m_column] = m.m_array[j + i * m_column];
        }
    }
    return *this;
}

BemComplexVector BemComplexMatrix::solve(BemComplexVector &v)
{
    // Todo: copy of original matrix is performed, decide if this is usefull or not.
    BemComplexVector result(v);
    BemComplexMatrix temp(*this);
    int * pivots = new int[m_row];
    LAPACKE_zgesv(LAPACK_ROW_MAJOR, m_row, 1, temp.m_array, m_row, pivots, result.m_array, 1);
    delete pivots;
    return result;
}


BemComplexVector::BemComplexVector(int n)
{
    this->m_n = n;
    m_array = new std::complex<double>[n];
    this->clear();
}

BemComplexVector::BemComplexVector(const BemComplexVector & v)
{
    m_n = v.n();
    m_array = new std::complex<double>[m_n];
    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = v.m_array[i];
    }
}

//double BemVector::length() const
//{
//    double length = 0;
//    for(int i =  0; i < m_n; i++)
//    {
//        length += m_array[i] * m_array[i];
//    }
//    return sqrt(length);
//}

//void BemComplexVector::operator+=(BemComplexVector & v)
//{
//    if(v.n() != m_n)
//        throw;

//    for(int i = 0; i < m_n; i++)
//    {
//        m_array[i] = m_array[i] + v(i).real();
//        m_array[i] = m_array[i + 1] + v(i).imag();
//    }
//}

//BemComplexVector BemComplexVector::operator* (const double & x) const
//{
//    BemComplexVector v(m_n);
//    for(int i = 0; i < m_n; i++)
//    {
//        v(i) = m_array[i] * x;
//    }
//    return v;
//}


//BemVector BemVector::operator* (BemMatrix & m) const
//{
//    {
//        BemVector vtemp(m_n);
//        vtemp.clear();

//        for(int i = 0; i  <  m_n; i++)
//        {
//            for(int j = 0; j  <  m.row(); j++)
//            {
//                vtemp(i) += m_array[j] * m(j, i);
//            }
//        }
//        return vtemp;
//    }
//}

//BemVector BemVector::operator+(const BemVector & m) const
//{
//    BemVector vtemp(m_n);

//    for(int i = 0; i  <  m_n; i++)
//    {
//        vtemp.m_array[i] = m.m_array[i] + m_array[i];

//    }
//    return vtemp;
//}

BemVector & BemVector::operator=(const BemVector & v)
{
    if(this == &v)
        return *this;

    delete[] m_array;

    m_n = v.m_n;
    m_array = new double[m_n];

    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = v.m_array[i];
    }
    return *this;
}


std::complex<double> & BemComplexVector::operator () (int i)
{
    return m_array[i];
}


void BemComplexVector::clear()
{
    for(int i = 0; i  <  m_n; i++)
    {
        m_array[i] = 0;
    }
}

QString BemComplexVector::toString()
{
    QString out = "";
    for(int i = 0; i < m_n; i++)
    {
        out += QString::number(m_array[i].real());
        out += " + " + QString::number(m_array[i].imag()) + "j\n";
    }
    return out;
}
