#include <QString>
#include <QDebug>

#include <assert.h>
#include <math.h>

#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <lapacke.h>
#include <cblas.h>

#include <complex>

#include "bem_matrix.h"

template<class Type>
BemMatrix<Type>::BemMatrix(int row, int column)
{    
    this->createArray(row,column);
    this->clear();
}

template<class Type>
BemMatrix<Type>::BemMatrix(Type array[], int row, int column)
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

template<class Type>
BemMatrix<Type>::BemMatrix(const BemMatrix & m)
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

template<class Type>
BemMatrix<Type>::~BemMatrix()
{
    delete[] m_array;
}


template <class Type>
BemMatrix<Type> & BemMatrix<Type>::operator=(const BemMatrix<Type> & m)
{
    if(this == &m)
        return *this;
    
    delete[] m_array;
    
    this->m_column = m.m_column;
    this->m_row = m.m_row;
    
    m_array = new Type[m_row * m_column];

    for(int i = 0; i < m_row; i++)
    {
        for(int j = 0; j < m_column; j++)
        {
            m_array[j + i * m_column] = m.m_array[j + i * m_column];
        }
    }
    return *this;
}

template <class Type>
bool BemMatrix<Type>::operator==(const BemMatrix & m) const
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

template <class Type>
BemMatrix<Type>  BemMatrix<Type>::operator+(const BemMatrix<Type> & m) const
{
    int i,j;
    BemMatrix<Type> mtemp(m.m_row, m.m_column);
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j < m_column; j++)
        {
            mtemp.m_array[j + i * m_column] = m.m_array[j + i * m_column] + m_array[j + i * m_column];
        }
    }
    return mtemp;
}

template <class Type>
BemMatrix<Type> BemMatrix<Type>::operator-(const BemMatrix<Type> & m)
{
    int i,j;
    BemMatrix<Type>  mtemp(m.m_row,m.m_column);
    
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j <  m_column; j++)
        {
            mtemp.m_array[j + i * m_column]=this->m_array[j + i * m_column]-m.m_array[j + i * m_column];
        }
    }
    return mtemp;
}


template <class Type>
BemMatrix<Type> BemMatrix<Type>::operator *(BemMatrix<Type> & m)
{
    int i,j,k;
    if (m_column != m.row())
        throw BemException();
    {
        BemMatrix<Type> mtemp(m_row, m.column());
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

template <class Type>
BemVector<Type> BemMatrix<Type>::operator *(BemVector<Type> & v)
{        
    BemVector<Type> vtemp(v.n());
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

template <class Type>
BemMatrix<Type> BemMatrix<Type>::operator *(const Type & x) const
{
    BemMatrix<Type> mtemp(m_row, m_column);
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

template <class Type>
BemMatrix<Type> operator *(Type x, BemMatrix<Type> m)
{
    BemMatrix<Type> mtemp(m.row(), m.column());
    for(int i = 0; i  <  m.row(); i++)
    {
        for(int j = 0; j  <  m.column(); j++)
        {
            mtemp(i, j) = m(i,j) * x;
        }
    }
    return mtemp;
}

template <class Type>
Type & BemMatrix<Type>::operator() (unsigned row, unsigned column) const
{
    return m_array[m_column * row + column];
}

template <class Type>
void BemMatrix<Type>::set(int row, int column, Type value)
{
    m_array[m_column * row + column] = value;
}

template <>
BemVector<double> BemMatrix<double>::solve(BemVector<double> &v)
{
    // Todo: copy of original matrix is performed, decide if this is usefull or not.
    BemVector<double> result(v);
    BemMatrix<double> temp(*this);
    int * pivots = new int[m_row];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, m_row, 1, temp.m_array, m_row, pivots, result.m_array, 1);
    delete pivots;
    return result;
}


template <>
BemVector<std::complex<double> > BemMatrix<std::complex<double> >::solve(BemVector<std::complex<double> > &v)
{
    // Todo: copy of original matrix is performed, decide if this is usefull or not.
    BemVector<std::complex<double> > result(v);
    BemMatrix<std::complex<double> > temp(*this);
    int * pivots = new int[m_row];
    LAPACKE_zgesv(LAPACK_ROW_MAJOR, m_row, 1, temp.m_array, m_row, pivots, result.m_array, 1);
    delete pivots;
    return result;
}

template <>
QString BemMatrix<double>::toString()
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

template <>
QString BemMatrix<std::complex<double> >::toString()
{
    int i,j;
    QString out = "";

    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j  <  m_column; j++)
        {
            out.append(QString("%1 + %2j ").arg(m_array[j + i * m_column].real())
                                            .arg(m_array[j + i * m_column].imag()));
        }
        out += "\n";
    }

    return out;
}


template <class Type>
void BemMatrix<Type>::createArray(int row, int column)
{
    this->m_column = column;
    this->m_row = row;
    m_array = new Type [row * column];
}

template <class Type>
void BemMatrix<Type>::clear()
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

template <class Type>
void BemMatrix<Type>::eye()
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

template <class Type>
BemVector<Type>::BemVector(int n)
{
    this->m_n = n;
    m_array = new Type[n];
    this->clear();
}

template <class Type>
BemVector<Type>::BemVector(const BemVector & v)
{
    m_n = v.n();
    m_array = new Type[m_n];
    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = v.m_array[i];
    }
}

//template <class Type>
//double BemVector<Type>::length() const
//{
//    double length = 0;
//    for(int i =  0; i < m_n; i++)
//    {
//        length += m_array[i] * m_array[i];
//    }
//    return sqrt(length);
//}

template <class Type>
void BemVector<Type>::operator+=(BemVector & v)
{
    if(v.n() != m_n)
        throw;
    for(int i = 0; i < m_n; i++)
    {
        m_array[i] = m_array[i] + v(i);
    }
}

template <class Type>
BemVector<Type> BemVector<Type>::operator* (const Type & x) const
{
    BemVector v(m_n);
    for(int i = 0; i < m_n; i++)
    {
        v(i) = m_array[i] * x;
    }
    return v;
}

template <class Type>
BemVector<Type> BemVector<Type>::operator* (BemMatrix<Type> & m) const
{
    {
        BemVector<Type> vtemp(m_n);
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

template <class Type>
BemVector<Type> BemVector<Type>::operator+(const BemVector & m) const
{
    BemVector<Type> vtemp(m_n);

    for(int i = 0; i  <  m_n; i++)
    {
        vtemp.m_array[i] = m.m_array[i] + m_array[i];

    }
    return vtemp;
}

template <class Type>
Type & BemVector<Type>::operator () (int i)
{
    return m_array[i];
}

template <class Type>
BemVector<Type> & BemVector<Type>::operator=(const BemVector<Type> & v)
{
    if(this == &v)
        return *this;

    delete[] m_array;

    m_n = v.n();
    m_array = new Type[m_n];

    for(int i = 0; i < m_n; i++)
    {
            m_array[i] = v.m_array[i];
    }
    return *this;
}

template <class Type>
void BemVector<Type>::clear()
{    
    for(int i = 0; i  <  m_n; i++)
    {
        m_array[i] = 0;
    }
}

template <>
QString BemVector<double>::toString()
{
    QString out = "";
    for(int i = 0; i < m_n; i++)
    {
        out += QString::number(m_array[i]);
        out += " ";
    }
    return out;
}

template <>
QString BemVector<std::complex<double> >::toString()
{
    QString out = "";
    for(int i = 0; i < m_n; i++)
    {
        out += QString("%1 + %2j \n").arg(m_array[i].real())
                .arg(m_array[i].imag());
    }
    return out;
}

template class BemMatrix<double>;
template class BemMatrix<std::complex<double> >;
template class BemVector<double>;
template class BemVector<std::complex<double> >;
