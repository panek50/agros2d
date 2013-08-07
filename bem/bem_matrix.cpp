#include <math.h>

#include  <QString>
#include  <QDebug>
#include  <assert.h>
#include <lapacke.h>

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

BemMatrix::BemMatrix(const BemMatrix & m) : QObject()
{
    int i,j;
    this->createArray(m.m_row,m.m_column);
    for(i = 0; i  <  this->m_row; i++)
    {
        for(j = 0; j < this->m_column; j++)
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

double & BemMatrix::operator() (unsigned row, unsigned column)
{
    return m_array[m_column * row + column];
}

void BemMatrix::set(int row, int column, double value)
{
    m_array[m_column * row + column] = value;
}

BemMatrix BemMatrix::solve(const BemMatrix & m)
{
    // Todo: copy of original matrix is performed, decide if this is usefull or not.
    BemMatrix result(m);
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

BemVector & BemVector::transpose()
{
    int temp = row();
    setRow(column());
    setColumn(temp);
    return * this;
}

double BemVector::length() const
{
   double length = 0;
   int n;
   if (m_type == VectorTypeColumn)
       n = column();
   else
       n= row();

   for(int i =  0; i < n; i++)
   {
       length += m_array[i] * m_array[i];
   }
   return sqrt(length);
}


Node & Node::operator = (BemVector & v)
{
    if(column() != 1)
        throw BemException();

    for(int i = 0; i < row(); i++)
    {
        (*this)(i) = v(i);
    }
    return (*this);
}

Node & Node::rotate(double phi)
{
    BemMatrix m_transform(2,2);
    m_transform(0, 0) = cos(phi);
    m_transform(0, 1) = -sin(phi);
    m_transform(1, 0) = sin(phi);
    m_transform(1, 1) = cos(phi);
    (*this) = m_transform * (*this);
    return (*this);
}

double Node::distanceOf(const Node & node)
{
    BemVector dist = (*this) - node;
    return dist.length();
}

