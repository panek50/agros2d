#include <math.h>

#include  <QString>
#include  <QDebug>
#include  <assert.h>

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

BemMatrix BemMatrix::operator+(const BemMatrix & m) const
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



BemMatrix BemMatrix::operator *(const BemMatrix & m)
{
    int i,j,k;
    assert(m_column == m.row());
    {
        BemMatrix mtemp(m_column, m.row());
        mtemp.clear();

        for(i = 0; i  <  m_row; i++)
        {
            for(j = 0; j  <  m_column; j++)
            {
                for(k = 0; k  <  m_column; k++)
                {
                    mtemp.m_array[j + i * m_column] = mtemp.m_array[j + i * m_column] + this->m_array[i * m_column + k] * m.m_array[k * m_column + j];
                }
            }
        }
        return mtemp;
    }
}

BemMatrix BemMatrix::operator *(const double & x) const
{
    BemMatrix mtemp(m_row, m_column);
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


BemVector::BemVector(int n) : BemMatrix(1, n) {}



BemVector::BemVector(BemMatrix m) : BemMatrix(1, m.column())
{
    if (!((m.row() == 1) || (m.column() == 1)))
    {
        assert(0);
    }

    for(int i = 0; i < this->column(); i++)
    {
        this->operator()(i) = m(0, i);
    }
}

BemVector::BemVector(const BemVector & v) : BemMatrix(1, v.column())
{
    for(int i = 0; i < this->column(); i++)
    {
        m_array[i] = v.m_array[i];
    }
}

BemVector operator *(double x, BemVector v)
{
    BemVector vtemp(v.column());
    for(int i = 0; i  <  v.column(); i++)
    {
        vtemp(i) = v(i) * x;
    }

    return vtemp;
}



BemVector BemVector::operator *(BemMatrix  m)
{
    int i, j;
    BemVector v = * this;

    assert(m.column() == this->column());
    {
        BemVector vtemp(this->column());
        vtemp.clear();

        for(i = 0; i  <  this->column(); i++)
        {
            for(j = 0; j  <  m.column(); j++)
            {
                vtemp(i) += m(j, i) * v(i);
            }
        }
        return vtemp;
    }
}


BemVector operator *(BemMatrix m, BemVector v)
{
    int i, j;
    assert(m.column() == v.column());
    {
        BemVector vtemp(v.column());
        vtemp.clear();

        for(i = 0; i  <  v.column(); i++)
        {
            for(j = 0; j  <  m.column(); j++)
            {
                vtemp(i) += m(i, j) * v(j);
            }
        }
        return vtemp;
    }
}

double BemVector::operator *(BemVector  v)
{
    assert(v.column() == this->column());

    int i;
    double result = 0;
    BemVector & a = (*this);

    for(i = 0; i  <  this->column(); i++)
    {
        result += v(i) * a(i);
    }

    return result;
}

BemVector & BemVector::operator =(const BemVector & v)
{
    if(this == &v)
        return *this;
    BemMatrix::operator =(v);
    return *this;
}

double BemVector::length()
{
    double length = 0;
    for(int i = 0; i < this->column(); i++)
    {
        length += this->operator ()(i);
    }
    length = sqrt(length);
    return length;
}

Node & Node::operator =(const Node & n)
{
    if(this == &n)
        return *this;
    BemVector::operator =(n);
    return *this;
}

Node::Node(double x, double y) : BemVector(2)
{
    BemVector & v = (*this);
    v(0) = x;
    v(1) = y;
}

void Node::shift(double x, double y)
{
    BemVector & v = (*this);
    v(0) = v(0) + x;
    v(1) = v(1) + y;

    qDebug() << v.toString();
}

void Node::rotate(double angle)
{
    BemVector & v = (*this);
    BemMatrix m(2,2);
    m(0, 0) = cos(angle);
    m(0, 1) = - sin(angle);
    m(1, 0) = sin(angle);
    m(1, 1) = cos(angle);
    v = m * v;
    qDebug() << v.toString();
}

Node::Node(BemVector v) : BemVector(v.column())
{
    for(int i = 0; i < this->column(); i++)
    {
        this->operator()(i) = v(i);
    }
}

