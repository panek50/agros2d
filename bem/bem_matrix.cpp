#include "bem_matrix.h"
#include <sstream>
#include <iostream>
#include <assert.h>



using namespace std;

BemMatrix::BemMatrix(int row, int column)
{    
    this->createArray(row,column);
}


BemMatrix::BemMatrix(const BemMatrix & m)
{
    int i,j;
    this->createArray(m.m_row,m.m_column);
    for(i = 0; i < this->m_row; i++)
    {
        for(j=0;j<this->m_column;j++)
        {
            this->m_array[i][j]=m.m_array[i][j];
        }
    }
}


BemMatrix::BemMatrix(double **array,int row, int column)
{
    int i,j;
    this->createArray(row, column);
    for(i=0;i<this->m_row;i++)
    {
        for(j=0;j<this->m_column;j++)
        {
            this->m_array[i][j]=array[i][j];
        }
    }
}

BemMatrix::~BemMatrix()
{
    int i;
    for(i=0;i<m_row;i++)
    {
        delete[] m_array[i];
    }
    delete[] m_array;
}


BemMatrix & BemMatrix::operator=(const BemMatrix & m)
{
    int i,j;
    if(this == &m)
        return *this;
    
    for(i=0;i<m_row;i++)
    {
        delete m_array[i];
    }
    delete m_array;
    
    this->m_column = m.m_column;
    this->m_row = m.m_row;
    
    m_array = new double* [m_row];
    for(i=0;i<m_row;i++)
    {
        m_array[i] = new double[m_column];
    }
    for(i=0;i<m_row;i++)
    {
        for(j=0;j<m_row;j++)
        {
            m_array[i][j]=m.m_array[i][j];
        }
    }
    return *this;
}

BemMatrix BemMatrix::operator+(const BemMatrix & m)
{
    int i,j;
    BemMatrix  mtemp(m.m_row,m.m_column);
    
    for(i=0;i<m_row;i++)
    {
        for(j=0;j<m_row;j++)
        {
            mtemp.m_array[i][j]=m.m_array[i][j]+this->m_array[i][j];
        }
    }
    return mtemp;
}


BemMatrix BemMatrix::operator-(const BemMatrix & m)
{
    int i,j;
    BemMatrix  mtemp(m.m_row,m.m_column);
    
    for(i = 0; i < m_row; i++)
    {
        for(j = 0; j< m_column; j++)
        {
            mtemp.m_array[i][j]=this->m_array[i][j]-m.m_array[i][j];
        }
    }
    return mtemp;
}


BemMatrix BemMatrix::operator *(const BemMatrix & m)
{
    int i,j,k;

    if (m.m_row != m_column)
    {
        if (m.m_row == 1)
        {
            BemMatrix mtemp(1, m.m_column);
            mtemp.clear();

            for(i = 0; i < m_row; i++)
            {
                for(j = 0; j < m_column; j++)
                {
                    mtemp.m_array[0][i] = mtemp.m_array[0][i] + this->m_array[i][j] * m.m_array[0][j];
                }
            }
            return mtemp;
        }
        else throw;
    }
    else
    {
        BemMatrix mtemp(m_row, m.m_column);
        mtemp.clear();

        for(i = 0; i < m_row; i++)
        {
            for(j = 0; j < m_column; j++)
            {
                for(k = 0; k < m_column; k++)
                {
                    mtemp.m_array[i][j]=mtemp.m_array[i][j] + this->m_array[i][k] * m.m_array[k][j];
                }
            }
        }
        return mtemp;
    }
    assert(0);
}

double BemMatrix::operator() (int row, int column)
{
    return m_array[row][column];
}

void BemMatrix::set(int row, int column, double value)
{
    m_array[row][column] = value;
}

string BemMatrix::toString()
{
    int i,j;
    stringstream out;
    for(i=0;i<m_row;i++)
    {
        for(j=0;j<m_column;j++)
        {
            out << m_array[i][j] << " ";
        }
        out << endl;
    }
    return out.str();
}


int BemMatrix::countColumn()
{
    return this->m_column;
}

int BemMatrix::countRow()
{
    return this->m_row;
}

void BemMatrix::createArray(int row, int column)
{
    int i;
    this->m_column = column;
    this->m_row = row;
    m_array = new double* [row];
    for(i=0;i<row;i++)
    {
        m_array[i] = new double[column];
    }
}

void BemMatrix::clear()
{
    int i,j;
    for(i = 0; i < m_row; i++)
    {
        for(j = 0; j < m_column; j++)
        {
            m_array[i][j] = 0;
        }
    }
}

void BemMatrix::eye()
{
    int i,j;
    for(i=0;i<this->m_row;i++)
    {
        for(j=0;j<this->m_column;j++)
        {
            this->m_array[i][j]=0;
            if(i==j)
                this->m_array[i][j]=1;
        }
    }
}


void BemMatrix::gaussElim()
{
    int i,j,k;
    int maxPosition;
    double max;
    double * ptemp;
    
    
    for(k=0;k<this->m_row;k++)
    {
        max=this->m_array[k][k];
        maxPosition=k;
        for(i=k;i<this->m_row;i++)
        {
            if((max<m_array[i][k])||(-max>m_array[i][k]))
            {
                max = m_array[i][k];
                maxPosition=i;
            }
            
        }
        ptemp=m_array[maxPosition];
        m_array[maxPosition]=m_array[k];
        m_array[k]=ptemp;
        
        for(j=k;j<this->m_column;j++)
        {
            this->m_array[k][j]= this->m_array[k][j]/max;
        }
        
        for(i=k+1;i<this->m_row;i++)
        {
            for(j=this->m_column-1;j>=k;j--)
            {
                this->m_array[i][j] =  this->m_array[i][j]-this->m_array[i][k]*this->m_array[k][j];
            }
            m_array[i][k]=0;
        }
    }
    
    
    for(k=m_row-1;k>=0;k--)
    {
        {
            for(i=k-1;i>=0;i--)
            {
                for(j=m_column-1;j>=i;j--)
                {
                    m_array[i][j] = m_array[i][j] - m_array[k][j]*m_array[i][k];
                    
                }
            }
            
        }
    }
    
}


void BemMatrix::luFactorisation()
{
    int i,j,k;
    int maxPosition;
    double max;
    double * ptemp;
    double delta;
    double sum;
    
    for(k=0;k<this->m_row-1;k++)
    {
        max=this->m_array[k][k];
        maxPosition=k;
        for(i=k+2;i<this->m_row;i++)
        {
            if((max<m_array[i][k])||(-max>m_array[i][k]))
            {
                max = m_array[i][k];
                maxPosition=i;
            }
            
        }
        ptemp=m_array[maxPosition];
        m_array[maxPosition]=m_array[k];
        m_array[k]=ptemp;
        for(i=k+1;i<this->m_row;i++)
        {
            for(j=k+1;j<m_column;j++)
            {
                
                this->m_array[i][j] =  this->m_array[i][j]-this->m_array[i][k]*this->m_array[k][j]/max;
            }
            m_array[i][k]=this->m_array[i][k]*this->m_array[k][m_column-(k+1)] /max/m_array[k][m_column-(k+1)];
        }
    }
    cout << this->toString();


    for(j=0;j<m_row;j++)
    {
        for(i=j;i<m_column;i++)
        {
            delta=(i==j);
            cout << delta << endl;
            sum=0;
            for(k=j;k<=(i-1);k++)
            {
                sum = sum + m_array[i][k]*m_array[k][j];
            }
            m_array[i][j]=(i==j)-sum;
        }
    }
}



Vector::Vector(int n) : BemMatrix(1, n) {}

