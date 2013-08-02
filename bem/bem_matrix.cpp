#include "bem_matrix.h"
#include  <QString>
#include  <assert.h>


Matrix::Matrix(int row, int column)
{    
    this->createArray(row,column);
    this->clear();
}

Matrix::Matrix(double array[], int row, int column)
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

Matrix::Matrix(const Matrix & m)
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


Matrix::~Matrix()
{
    delete[] m_array;
}


Matrix & Matrix::operator=(const Matrix & m)
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

bool Matrix::operator==(const Matrix & m) const
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

Matrix Matrix::operator+(const Matrix & m)
{
    int i,j;
    Matrix  mtemp(m.m_row, m.m_column);
    
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j < m_column; j++)
        {
            mtemp.m_array[j + i * m_column] = m.m_array[j + i * m_column] + m_array[j + i * m_column];
        }
    }
    return mtemp;
}


Matrix Matrix::operator-(const Matrix & m)
{
    int i,j;
    Matrix  mtemp(m.m_row,m.m_column);
    
    for(i = 0; i  <  m_row; i++)
    {
        for(j = 0; j <  m_column; j++)
        {
            mtemp.m_array[j + i * m_column]=this->m_array[j + i * m_column]-m.m_array[j + i * m_column];
        }
    }
    return mtemp;
}


Matrix Matrix::operator *(const Matrix & m)
{
    int i,j,k;

    if (m.m_row != m_column)
    {
        if (m.m_row == 1)
        {
            Matrix mtemp(1, m.m_column);
            mtemp.clear();

            for(i = 0; i  <  m_row; i++)
            {
                for(j = 0; j  <  m_column; j++)
                {
                    mtemp.m_array[i] = mtemp.m_array[i] + this->m_array[j + i * m_column] * m.m_array[j];
                }
            }
            return mtemp;
        }
        else throw;
    }
    else
    {
        Matrix mtemp(m_row, m.m_column);
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
    assert(0);
}

double & Matrix::operator() (unsigned row, unsigned column)
{
    return m_array[m_column * row + column];
}

void Matrix::set(int row, int column, double value)
{
    m_array[m_column * row + column] = value;
}

QString Matrix::toString()
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


int Matrix::countColumn()
{
    return this->m_column;
}

int Matrix::countRow()
{
    return this->m_row;
}

void Matrix::createArray(int row, int column)
{
    int i;
    this->m_column = column;
    this->m_row = row;
    m_array = new double [row * column];
}

void Matrix::clear()
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

void Matrix::eye()
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


//void Matrix::gaussElim()
//{
//    int i,j,k;
//    int maxPosition;
//    double max;
//    double * ptemp;


//    for(k=0;k < this->m_row;k++)
//    {
//        max = this->m_array[k][k];
//        maxPosition=k;
//        for(i=k; i < this->m_row; i++)
//        {
//            if((max < m_array[i][k])||(-max>m_array[i][k]))
//            {
//                max = m_array[i][k];
//                maxPosition=i;
//            }

//        }
//        ptemp=m_array[maxPosition];
//        m_array[maxPosition]=m_array[k];
//        m_array[k]=ptemp;

//        for(j=k; j < this->m_column; j++)
//        {
//            this->m_array[k][j]= this->m_array[k][j]/max;
//        }

//        for(i=k+1; i < this->m_row; i++)
//        {
//            for(j=this->m_column-1; j>=k; j--)
//            {
//                this->m_array[j + i * m_column] =  this->m_array[j + i * m_column]-this->m_array[i][k]*this->m_array[k][j];
//            }
//            m_array[i][k]=0;
//        }
//    }


//    for(k=m_row-1;k>=0;k--)
//    {
//        {
//            for(i=k-1; i>=0; i--)
//            {
//                for(j=m_column-1; j>=i; j--)
//                {
//                    m_array[j + i * m_column] = m_array[j + i * m_column] - m_array[k][j]*m_array[i][k];

//                }
//            }

//        }
//    }

//}


//void Matrix::luFactorisation()
//{
//    int i,j,k;
//    int maxPosition;
//    double max;
//    double * ptemp;
//    double sum;

//    for(k=0;k < this->m_row-1;k++)
//    {
//        max=this->m_array[k][k];
//        maxPosition=k;
//        for(i=k+2; i < this->m_row; i++)
//        {
//            if((max < m_array[i][k])||(-max>m_array[i][k]))
//            {
//                max = m_array[i][k];
//                maxPosition=i;
//            }

//        }
//        ptemp=m_array[maxPosition];
//        m_array[maxPosition]=m_array[k];
//        m_array[k]=ptemp;
//        for(i=k+1; i < this->m_row; i++)
//        {
//            for(j=k+1; j < m_column; j++)
//            {

//                this->m_array[j + i * m_column] =  this->m_array[j + i * m_column]-this->m_array[i][k]*this->m_array[k][j]/max;
//            }
//            m_array[i][k]=this->m_array[i][k]*this->m_array[k][m_column-(k+1)] /max/m_array[k][m_column-(k+1)];
//        }
//    }


//    for(j = 0; j < m_row; j++)
//    {
//        for(i = j; i < m_column; i++)
//        {
//            sum = 0;
//            for(k=j;k < =(i-1);k++)
//            {
//                sum = sum + m_array[i * m_column + k] * m_array[k * m_column + j];
//            }
//            m_array[j + i * m_column]=(i==j)-sum;
//        }
//    }
//}



Vector::Vector(int n) : Matrix(1, n) {}

