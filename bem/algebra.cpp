#include "algebra.h"
#include <QTextStream>

namespace Algebra
{
Matrix::Matrix()
{
    m_rows = 0;
    m_columns = 0;
    m_array = 0;
}

Matrix::Matrix(int rows, int columns)
{
    m_rows = rows;
    m_columns = columns;

    allocate();

    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            m_array[i][j] =  i+j;
        }
    }
}

Matrix::Matrix(const Matrix & matrix)
{
    {
        m_rows = matrix.rows();
        m_columns = matrix.columns();

        allocate();

        for(int i = 0; i < m_rows; i++)
        {
            for(int j = 0; j < m_columns; j++)
            {
                m_array[i][j] =  matrix.m_array[i][j];
            }
        }
    }
}

Matrix::Matrix(double** array, int rows, int columns)
{
    m_rows = rows;
    m_columns = columns;

    allocate();

    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            m_array[i][j] =  array[i][j];
        }
    }
}

Matrix::~Matrix()
{
    release();
}


Matrix Matrix::sum(const Matrix & matrix, char sign)
{
    Matrix output(m_rows, m_columns);
    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            output.m_array[m_ri[i]][m_ci[j]] =  sign * matrix.m_array[m_ri[i]][m_ci[j]] +  m_array[m_ri[i]][m_ci[j]];
        }
    }
    return output;
}

Matrix Matrix::operator +(const Matrix & matrix)
{
    return(sum(matrix, 1));
}

Matrix Matrix::operator -(const Matrix & matrix)
{
    return(sum(matrix, -1));
}

Matrix Matrix::operator *(const Matrix & matrix)
{
    Matrix output(m_rows, m_columns);
    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            output.m_array[m_ri[i]][m_ci[j]] = 0;
            for(int k = 0; k < m_columns; k++)
            {
                output.m_array[m_ri[i]][m_ci[j]] +=  matrix.m_array[m_ri[i]][m_ci[j]] * m_array[m_ri[i]][m_ci[j]];
            }
        }
    }
    return output;
}


Matrix Matrix::operator *(double x)
{
    Matrix output(m_rows, m_columns);
    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            output.m_array[i][j] = x * m_array[i][j];
        }
    }
    return output;
}

Matrix operator *(double x, const Matrix & matrix)
{
    int rows = matrix.m_columns;
    int columns = matrix.m_rows;
    Matrix output(rows, columns);
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            output.m_array[i][j] = x * matrix.m_array[i][j];
        }
    }
    return output;
}

void Matrix::operator =(const Matrix & matrix)
{
    release();

    m_rows = matrix.m_rows;
    m_columns = matrix.m_columns;
    allocate();

    for(int i = 0; i < m_rows; i++)
    {
        for(int j = 0; j < m_columns; j++)
        {
            m_array[i][j] =  matrix.m_array[i][j];
        }
    }
}


QString Matrix::toString()
{
    QString output = "";

    // ToDo: assert could be here
    if (m_array != 0)
    {
        for(int i = 0; i < m_rows; i++)
        {
            for(int j = 0; j < m_columns; j++)
            {
                output +=  QString::number(m_array[m_ri[i]][m_ci[j]]) + ", " ;
            }
            output += "\n";
        }
    }

    return output;
}

void Matrix::switchRow(int first, int second)
{
    // ToDo: Check array boundaries
    int index;

    if (first != second)
    {
        index = m_ri[first];
        m_ri[first] = m_ri[second];
        m_ri[second] = index;
    }
}

void Matrix::switchColumn(int first, int second)
{
    // ToDo: Check array boundaries
    int index;

    if (first != second)
    {
        index = m_ci[first];
        m_ci[first] = m_ci[second];
        m_ci[second] = index;
    }
}

// Private
void Matrix::allocate()
{
    m_array = new double*[m_rows];
    m_ri = new int[m_rows];
    m_ci = new int[m_columns];

    for(int i = 0; i < m_rows; i++)
    {
        m_array[i] = new double[m_columns];
        m_ri[i] = i;
    }

    for(int j = 0; j < m_columns; j++)
    {
        m_ci[j] = j;
    }
}

void Matrix::release()
{
    if(m_array != 0)
    {
        for(int i = 0; i < m_rows; i++)
        {
            delete[] m_array[i];
        }
        delete[] m_array;
        delete[] m_ri;
        delete[] m_ci;

        m_array = 0;
        m_rows = 0;
        m_columns = 0;
    }
}

}
