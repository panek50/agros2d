#include <QtCore/QCoreApplication>
#include <QTextStream>
#include "bem_matrix.h"

using namespace std;

int main()
{

    QTextStream qout(stdout);
    int row = 3;
    int column = 3;

    Matrix b(row,column);
    int k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            b.set(i, j, k);
            k++;
        }
    }


    Matrix c = b + b;
    qout << b.toString();
    qout << c.toString();
//    Vector x(4);
//    x.set(0, 1);

//      b(1, 1) = 10;
//      qout << b(1, 1);

//      qout << "\n";
//      qout << (x == b);
//      qout << "\n";
//      qout << "-----\n";
//      qout << b.toString();
//      qout << "-----\n";
//      qout << (b * x).toString();
      qout.flush();
    return 0;
}
