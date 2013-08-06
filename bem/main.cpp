#include <QtCore/QCoreApplication>
#include <QTextStream>
#include "bem_matrix.h"

using namespace std;

int main()
{

    QTextStream qout(stdout);
    int row = 3;
    int column = 3;

    BemMatrix b(row,column);
    int k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            b(i, j) = k;
            k++;
        }
    }


//    BemMatrix c = 3 * b * 3;
    BemVector a(3);
    a(0) = 0;
    a(1) = 1;
    a(2) = 2;
    BemVector d = a * a;
//    qout << "a = " << a.toString();
//    qout << "c = " << c.toString();
    qout << "result = " << d.toString();
//    Node x;
//    x(0) = 1;
//    Node y;
//    y(0) = 1;
//    x = y * 5;
//    qout << "Node:" << x.toString();
    qout.flush();
    return 0;
}
