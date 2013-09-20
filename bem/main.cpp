#include <QtCore/QCoreApplication>
#include <QTextStream>
#include "bem_matrix.h"

using namespace std;

int main()
{

    QTextStream qout(stdout);
//    int row = 3;
//    int column = 3;

//    BemMatrix b(row,column);
//    int k = 0;
//    for (int i = 0; i < row; i++)
//    {
//        for (int j = 0; j < column; j++)
//        {
//            b(i, j) = k;
//            k++;
//        }
//    }

//    BemMatrix a(3, 1);




//    a(0, 0) = 1;
//    a(1, 0) = 2;
//    a(2, 0) = 3;

//    BemVector l(a);

//    c(0) = 1;
//    c(1) = 2;
//    c(2) = 3;
//    try{
//        BemVector d = (a * c);
//    }
//    catch (BemException e)
//    {
//        qout << e.what() << " Exception on line: " << __LINE__ << "." << " File: " << __FILE__ "\n";
//    }


//    Node y,z;
//    z(0) = 2;
//    z(1) = 3;
//    y(0) = 2;
//    y(1) = 3;
    BemMatrix a(2,2);
    a(0,0) = 1;
    a(0,1) = 0;
    a(1,0) = 2;
    a(1,1) = 3;

    BemVector b(2);
    b(0) = 3;
    b(1) = 2;

    BemVector d(3);
    d = a * b;
    d += b;
    d = a.solve(b);

    qout << d.toString();


    // qout << (y.rotate(3.141592654).toString());
    return 0;
}
