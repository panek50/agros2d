#include <QtCore/QCoreApplication>
#include <iostream>
#include "bem_matrix.h"

using namespace std;

int main()
{
    int row=4;
    int column=4;

    BemMatrix b(row,column);
    int k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
             b.set(i, j, k);
             k++;
        }
    }

    // Matrix c = b + b + b;
    Vector x(4);
    x.set(0, 1);


//    cout << b.toString();
//    b.luFactorisation();
//      cout << b.toString();
//     cout << b(1,1) << "\n";
//     cout << x.toString();
//     cout << (b*b).toString();
     cout << (x*b).toString();
     cout << "-----\n";
     cout << b.toString();
     cout << "-----\n";
     cout << (b*x).toString();

}
