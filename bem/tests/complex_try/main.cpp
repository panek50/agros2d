#include <QCoreApplication>
#include <QDebug>
#include <boost/math/special_functions/hankel.hpp>
#include "../../bem_matrix.h"
#include <cblas.h>

int main()
{    
    // Experimnet with complex numbers

    // std::complex<double> x(1, 0);
    // x.imag() = 10 ;
    // x.real() = 10;
    // x = cos(x);
    // qDebug() << x.real();
    // qDebug() << x.imag();

    // Cyclic Hankel function
    // std::complex<double> z;
    // z = boost::math::cyl_hankel_1(10, 10);
    // qDebug() << z.imag();

    BemComplexMatrix a(2, 2);
    BemComplexVector b(2);
    BemComplexVector c(2);
    a(0, 0) = 5 + 1J;

//    a(0, 0).real() = 1;
//    a(0, 0).imag() = 1;

    a(1, 1).real() = 1;
    a(1, 1).imag() = 1;

    a(0, 1).real() = 0;
    a(0, 1).imag() = 1;

    a(1, 0).real() = 0;
    a(1, 0).imag() = 1;


    b(0).real() = 1;
    b(1).real() = 1;

    c(0).real() = 2;
    c(0).imag() = 1;


    std::complex<double> d;
    cblas_zdotu_sub(2, (std::complex<double> *) b.m_array, 2, (std::complex<double> *) c.m_array, 2, (double *) &d);

    qDebug() << a.toString();
    qDebug() << b.toString();

    qDebug() << a.solve(b).toString();

    qDebug() << d.imag();
    qDebug() << d.real();
    d = exp(d);
    qDebug() << d.imag();
    qDebug() << d.real();
}
