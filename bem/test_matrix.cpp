#include <QtTest/QtTest>

#include "bem_matrix.h"

class TestMatrix: public QObject
{
    Q_OBJECT
private slots:
    void testAdd();
    void testAdd_data();
};

void TestMatrix::testAdd()
{
    QFETCH(Matrix, a);
    QFETCH(Matrix, b);
    QFETCH(Matrix, result);

    QCOMPARE(a + b, result);
}

void TestMatrix::testAdd_data()
{
    double array[3 * 3] = {1, 2, 3, 4, 1, 0, 0, 0, 0};
    double res_array[3 * 3] = {2, 4, 6, 8, 2, 0, 0, 0, 0};
    Matrix ma(array, 3, 3);
    Matrix mb(array, 3, 3);
    Matrix mresult(res_array, 3, 3);

    QTest::addColumn<Matrix>("a");
    QTest::addColumn<Matrix>("b");
    QTest::addColumn<Matrix>("result");

    QTest::newRow("a + b") << ma << mb << mresult;
}


QTEST_MAIN(TestMatrix)
#include "test_matrix.moc"

