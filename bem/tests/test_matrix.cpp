#include <qt5/QtTest/QTest>
#include "../bem_matrix.h"

/// ToDo: Create test for zero order
/*
0.269473 -0.0707733 0.0532936 -0.111945 0.0532936 -0.0707733
-0.0750147 0.269473 0.00719838 0.0532936 -0.0552827 -0.00618587
0.0532936 0.00719838 0.269473 -0.0750147 -0.00618587 -0.0552827
-0.111945 0.0532936 -0.0707733 0.269473 -0.0707733 0.0532936
0.0532936 -0.0552827 -0.00618587 -0.0750147 0.269473 0.00719838
-0.0750147 -0.00618587 -0.0552827 0.0532936 0.00719838 0.269473
"
"0.5 -0.0348022 -0.176208 -0.0779791 -0.176208 -0.0348022
-0.0935835 0.5 0 -0.176208 -0.0826247 -0.147584
-0.176208 0 0.5 -0.0935835 -0.147584 -0.0826247
-0.0779791 -0.176208 -0.0348022 0.5 -0.0348022 -0.176208
-0.176208 -0.0826247 -0.147584 -0.0935835 0.5 0
-0.0935835 -0.147584 -0.0826247 -0.176208 0 0.5
*/

class TestMatrix: public QObject
{
    Q_OBJECT
private slots:
    void testAdd();
    void testAdd_data();
};

void TestMatrix::testAdd()
{
    QFETCH(BemMatrix, a);
    QFETCH(BemMatrix, b);
    QFETCH(BemMatrix, result);

    QCOMPARE(a + b, result);
}

void TestMatrix::testAdd_data()
{
    double array[3 * 3] = {1, 2, 3, 4, 1, 0, 0, 0, 0};
    double res_array[3 * 3] = {2, 4, 6, 8, 2, 0, 0, 0, 0};
    BemMatrix ma(array, 3, 3);
    BemMatrix mb(array, 3, 3);
    BemMatrix mresult(res_array, 3, 3);

    QTest::addColumn<BemMatrix>("a");
    QTest::addColumn<BemMatrix>("b");
    QTest::addColumn<BemMatrix>("result");

    QTest::newRow("a + b") << ma << mb << mresult;
}


QTEST_MAIN(TestMatrix)
#include "test_matrix.moc"

