from unittest import TestLoader, TestResult
import test_script_geometry

result = TestResult()

""" geometry """
TestLoader().loadTestsFromTestCase(test_script_geometry.TestGeometry).run(result)
print('Basic geometry test: {0}'.format(result.wasSuccessful()))

TestLoader().loadTestsFromTestCase(test_script_geometry.TestGeometryTransformations).run(result)
print('Geometry transformations test: {0}'.format(result.wasSuccessful()))

if (result.wasSuccessful() != True):
    print(result.errors)