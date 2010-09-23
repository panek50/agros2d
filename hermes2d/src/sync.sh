#!/bin/sh

rm -f *.cpp
rm -f *.h

cp ../../../hermes/hermes2d/src/*.cpp .
cp ../../../hermes/hermes2d/src/*.h .
cp ../../../hermes/hermes2d/src/views/*.cpp ./views
cp ../../../hermes/hermes2d/src/views/*.h ./views
cp ../../../hermes/hermes2d/src/ref_selectors/*.cpp ./ref_selectors
cp ../../../hermes/hermes2d/src/ref_selectors/*.h ./ref_selectors
cp ../../../hermes/hermes2d/src/gen/*.cpp ./gen
cp ../../../hermes/hermes2d/src/solver/*.cpp ./solver
cp ../../../hermes/hermes2d/src/solver/*.h ./solver
cp ../../../hermes/hermes2d/src/shapeset/*.cpp ./shapeset
cp ../../../hermes/hermes2d/src/shapeset/*.h ./shapeset
cp ../../../hermes/hermes2d/src/space/*.cpp ./space
cp ../../../hermes/hermes2d/src/space/*.h ./space
cp ../../../hermes/hermes2d/src/compat/*.cpp ./compat
cp ../../../hermes/hermes2d/src/compat/*.h ./compat

cp ../../../hermes/hermes2d/common/*.cpp ../common
cp ../../../hermes/hermes2d/common/*.h ../common

#cp ../../../hermes/hermes2d/hermes_common/common_time_period.cpp .
#cp ../../../hermes/hermes2d/hermes_common/common_time_period.h .
#cp ../../../hermes/hermes2d/hermes_common/tuple.h .

#cp ../../../hermes/hermes2d/hermes_common/matrix.cpp .
#cp ../../../hermes/hermes2d/hermes_common/matrix.h .
#cp ../../../hermes/hermes2d/hermes_common/solvers.cpp .
#cp ../../../hermes/hermes2d/hermes_common/solvers.h .
#cp ../../../hermes/hermes2d/hermes_common/superlu_solver.cpp .
#cp ../../../hermes/hermes2d/hermes_common/umfpack_solver.cpp .
#cp ../../../hermes/hermes2d/hermes_common/sparselib_solver.cpp .
#cp ../../../hermes/hermes2d/hermes_common/python_solvers.cpp .
#cp ../../../hermes/hermes2d/hermes_common/tuple.h .
