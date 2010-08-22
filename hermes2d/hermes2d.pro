QT -= GUI
TARGET = lib/hermes2d
TEMPLATE = lib
OBJECTS_DIR = build
CONFIG = += staticlib
DEFINES += NOGLUT

DEFINES += COMMON_WITH_UMFPACK
DEFINES += COMMON_WITH_SUPERLU
DEFINES += COMPLEX=std::complex\<double\>

INCLUDEPATH += src \
        src/compat \
        src/sparselib/mv \
        src/sparselib/iml \
        src/sparselib

# sparselib
SOURCES += src/sparselib/compcol_double.cc \
        src/sparselib/coord_double.cc \
        src/sparselib/icpre_double.cc \
        src/sparselib/iohb_double.cc \
        src/sparselib/qsort_double.cc \
        src/sparselib/comprow_double.cc \
        src/sparselib/diagpre_double.cc \
        src/sparselib/ilupre_double.cc \
        src/sparselib/iotext_double.cc \
        src/sparselib/qsort_int.cc \
        src/sparselib/mv/mvblasc.cc \
        src/sparselib/mv/mvblasf.cc \
        src/sparselib/mv/mvmc.cc \
        src/sparselib/mv/mvmf.cc \
        src/sparselib/mv/mvvc.cc \
        src/sparselib/mv/mvvd.cc \
        src/sparselib/mv/mvvf.cc \
        src/sparselib/mv/mvblasd.cc \
        src/sparselib/mv/mvblasi.cc \
        src/sparselib/mv/mvmd.cc \
        src/sparselib/mv/mvmi.cc \
        src/sparselib/mv/mvvcio.cc \
        src/sparselib/mv/mvvdio.cc \
        src/sparselib/mv/mvvi.cc \
        src/sparselib/spblas/spmm.cc \
        src/sparselib/spblas/spsm.cc

# hermes
SOURCES += src/compat/fmemopen.cpp \
        src/compat/c99_functions.cpp \
        src/common.cpp \
        src/hash.cpp \
        src/mesh.cpp \
        src/regul.cpp \
        src/refmap.cpp \
        src/curved.cpp \
        src/transform.cpp \
        src/traverse.cpp \
        src/shapeset.cpp \
        # src/shapeset_hc_gradeigen.cpp \
        # src/shapeset_hc_legendre.cpp \
        # src/shapeset_l2_legendre.cpp \
        # src/shapeset_hc_gradleg.cpp \
        # src/shapeset_hd_legendre.cpp \
        src/shapeset_h1_ortho.cpp \
        src/shapeset_h1_jacobi.cpp \
        src/shapeset_h1_quad.cpp \
        src/precalc.cpp \
        src/solution.cpp \
        src/limit_order.cpp \
        src/discrete_problem.cpp \
        src/linear_problem.cpp \
        src/filter.cpp \
        src/space.cpp \
        src/space_h1.cpp \
        src/space_hcurl.cpp \
        src/space_l2.cpp \
        src/space_hdiv.cpp \
        src/linear1.cpp \
        src/linear2.cpp \
        src/linear3.cpp \
        src/graph.cpp \
        src/quad_std.cpp \
        src/qsort.cpp \
        src/norm.cpp \
        src/refinement_type.cpp \
        src/element_to_refine.cpp \
        src/ref_selectors/hcurl_proj_based_selector.cpp \
        src/ref_selectors/h1_proj_based_selector.cpp \
        src/ref_selectors/l2_proj_based_selector.cpp \
        src/ref_selectors/optimum_selector.cpp \
        src/ref_selectors/order_permutator.cpp \
        src/ref_selectors/proj_based_selector.cpp \
        src/ref_selectors/selector.cpp \
        src/adapt.cpp \
        src/matrix.cpp \
        src/matrix_old.cpp \
        src/hermes2d.cpp \
        src/weakform.cpp \
        src/solvers.cpp \
        src/python_solvers.cpp \
        src/umfpack_solver.cpp \
        src/superlu_solver.cpp \
        src/sparselib_solver.cpp \
        src/mumps_solver.cpp \
        src/precond_ml.cpp \
        src/precond_ifpack.cpp \
        src/forms.cpp \
        src/mesh_parser.cpp \
        src/mesh_lexer.cpp \
        src/exodusii.cpp \
        src/h2d_reader.cpp \
        src/views/base_view.cpp \
        src/views/mesh_view.cpp \
        src/views/order_view.cpp \
        src/views/scalar_view.cpp \
        src/views/stream_view.cpp \
        src/views/vector_base_view.cpp \
        src/views/vector_view.cpp \
        src/views/view.cpp \
        src/views/view_data.cpp \
        src/views/view_support.cpp \
        src/common_time_period.cpp \
        src/data_table.cpp

HEADERS = += src/common.h

linux-g++ {
    INCLUDEPATH += /usr/include
    INCLUDEPATH += /usr/include/suitesparse
    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
    LIBS += -ldmumps_seq
    LIBS += -llapack
}
win32-g++ {
    INCLUDEPATH += c:/qt/mingw/include
    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
}
macx-g++ {
    INCLUDEPATH += /opt/local/include
    INCLUDEPATH += /opt/local/include/ufsparse

    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
}
