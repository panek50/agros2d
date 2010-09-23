QT -= GUI
TARGET = lib/hermes2d
TEMPLATE = lib
OBJECTS_DIR = build
CONFIG = += staticlib
DEFINES += NOGLUT
DEFINES += WITH_UMFPACK

linux-g++ {
    INCLUDEPATH += /usr/include/suitesparse
    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
}
win32-g++ {
    INCLUDEPATH += c:/qt/mingw/include
    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
}
macx-g++ {
    INCLUDEPATH += /opt/local/include/ufsparse
    LIBS += -lumfpack
    LIBS += -lamd
    LIBS += -lblas
}

INCLUDEPATH += src \
        src/compat
SOURCES += src/compat/fmemopen.cpp \
        src/compat/c99_functions.cpp \
        common/callstack.cpp \
        common/error.cpp \
        common/timer.cpp \
        common/trace.cpp \
        common/utils.cpp \
        src/solver/amesos.cpp \
        src/solver/aztecoo.cpp \
        src/solver/epetra.cpp \
        src/solver/mumps.cpp \
        src/solver/nox.cpp \
        src/solver/pardiso.cpp \
        src/solver/petsc.cpp \
        src/solver/precond_ifpack.cpp \
        src/solver/precond_ml.cpp \
        src/solver/umfpack_solver.cpp \
        src/common.cpp \
        src/data_table.cpp \
        src/hash.cpp \
        src/mesh.cpp \
        src/regul.cpp \
        src/refmap.cpp \
        src/curved.cpp \
        src/transform.cpp \
        src/traverse.cpp \
        src/precalc.cpp \
        src/solution.cpp \
        src/filter.cpp \
        src/space/space.cpp \
        src/space/space_h1.cpp \
        src/space/space_hcurl.cpp \
        src/space/space_l2.cpp \
        src/space/space_hdiv.cpp \
        src/linear1.cpp \
        src/linear2.cpp \
        src/linear3.cpp \
        src/graph.cpp \
        src/quad_std.cpp \
        src/shapeset/shapeset.cpp \
        # src/shapeset/shapeset_hc_gradeigen.cpp \
        src/shapeset/shapeset_hc_legendre.cpp \
        # src/shapeset/shapeset_h1_eigen.cpp \
        src/shapeset/shapeset_h1_ortho.cpp \
        src/shapeset/shapeset_l2_legendre.cpp \
        # src/shapeset/shapeset_hc_eigen2.cpp \
        # src/shapeset/shapeset_hc_gradleg.cpp \
        src/shapeset/shapeset_hd_legendre.cpp \
        src/shapeset/shapeset_h1_jacobi.cpp \
        src/shapeset/shapeset_h1_quad.cpp \
        src/qsort.cpp \
        src/norm.cpp \
        src/limit_order.cpp \
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
        src/common_time_period.cpp \
        src/matrix.cpp \
        src/hermes2d.cpp \
        src/weakform.cpp \
        src/feproblem.cpp \
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
        src/views/view_support.cpp

HEADERS = += src/common.h
