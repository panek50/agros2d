QT -= GUI
TARGET = lib/hermes2d
TEMPLATE = lib
OBJECTS_DIR = build
CONFIG = += staticlib
DEFINES += NOGLUT

INCLUDEPATH += src \
        src/compat
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
        src/shapeset_hc_gradeigen.cpp \
        src/shapeset_hc_legendre.cpp \
        src/shapeset_h1_ortho.cpp \
        src/shapeset_l2_legendre.cpp \
        src/shapeset_hc_gradleg.cpp \
        src/shapeset_hd_legendre.cpp \
        src/shapeset_h1_jacobi.cpp \
        src/shapeset_h1_quad.cpp \
        src/precalc.cpp \
        src/solution.cpp \
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
        src/common_time_period.cpp \
        src/matrix.cpp \
        src/hermes2d.cpp \
        src/weakform.cpp \
        src/solver_nox.cpp \
        src/solver_epetra.cpp \
        src/solver_aztecoo.cpp \
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
        src/views/view_support.cpp

HEADERS = += src/common.h
