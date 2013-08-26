# agros2d - hp-FEM multiphysics application based on Hermes2D library
OBJECTS_DIR = build
MOC_DIR = build

TEMPLATE = lib

DEFINES += BEM_LIBRARY

SOURCES += bem.cpp
SOURCES += algebra.cpp

HEADERS += bem.h
HEADERS += algebra.h

INCLUDEPATH += ../hermes2d/include
INCLUDEPATH += ../hermes2d/include/mesh/
INCLUDEPATH += ../hermes_common/include
INCLUDEPATH += ../pythonlab-library
INCLUDEPATH += ../util
INCLUDEPATH += ../agros2d-library/

LIBS += -lagros_3rdparty
LIBS += -lagros_util
LIBS += -lagros_pythonlab_library
LIBS += -lagros2d_hermes2d
LIBS += -lagros2d_library

linux-g++|linux-g++-64|linux-g++-32|linux-clang {
    # DEFINES += WITH_UNITY
    TARGET = ../libs/bem

    CONFIG += warn_off
    # QMAKE_CXXFLAGS += -Wun
}

include(../agros2d.pri)
include(../agros2d_version.pri)
