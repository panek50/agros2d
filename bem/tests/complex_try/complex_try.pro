#-------------------------------------------------
#
# Project created by QtCreator 2014-01-21T09:47:19
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = complex_try
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    ../../bem_matrix.cpp

HEADERS += \
    ../../bem_matrix.h

LIBS += -llapacke
LIBS += -lboost_math_c99
LIBS += -lblas
