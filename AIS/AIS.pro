#-------------------------------------------------
#
# Project created by QtCreator 2016-06-09T20:41:23
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = AIS
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    ais.cpp \
    alglibinternal.cpp \
    alglibmisc.cpp \
    ap.cpp \
    dataanalysis.cpp \
    diffequations.cpp \
    fasttransforms.cpp \
    integration.cpp \
    interpolation.cpp \
    linalg.cpp \
    optimization.cpp \
    solvers.cpp \
    specialfunctions.cpp \
    statistics.cpp

HEADERS += \
    ais.h \
    alglibinternal.h \
    alglibmisc.h \
    ap.h \
    dataanalysis.h \
    diffequations.h \
    fasttransforms.h \
    integration.h \
    interpolation.h \
    linalg.h \
    optimization.h \
    solvers.h \
    specialfunctions.h \
    statistics.h \
    stdafx.h

OTHER_FILES += \
    fft_synchro.txt \
    ais263.dat \
    ais535.dat \
    ais553.dat \
    ../build-AIS-Qt_4_8_5-Debug/ais138.dat
