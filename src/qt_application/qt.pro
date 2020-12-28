#-------------------------------------------------
#
# Project created by QtCreator 2017-04-05T14:30:09
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = qtIsing
TEMPLATE = app


SOURCES += qt.cpp\
           app.cpp \
           ising_window.cpp \
           ../model/ising/ising.cpp \
           ../simulation.cpp \
           ../grid/grid.cpp

HEADERS  += ../../include/qt_application/app.h \
            ../../include/qt_application/ising_window.h \
            ../../include/model/ising/ising.h \
            ../../include/simulation.h \
            ../../include/grid/grid.h

CONFIG += c++17
unix:INCLUDEPATH += "../../include"

QMAKE_CXXFLAGS_RELEASE *= -O3
QMAKE_CXXFLAGS_RELEASE *= -Wall -Wextra
