TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    test_header.h

INCLUDEPATH += /usr/local/include
LIB += -L/usr/local/lib
LIBS += -llapack -lblas -larmadillo
