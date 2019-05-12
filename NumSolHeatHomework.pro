#-------------------------------------------------
#
# Project created by QtCreator 2019-05-09T16:49:50
#
#-------------------------------------------------

QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = NumSolHeatHomework
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11     # Delete 'console'
#CONFIG -= app_bundle        # Delete for using Qt libraries
#CONFIG -= qt                # .............................

SOURCES += \
        main.cpp \
    types.cpp \
    implicit_diff_scheme_cyl.cpp \
    plotter.cpp \
    mainwindow.cpp

HEADERS += \
    types.h \
    implicit_diff_scheme_cyl.h \
    plotter.h \
    err.h \
    mainwindow.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

unix:!macx: LIBS += -L$$PWD/../../../../usr/local/lib/ -lgsl
unix:!macx: LIBS += -L$$PWD/../../../../usr/local/lib/ -lgslcblas

INCLUDEPATH += $$PWD/../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../usr/local/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../usr/local/lib/libgsl.a
unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../usr/local/lib/libgslcblas.a

FORMS += \
    mainwindow.ui
