#!/bin/sh

CXX=${1}
CXXFLAGS=${2}
CPPFLAGS=${3}
LDFLAGS=${4}
LIBS=${5}
MPICXX=${6}
MPILIBS=${7}
PYPATH=${8}

sed -e "s#@CXX@#${CXX}#g" -e "s#@CXXFLAGS@#${CXXFLAGS}#g" -e "s#@CPPFLAGS@#${CPPFLAGS}#g" -e "s#@LDFLAGS@#${LDFLAGS}#g" -e "s#@LIBS@#${LIBS}#g" -e "s#@MPICXX@#${MPICXX}#g" -e "s#@MPILIBS@#${MPILIBS}#g" -e "s#@PYPATH@#${PYPATH}#g" harpconfig.in > harpconfig

