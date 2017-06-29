#!/bin/sh

CXX=${1}
CXXFLAGS=${2}
CPPFLAGS=${3}
LDFLAGS=${4}
LIBS=${5}

MPICXX=${6}
MPICXXFLAGS=${7}
MPICPPFLAGS=${8}
MPILDFLAGS=${9}
MPILIBS=${10}

PYPATH=${11}

PLUGCOMP=${12}
PLUGLINK=${13}
PLUGEXT=${14}

LIBEXT=${15}

sed \
-e "s#@MPICXX@#${MPICXX}#g" \
-e "s#@MPICXXFLAGS@#${MPICXXFLAGS}#g" \
-e "s#@MPICPPFLAGS@#${MPICPPFLAGS}#g" \
-e "s#@MPILDFLAGS@#${MPILDFLAGS}#g" \
-e "s#@MPILIBS@#${MPILIBS}#g" \
-e "s#@CXX@#${CXX}#g" \
-e "s#@CXXFLAGS@#${CXXFLAGS}#g" \
-e "s#@CPPFLAGS@#${CPPFLAGS}#g" \
-e "s#@LDFLAGS@#${LDFLAGS}#g" \
-e "s#@LIBS@#${LIBS}#g" \
-e "s#@PYPATH@#${PYPATH}#g" \
-e "s#@PLUGCOMP@#${PLUGCOMP}#g" \
-e "s#@PLUGLINK@#${PLUGLINK}#g" \
-e "s#@PLUGEXT@#${PLUGEXT}#g" \
-e "s#@LIBEXT@#${LIBEXT}#g" \
harpconfig.in > harpconfig

