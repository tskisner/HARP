#!/bin/sh

INTERP=$1
EXEC=$2
IN=$3

cat ${IN} | sed -e "s#@PYTHON@#${INTERP}#g" > ${EXEC}; chmod +x ${EXEC}

