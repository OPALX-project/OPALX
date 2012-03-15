#!/bin/tcsh

echo "Running dcfdiff ..."
dcfdiff $*
set rc = $?

echo "-------------------------"
echo "Return code = $rc"

