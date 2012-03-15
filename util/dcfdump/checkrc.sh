#!/bin/tcsh

echo "Running dcfdump ..."
dcfdump $*
set rc = $?

echo "-------------------------"
echo "Return code = $rc"

