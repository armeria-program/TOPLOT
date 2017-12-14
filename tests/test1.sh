#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure"
echo "--------------------------------------------------------------"
../src/tplot --pdb 1f3r.pdb || exit 1

