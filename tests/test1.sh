#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure"
echo "--------------------------------------------------------------"
../src/toplot --pdb 1cfd.pdb || exit 1

