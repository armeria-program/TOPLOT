#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure 1f3r"
echo "--------------------------------------------------------------"
#valgrind ../src/toplot --pdb 1f3r.pdb || exit 1
../src/toplot --pdb 1f3r.pdb || exit 1

