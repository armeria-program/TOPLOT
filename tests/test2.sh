#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure 3uq4"
echo "--------------------------------------------------------------"
#valgrind ../src/toplot --pdb 3uq4.pdb || exit 1
../src/toplot --pdb 3uq4.pdb || exit 1

