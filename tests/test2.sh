#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure 3UQ4"
echo "--------------------------------------------------------------"
#valgrind ../src/toplot --pdb 3UQ4.pdb || exit 1
../src/toplot --pdb 3UQ4.pdb || exit 1

