#! /bin/sh

echo "--------------------------------------------------------------"
echo "TOPLOT running on test structure 1cfd"
echo "--------------------------------------------------------------"
#valgrind --leak-check=full ../src/toplot --pdb 1cfd.pdb || exit 1
../src/toplot --pdb 1cfd.pdb || exit 1

