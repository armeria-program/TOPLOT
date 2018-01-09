/*===============================================================================
getpdbs.h : Read PDB structures 
Copyright (C) 2004-2018 Jens Kleinjung
================================================================================*/

#ifndef GETPDBS_H
#define GETPDBS_H

#include "structure.h"

void read_pdb(FILE *pdbInFile, Str *str);
void backbone_completeness(Str *str, int ca);

#endif
