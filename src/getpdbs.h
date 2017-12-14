/*===============================================================================
 $Id: getpdbs.h,v 1.3 2006/10/27 09:19:00 jkleinj Exp $
 getpdbs.h : Read PDB structures 
 Copyright (C) 2004 Jens Kleinjung
 GNU GPL License applies
================================================================================*/

#ifndef GETPDBS_H
#define GETPDBS_H

#include "structure.h"

void read_pdb(FILE *pdbfile, Str *str, int *allatom);
int read_all_pdb(FILE *pdbfile, Str *str, int *allatom);

#endif
