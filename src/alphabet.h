/*===============================================================================
alphabet.h : topology alphabet
Copyright (C) 2006-2018 Jens Kleinjung
================================================================================*/

#ifndef ALPHABET_H
#define ALPHABET_H

#include "structure.h"

/* prototypes */
char topo_state(Str *str, int seg);
void topo_sequence(Str *str, char *topseq, char *pdbfilename);
void angle_array(Str *str, char *angleFileName);

#endif
