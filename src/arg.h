/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2006-2018 Jens Kleinjung
Copyright (C) 2006 Alessandro Pandini
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include "toplot.h"

void usage( void );
void parse_args(int argc, char **argv, char *pdbfilename, int *angle);

#endif
