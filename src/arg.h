/*==============================================================================
arg.h : parse command line arguments
(C) 2006 Alsseandro Pandini and Jens Kleinjung, GNU GPL License applies
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include "toplot.h"

void usage( void );
void parse_args(int argc, char **argv, char *pdbfilename);

#endif
