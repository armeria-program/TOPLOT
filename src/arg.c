/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2006-2018 Jens Kleinjung
Copyright (C) 2006 Alessandro Pandini
==============================================================================*/

#include <getopt.h>
#include "toplot.h"

/*____________________________________________________________________________*/
/* parse GA command line long_options */
void usage( void )
{
	fprintf(stdout, "\nTOPLOT: TOPology pLOT\n");

	fprintf(stdout, "\ntoplot [--pdb ...] [OPTIONS ...]\n"
			"\tOPTIONS\n"
			"\t  --pdb <PDB input>\n"
			"\t  --angle\n"
			"\t  --help\n\n");

	fprintf(stdout,	"(C) 2006-2018 Jens Kleinjung\n\n");

	exit(0);
}

/*____________________________________________________________________________*/
void parse_args(int argc, char **argv, char *pdbFileName, int *angle)
{
	int c;

	static struct option long_options[] =
	{
		{"pdb", required_argument, 0, 1},
        	{"angle", no_argument, 0, 2},
        	{"help", no_argument, 0, 11},
		{0, 0, 0, 0}
	};

	if (argc < 2)
	{
		usage();
	}

	while ((c = getopt_long (argc, argv, "1:2 11", long_options, NULL)) != -1)
	{
		switch(c)
		{
            case 1:
                strcpy(pdbFileName, optarg);
                break;
            case 2:
                *angle = 1;
                break;
            case 11:
				usage();
                exit(0);
			default:
				usage();
				break;	
		}
	}
}


