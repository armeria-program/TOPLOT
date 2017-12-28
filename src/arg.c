/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2006-2017 Jens Kleinjung
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
			"\t  --help\n\n");

	fprintf(stdout,	"(C) 2006-2017 Jens Kleinjung\n\n");

	exit(0);
}

/*____________________________________________________________________________*/
void parse_args(int argc, char **argv, char *pdbFileName)
{
	int c;

	static struct option long_options[] =
	{
		{"pdb", required_argument, 0, 1},
        {"help", no_argument, 0, 11},
		{0, 0, 0, 0}
	};

	if (argc < 2)
	{
		usage();
	}

	while ((c = getopt_long (argc, argv, "1:11", long_options, NULL)) != -1)
	{
		switch(c)
		{
            case 1:
                strcpy(pdbFileName, optarg);
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


