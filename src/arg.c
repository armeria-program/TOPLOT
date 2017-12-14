/*==============================================================================
$Id: parse_args.c,v 1.3 2006/10/25 23:28:42 jkleinj Exp $
arg.c : parse command line arguments
(C) 2006 Alsseandro Pandini and Jens Kleinjung
==============================================================================*/

#include <getopt.h>
#include "toplot.h"

/*____________________________________________________________________________*/
/* parse GA command line long_options */
void usage( void )
{
	fprintf(stdout, "\nTOPLOT: TOPology pLOT\n"
					"(C) 2006-2010 Jens Kleinjung\n");
	fprintf(stdout, "\nUsage:  \ttoplot --protein <protein file>\n");
	fprintf(stdout, "\nOptions:\t--protein\tprotein input file\n"
					"        \t--help   \tthis help output\n\n"
		);
}
void parse_args(int argc, char **argv, char *pdbfilename)
{
	int c;

	static struct option long_options[] =
	{
		{"protein", required_argument, 0, 1},
        {"help", no_argument, 0, 11},
		{0, 0, 0, 0}
	};

	if (argc < 2)
	{
		usage();
		exit(1);
	}

	while ((c = getopt_long (argc, argv, "1:11", long_options, NULL)) != -1)
	{
		switch(c)
		{
            case 1:
                strcpy(pdbfilename, optarg);
                /*fprintf(stdout, "PROTEIN set to name %s\n", pdbfilename);*/
                break;
            case 11:
				usage();
                exit(1);
			default:
				usage();
				break;	
		}
	}
}


