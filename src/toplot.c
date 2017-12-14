/*=============================================================================
 $Id: toplot.c,v 1.18 2007/04/03 15:50:20 jkleinj Exp $
 $Name:  $

 TOPLOT: TOPology pLOT
 ~~~~~~~~~~~~~~~~~~~~~
 Program for topology plotting and analysis
 Copyright (C) 2006 Jens Kleinjung
 GNU GPL License applies

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
==============================================================================*/

#include "arg.h"
#include "alphabet.h"
#include "contact.h"
#include "getpdbs.h"
#include "structure.h"
#include "toplot.h"

/*____________________________________________________________________________*/
/* file */
FILE *safe_open(const char *name, const char *mode)
{
    FILE *file = fopen(name, mode);
	assert(file != 0);
    return file;
}

/*____________________________________________________________________________*/
/* allocation */
static void *check_non_null(void *ptr)
{
	assert (ptr != 0);
    return ptr;
}

void *safe_malloc(size_t size)
{
    return check_non_null(malloc(size));
}

void *safe_realloc(void *ptr, size_t size)
{
    return check_non_null(realloc(ptr, size));
}

/*____________________________________________________________________________*/
void error_exit(char *err_msg)
{
	fprintf(stderr, "Exiting: Error: %s\n", err_msg);
	exit(1);
}

/*____________________________________________________________________________*/
/* print Ramachandran data: plot PHI/PSI */
void ramachandran(char *filename, Str *str)
{
	FILE *file;
	unsigned int i;

	file = safe_open(filename, "w");

    for (i = 0; i < str->natom; ++ i)
        fprintf(file, "%f\t%f\n", str->phi[i][0], str->psi[i][0]);

	fclose(file);
}

/*____________________________________________________________________________*/
/* print topology data: plot PHIt/PSIt */
void topology(char *filename, Str *str)
{
	FILE *file;
	unsigned int i;

	file = safe_open(filename, "w");

    for (i = 1; i < str->nseg; ++ i)
        fprintf(file, "%f\n", str->phit[i][0]);

	fclose(file);
}

/*____________________________________________________________________________*/
/* print topology angles */
void print_angles(Str *str, FILE *angFile)
{
	unsigned int i;
	float state;
	float angle;

    for (i = 0; i < str->nseg; ++ i)
	{
		angle = str->phit[i][0];
		state = str->phit[i][1];

		if (i == 0)
		{
			/*fprintf(stderr, "seg %d, atoms %d to %d, type %d\n", 
					i, str->seg[i][0], (str->seg[i][0] + str->seg[i][1] - 1), str->seg[i][2]);*/
			continue;
		}
	}
}

/* free memory of structure */
/*____________________________________________________________________________*/
void free_pdb(Str *str, int allatom)
{
	free(str->atom);
	free(str->sequence.res);
	if (allatom)
	{
		free(str->allatom);
		free(str->phi);
		free(str->psi);
		free(str->ss);
		free(str->seg);
		free(str->phit);
		free(str->axis);
		free(str->axispoint);
	}
}

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	/* files */
    FILE *pdbfile; char *pdbfilename;
	/*char *ramafilename = "rama.dat";*/ /* output of Ramachandran plot data*/
	/*char *topofilename = "topo.dat";*/ /* output of topology plot data*/
    /*FILE *outfile; char *outfilename = "toplot.log";*/ /* TOPLOT log file */
    /*FILE *aaAngFile; char *aaAngFilename = "aa.ang";*/ /* angles beta/beta */
    /*FILE *abAngFile; char *abAngFilename = "ab.ang";*/ /* angles alpha/beta */
    /*FILE *bbAngFile; char *bbAngFilename = "bb.ang";*/ /* angles alpha/alpha */

	/* counters */
    unsigned int i = 0; /*unsigned int j = 0; unsigned int k = 0; unsigned int l = 0;*/

	int allatom = ALLATOM;

	Str pdb; /* protein structure data structure */
	char *topseq = 0; /* topology sequence */
	char err_msg[100];

	/*____________________________________________________________________________*/
	/* parse command line arguments */
	pdbfilename = safe_malloc(200 * sizeof(char));
	parse_args(argc, argv, pdbfilename);

	/*____________________________________________________________________________*/
	/* read PDB file */
	pdbfile = safe_open(pdbfilename, "r");
	read_pdb(pdbfile, &pdb, &allatom); /* read Calpha atoms */

	if (allatom) { /* read_pdb : allatom = 0 if incomplete backbone */
		allatom = read_all_pdb(pdbfile, &pdb, &allatom); /* read all atom structure */
	} else {
		strcpy(&err_msg[0], "Protein backbone incomplete!");
		error_exit(&err_msg[0]);
	}
	fclose(pdbfile);

	/*____________________________________________________________________________*/
	/* calculate topology */

	/* angle output file */
	angFile = safe_open(aaAngFilename, "w");

	if (allatom) {
		/* dihedral angles (PHI, PSI) */
		for (i = 0; i < pdb.sequence.length; ++ i)
			phi_psi(&pdb, i);

		/* Ramachandran plot */
		/*ramachandran(ramafilename, &pdb);*/

		/* sec.str. chain segments according to PHI/PSI angles */
		ss_segments(&pdb);

		/* axis angles of sec.str. segments */
		pdb.axis = safe_malloc(pdb.nseg * sizeof(Vec));
		pdb.axispoint = safe_malloc(pdb.nseg * sizeof(Vec [3]));

		segment_angle(&pdb);

		/* topology plot */
		/*topology(topofilename, &pdb);*/

		/* contacts between segments */
		contact(&pdb);

		/* topology sequence */
		topseq = safe_malloc((pdb.nseg + 1) * sizeof(char));
		topo_sequence(&pdb, &topseq[0]);

		print_angles(&pdb, angFile);
	}

	fclose(angFile);

	/*____________________________________________________________________________*/
	/* capture errors */
	else
	{
		strcpy(&err_msg[0], "Protein backbone incomplete!");
		error_exit(&err_msg[0]);
	}
    /*____________________________________________________________________________*/
	/* free memory */
	free(topseq);
	free_pdb(&pdb, allatom);
	free(pdbfilename);

	return 0;
}

