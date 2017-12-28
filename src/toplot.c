/*=============================================================================
TOPLOT: TOPology pLOT
Encode protein topology into string
Copyright (C) 2006-2017 Jens Kleinjung

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
void ramachandran(char *fileName, Str *str)
{
	FILE *file;
	unsigned int i;

	file = safe_open(fileName, "w");

    for (i = 0; i < str->nAtom; ++ i)
        fprintf(file, "%f\t%f\n", str->phi[i][0], str->psi[i][0]);

	fclose(file);
}

/*____________________________________________________________________________*/
/* print topology data: plot PHIt/PSIt */
void topology(char *fileName, Str *str)
{
	FILE *file;
	unsigned int i;

	file = safe_open(fileName, "w");

    for (i = 1; i < str->nseg; ++ i)
        fprintf(file, "%f\n", str->phit[i][0]);

	fclose(file);
}

/*____________________________________________________________________________*/
/* print topology angles */
void print_angles(Str *str, FILE *angFile)
{
	unsigned int i;
	//float state;
	//float angle;

    for (i = 0; i < str->nseg; ++ i) {
		//angle = str->phit[i][0];
		//state = str->phit[i][1];

		if (i == 0) {
			/*fprintf(stderr, "seg %d, atoms %d to %d, type %d\n", 
					i, str->seg[i][0], (str->seg[i][0] + str->seg[i][1] - 1), str->seg[i][2]);*/
			continue;
		}
	}
}

/* free memory of structure */
/*____________________________________________________________________________*/
void free_pdb(Str *str)
{
	free(str->atom);
	free(str->sequence.res);
	free(str->phi);
	free(str->psi);
	free(str->ss);
	free(str->seg);
	free(str->phit);
	free(str->axis);
	free(str->axispoint);
}

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	/* files */
    FILE *pdbfile; char *pdbFileName;
    /*FILE *angFile; char *angfileName = "ang.dat";*/

	/* counters */
    unsigned int i = 0;

	Str pdb; /* protein structure data structure */
	char *topseq = 0; /* topology sequence */

	/*____________________________________________________________________________*/
	/* parse command line arguments */
	pdbFileName = safe_malloc(200 * sizeof(char));
	parse_args(argc, argv, pdbFileName);

	/*____________________________________________________________________________*/
	/* read PDB file */
	pdbfile = safe_open(pdbFileName, "r");
	read_pdb(pdbfile, &pdb);
	fclose(pdbfile);

	/*____________________________________________________________________________*/
	/* calculate topology */
	/* dihedral angles (PHI, PSI) */
	for (i = 0; i < pdb.sequence.length; ++ i) {
		phi_psi(&pdb, i);
	}

	/* sec.str. chain segments according to PHI/PSI angles */
	ss_segments(&pdb);

	/* axis angles of sec.str. segments */
	pdb.axis = safe_malloc(pdb.nseg * sizeof(Vec));
	pdb.axispoint = safe_malloc(pdb.nseg * sizeof(Vec [3]));

	segment_angle(&pdb);

	/* contacts between segments */
	contact(&pdb);

	/* topology sequence */
	topseq = safe_malloc((pdb.nseg + 1) * sizeof(char));
	topo_sequence(&pdb, &(topseq[0]), &(pdbFileName[0]));

	/* angle output file */
	/*
	angFile = safe_open(angfileName, "w");
	print_angles(&pdb, angFile);
	fclose(angFile);
	*/

    /*____________________________________________________________________________*/
	/* free memory */
	free(topseq);
	free_pdb(&pdb);
	free(pdbFileName);

	return 0;
}

