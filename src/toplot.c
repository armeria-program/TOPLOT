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
	fprintf(stdout, "Reading input structure %s\n", pdbFileName);
	pdbfile = safe_open(pdbFileName, "r");
	read_pdb(pdbfile, &pdb);
	fclose(pdbfile);

	backbone_completeness(&pdb, 0);

	/*____________________________________________________________________________*/
	/* calculate topology */
	fprintf(stdout, "Assessing topology\n");
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

	/*____________________________________________________________________________*/
	/* topology sequence */
	fprintf(stdout, "Writing topology string %s.fastt\n", pdbFileName);
	topseq = safe_malloc((pdb.nseg + 1) * sizeof(char));
	topo_sequence(&pdb, &(topseq[0]), &(pdbFileName[0]));

    /*____________________________________________________________________________*/
	/* free memory */
	free(pdb.sequence.res);
	free(pdb.atom);
	free(pdb.phi);
	free(pdb.psi);
	free(pdb.ss);
	free(pdb.seg);
	free(pdb.phit);
	free(pdb.axis);
	free(pdb.axispoint);
	free(pdbFileName);
	free(topseq);

	return 0;
}

