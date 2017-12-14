/*==============================================================================
 $Id: getpdbs.c,v 1.5 2007/04/03 15:50:20 jkleinj Exp $ 
 getpdbs.c : Routines for reading PDB structures
 Copyright (C) 2004 Jens Kleinjung
 GNU GPL Licence applies
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getpdbs.h"
#include "structure.h"

/*____________________________________________________________________________*/
/* amino acid 3-letter to 1-letter code conversion */
__inline__ static char aacode(char *code3)
{
	unsigned int i;
	char aa1 = 'O';
	char *aa3[] = {"ALA","---","CYS","ASP","GLU","PHE","GLY","HIS","ILE","---","LYS","LEU","MET","ASN","---","PRO","GLN","ARG","SER","THR","---","VAL","TRP","---","TYR","---"};

	for (i = 0; i < 26; ++i)
		if (strncmp(code3, aa3[i], 3) == 0)
		{
			aa1 = i + 65;
			break;
		}

	if (aa1 == 'O')
	{
		fprintf(stderr, "Exiting: Unkown residue type '%3s' in protein structure. Changing to standard residue \"GLY\".\n", code3);
		exit(2);
		aa1 = 'G';
	}
	
	return aa1;
}

/*____________________________________________________________________________*/
/* read PDB file */

/*
Definition of PDB format:
http://pdbbeta.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html
The ATOM record:

COLUMNS        DATA TYPE       FIELD         DEFINITION
________________________________________________________________________________-
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X
                                               in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
                                               in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
                                               in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier,
                                               left-justified.
77 - 78        LString(2)      element       Element symbol,
                                               right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
*/

void read_pdb(FILE *pdbfile, Str *str, int *allatom)
{
	unsigned int i, j, l;
	unsigned int k = 0;
	char line[80];
	unsigned int allocated = 64;

	/* initialise/allocate memory for set of (64) selected (CA) atom entries */
	str->natom = 0;
	str->atom = safe_malloc(allocated * sizeof(Atom));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc(allocated * sizeof(char));
	if (allatom) 
	{
		str->phi  = safe_malloc(allocated * sizeof(float [6]));
		str->psi  = safe_malloc(allocated * sizeof(float [6]));
		str->ss   = safe_malloc(allocated * sizeof(int [2]));
	}

	/* not all PDB data types are used in this program to save resources */
    while(! feof(pdbfile))
    {
        fgets(line, 80, pdbfile); /* read line */

		/* record name */
		if(strncmp(line, "ATOM  ", 6) != 0)
			continue;

		/* atom number */
		str->atom[str->natom].atom_nr = atoi(&line[6]);

		/* atom name */
		for (i = 12, j = 0; i < 16; )
			str->atom[str->natom].atom_ne[j++] = line[i++];
		str->atom[str->natom].atom_ne[j] = '\0';

		/* decision about recording selected (CA) atoms */
		if(strncmp(str->atom[str->natom].atom_ne, " CA ", 4) != 0)
			continue;

		/* alternative location */
		/*str->atom[str->natom].alt_loc[0] = line[16];	
		str->atom[str->natom].alt_loc[1] = '\0';*/

		/* residue name */
		for (i = 17, j = 0; i < 20; )
			str->atom[str->natom].res_ne[j++] = line[i++];
		str->atom[str->natom].res_ne[j] = '\0';

		/* chain identifier */
		str->atom[str->natom].chain_id[0] = line[21];
		str->atom[str->natom].chain_id[1] = '\0';

		/* residue number */
		str->atom[str->natom].res_nr = atoi(&line[22]);

		/* skip duplicate assignments, 
		e.g. when alternative locations are specified */
		/* scan backwards for identical atom name in same residue */
		for (l = 1; l < str->natom; ++ l)
		{
			/*fprintf(stderr, "l=%d res0=%d, res1=%d\n", 
				l,  str->atom[str->natom].res_nr, str->atom[str->natom - l].res_nr);*/
			if (str->atom[str->natom].res_nr == str->atom[str->natom - l].res_nr)
			{
				if (strcmp(str->atom[str->natom].atom_ne, str->atom[str->natom - l].atom_ne) == 0)
				{
					l = -1; /* flag up duplicate assignment */
					break;
				}
			}
		}

		if (l == -1) /* skip duplicate assignment */
			continue;

		/* code for insertion of residues */
		/*str->atom[str->natom].icode[0] = line[26];
		str->atom[str->natom].icode[1] = '\0';*/

		/* coordinates */
		str->atom[str->natom].x = atof(&line[30]);
		str->atom[str->natom].y = atof(&line[38]);
		str->atom[str->natom].z = atof(&line[46]);

		/*printf("x %6.4f, y %6.4f, z %6.4f\n", str->atom[str->natom].x,
			str->atom[str->natom].y, str->atom[str->natom].z);*/

		/* occupancy */
		/*str->atom[str->natom].occupancy = atof(&line[54]);*/

		/* temperature factor */
		/*str->atom[str->natom].temp_f = atof(&line[60]);*/

		/* segment identifier */
		/*for (i = 72, j = 0; i < 76; )
			str->atom[str->natom].seg_id[j++] = line[i++];
		str->atom[str->natom].seg_id[j] = '\0';*/

		/* element */
		/*for (i = 76, j = 0; i < 78; )
			str->atom[str->natom].element[j++] = line[i++];
		str->atom[str->natom].element[j] = '\0';*/

		/* charge */
		/*for (i = 78, j = 0; i < 80; )
			str->atom[str->natom].charge[j++] = line[i++];
		str->atom[str->natom].charge[j] = '\0';*/

		/* description: everything before coordinates */
		for (i = 0, j = 0; i < 30; )
			str->atom[str->natom].descrip[j++] = line[i++];
		str->atom[str->natom].descrip[j] = '\0';

		/* assign residue to sequence and count selected (CA) atoms */
		/* transform to 1-letter code and store sequence */
		str->sequence.res[k++] = aacode(str->atom[str->natom].res_ne);

		/* initialise other values */
		str->atom[str->natom].seg = 0;

		++ str->natom;

		/* allocate more memory if needed */
		if (str->natom == allocated)
		{
			allocated += 64;
			str->atom = safe_realloc(str->atom, allocated * sizeof(Atom));
			str->sequence.res = safe_realloc(str->sequence.res, allocated * sizeof(char));
			if (allatom) 
			{
				str->phi = safe_realloc(str->phi, allocated * sizeof(float [6]));
				str->psi = safe_realloc(str->psi, allocated * sizeof(float [6]));
				str->ss  = safe_realloc(str->ss,  allocated * sizeof(int [2]));
			}
		}
	}
	str->sequence.res[k] = '\0';
	str->sequence.length = str->natom;
}

/*____________________________________________________________________________*/
/* read all atoms from PDB file */
int read_all_pdb(FILE *pdbfile, Str *str, int *allatom)
{
	unsigned int i, j, l;
	char line[80];
	unsigned int allocated = 64;
	unsigned int res = 0; /* (internal) residue number */
	unsigned int bb_count = 0; /* count backbone atoms : completion check */

	/* initialise/allocate memory for set of (64) all atom entries */
	str->nallatom = 0;
	str->allatom = safe_malloc(allocated * sizeof(Atom));

	/* reset file pointer */
	fseek(pdbfile, 0, 0);

	/* for comments see routine above */
    while(! feof(pdbfile))
    {
        fgets(line, 80, pdbfile); /* read line */

		if(strncmp(line, "ATOM  ", 6) != 0)
			continue;

		str->allatom[str->nallatom].atom_nr = atoi(&line[6]);

		/* alternative location */
		/*str->allatom[str->nallatom].alt_loc[0] = line[16];
		str->allatom[str->nallatom].alt_loc[1] = '\0';
		if ((line[16] != ' ') && (line[16] != 'A'))
			continue;*/

		for (i = 12, j = 0; i < 16; )
			str->allatom[str->nallatom].atom_ne[j++] = line[i++];
		str->allatom[str->nallatom].atom_ne[j] = '\0';

		for (i = 17, j = 0; i < 20; )
			str->allatom[str->nallatom].res_ne[j++] = line[i++];
		str->allatom[str->nallatom].res_ne[j] = '\0';

		str->allatom[str->nallatom].chain_id[0] = line[21];
		str->allatom[str->nallatom].chain_id[1] = '\0';

		str->allatom[str->nallatom].res_nr = atoi(&line[22]);

		/* skip duplicate assignments, 
		e.g. when alternative locations are specified */
		/* scan backwards for identical atom name in same residue */
		for (l = 1; l <= str->nallatom; ++ l)
		{
			/*fprintf(stderr, "l=%d allres0=%d, allres1=%d\n", 
				l,  str->allatom[str->nallatom].res_nr, str->allatom[str->nallatom - l].res_nr);*/
			if (str->allatom[str->nallatom].res_nr == str->allatom[str->nallatom - l].res_nr)
			{
				if (strcmp(str->allatom[str->nallatom].atom_ne, str->allatom[str->nallatom - l].atom_ne) == 0)
				{
					l = -1; /* flag up duplicate assignment */
					break;
				}
			}
		}

		if (l == -1) /* skip duplicate assignment */
			continue;

		str->allatom[str->nallatom].x = atof(&line[30]);
		str->allatom[str->nallatom].y = atof(&line[38]);
		str->allatom[str->nallatom].z = atof(&line[46]);

		for (i = 0, j = 0; i < 30; )
			str->allatom[str->nallatom].descrip[j++] = line[i++];
		str->allatom[str->nallatom].descrip[j] = '\0';

		/* recording atom numbers of atoms constituting PHI/PSI angles */
		/* N */
		if(strncmp(str->allatom[str->nallatom].atom_ne, " N  ", 4) == 0)
		{
			str->phi[res][2] = (float)str->nallatom;
			str->psi[res][1] = (float)str->nallatom;
			if (res > 0) 
				str->psi[res - 1][4] = (float)str->nallatom;
			++ bb_count;
		}
		/* CA */
		if(strncmp(str->allatom[str->nallatom].atom_ne, " CA ", 4) == 0)
		{
            str->phi[res][3] = (float)str->nallatom;
            str->psi[res][2] = (float)str->nallatom;
			/* preset dihedral angles */
			str->phi[res][0] = 999.;
			str->psi[res][0] = 999.;
			++ bb_count;
		}
		/* C */
		if(strncmp(str->allatom[str->nallatom].atom_ne, " C  ", 4) == 0)
		{
            str->phi[res][4] = (float)str->nallatom;
            str->psi[res][3] = (float)str->nallatom;
			if (res < str->natom - 1)
				str->phi[res + 1][1] = (float)str->nallatom;
			++ bb_count;
			/* check for backbone completeness */
			if ((bb_count % 3) == 0)
				++res;
			else
				return 0; /* resturn allatom = 0 */
		}

		++ str->nallatom;

		if (str->nallatom == allocated)
		{
			allocated += 64;
			str->allatom = safe_realloc(str->allatom, allocated * sizeof(Atom));
		}

	}

	/*if (allatom)
		fprintf(stdout, "\tall atom number = %5d\n", str->nallatom);*/

	return 1; /* return allatom = 1 */
}


