/*==============================================================================
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

void read_pdb(FILE *pdbInFile, Str *str)
{
	unsigned int i, j, l;
	char line[80];
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
    char stopline[80] = "";
    int stopflag = 0;

   /*____________________________________________________________________________*/
    /* initialise/allocate memory for set of (64) selected (CA) atom entries */
    str->nAtom = 0;
    str->nResidue = 0;
    str->nChain = 0;

    str->atom = safe_malloc(allocated_atom * sizeof(Atom));

    /* allocate memory for sequence residues */
    str->sequence.res = safe_malloc(allocated_residue * sizeof(char));

    /*____________________________________________________________________________*/
    /* count the number of models */
    while(fgets(line, 80, pdbInFile) != 0) {
        if (strncmp(line, "MODEL ", 6) == 0) {
            if (stopflag == 0) {
                stopflag = 1;
                continue;
            } else {
                strcpy(stopline, line);
                break;
            }
        }
    }

    /* rewind the file handle to the start */
    if (fseek(pdbInFile, 0L, SEEK_SET) != 0) {
        /* handle repositioning error */
    }

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc(allocated_residue * sizeof(char));
	str->phi  = safe_malloc(allocated_residue * sizeof(float [6]));
	str->psi  = safe_malloc(allocated_residue * sizeof(float [6]));
	str->ss   = safe_malloc(allocated_residue * sizeof(int [2]));

	/* not all PDB data types are used in this program to save resources */
    while(fgets(line, 80, pdbInFile)) { /* read line */

        /*____________________________________________________________________________*/
        /* check conditions to start assigning this entry */
        /* skip other models */
        if((strcmp(line, stopline) == 0) && (stopflag == 1))
            break;

        /* read only ATOM/HETATM records */
        if((strncmp(line, "ATOM  ", 6) != 0) && (strncmp(line, "HETATM", 6) != 0))
            continue;

        /* skip alternative locations except for location 'A' */
        if (line[16] != 32 && line[16] != 65) {
            /*fprintf(stderr, "Warning: Skipping atom %d in alternative location %c\n",
                 atoi(&line[6]), line[16]);*/
            continue;
        }

        /*____________________________________________________________________________*/
		/* atom number */
		str->atom[str->nAtom].atom_nr = atoi(&line[6]);

		/* atom name */
		for (i = 12, j = 0; i < 16; )
			str->atom[str->nAtom].atom_ne[j++] = line[i++];
		str->atom[str->nAtom].atom_ne[j] = '\0';

		/* decision about recording selected (CA) atoms */
		if	((strncmp(str->atom[str->nAtom].atom_ne, " N  ", 4) != 0) && 
			 (strncmp(str->atom[str->nAtom].atom_ne, " CA ", 4) != 0) &&
			 (strncmp(str->atom[str->nAtom].atom_ne, " C  ", 3) != 0)) { 
			continue;
		}

        /* skip alternative locations except for location 'A' */
        if (line[16] != 32 && line[16] != 65) {
            /*fprintf(stderr, "Warning: Skipping atom %d in alternative location %c\n",
                 atoi(&line[6]), line[16]);*/
            continue;
        }

		/* alternative location */
		/*str->atom[str->nAtom].alt_loc[0] = line[16];	
		str->atom[str->nAtom].alt_loc[1] = '\0';*/

		/* residue name */
		for (i = 17, j = 0; i < 20; )
			str->atom[str->nAtom].res_ne[j++] = line[i++];
		str->atom[str->nAtom].res_ne[j] = '\0';

		/* chain identifier */
		str->atom[str->nAtom].chain_id[0] = line[21];
		str->atom[str->nAtom].chain_id[1] = '\0';

		/* residue number */
		str->atom[str->nAtom].res_nr = atoi(&line[22]);

		/* skip duplicate assignments, 
			e.g. when alternative locations are specified */
		/* scan backwards for identical atom name in same residue */
		for (l = 1; l < str->nAtom; ++ l) {
			/*fprintf(stderr, "l=%d res0=%d, res1=%d\n", 
				l,  str->atom[str->nAtom].res_nr, str->atom[str->nAtom - l].res_nr);*/
			if (str->atom[str->nAtom].res_nr == str->atom[str->nAtom - l].res_nr) {
				if (strcmp(str->atom[str->nAtom].atom_ne, str->atom[str->nAtom - l].atom_ne) == 0) {
					l = -1; /* flag up duplicate assignment */
					break;
				}
			}
		}

		if (l == -1) /* skip duplicate assignment */
			continue;

		/* code for insertion of residues */
		/*str->atom[str->nAtom].icode[0] = line[26];
		str->atom[str->nAtom].icode[1] = '\0';*/

		/* coordinates */
		str->atom[str->nAtom].x = atof(&line[30]);
		str->atom[str->nAtom].y = atof(&line[38]);
		str->atom[str->nAtom].z = atof(&line[46]);

		/*printf("x %6.4f, y %6.4f, z %6.4f\n", str->atom[str->nAtom].x,
			str->atom[str->nAtom].y, str->atom[str->nAtom].z);*/

		/* occupancy */
		/*str->atom[str->nAtom].occupancy = atof(&line[54]);*/

		/* temperature factor */
		/*str->atom[str->nAtom].temp_f = atof(&line[60]);*/

		/* segment identifier */
		/*for (i = 72, j = 0; i < 76; )
			str->atom[str->nAtom].seg_id[j++] = line[i++];
		str->atom[str->nAtom].seg_id[j] = '\0';*/

		/* element */
		/*for (i = 76, j = 0; i < 78; )
			str->atom[str->nAtom].element[j++] = line[i++];
		str->atom[str->nAtom].element[j] = '\0';*/

		/* charge */
		/*for (i = 78, j = 0; i < 80; )
			str->atom[str->nAtom].charge[j++] = line[i++];
		str->atom[str->nAtom].charge[j] = '\0';*/

		/* description: everything before coordinates */
		for (i = 0, j = 0; i < 30; )
			str->atom[str->nAtom].descrip[j++] = line[i++];
		str->atom[str->nAtom].descrip[j] = '\0';

		/* assign residue to sequence and count selected (CA) atoms */
		/* transform to 1-letter code and store sequence */
		/* PHI: C N CA C */
		/* PSI:   N CA C N */
		if (strncmp(str->atom[str->nAtom].atom_ne, " N  ", 4) == 0) {
			if (str->nResidue > 0) {
				str->phi[str->nResidue - 1][4] = str->phi[str->nResidue][1];
			}
			str->phi[str->nResidue][2] = str->nAtom;
			str->psi[str->nResidue][1] = str->nAtom;
		} else if (strncmp(str->atom[str->nAtom].atom_ne, " CA ", 4) == 0) {
			str->phi[str->nResidue][3] = str->nAtom;
			str->psi[str->nResidue][2] = str->nAtom;
		} else if (strncmp(str->atom[str->nAtom].atom_ne, " C  ", 4) == 0) {
			str->phi[str->nResidue][4] = str->nAtom;
			str->psi[str->nResidue][3] = str->nAtom;
			if (str->nResidue > 0) {
				str->phi[str->nResidue][1] = str->phi[str->nResidue - 1][4];
			}

			str->sequence.res[str->nResidue] = aacode(str->atom[str->nAtom].res_ne);
			++ str->nResidue;
		}

		/* increment residues */
		if (str->nResidue == allocated_residue) {
            str->sequence.res = safe_realloc(str->sequence.res, (allocated_residue += 64) * sizeof(char));
			str->phi = safe_realloc(str->phi, allocated_residue * sizeof(float [6]));
			str->psi = safe_realloc(str->psi, allocated_residue * sizeof(float [6]));
			str->ss  = safe_realloc(str->ss,  allocated_residue * sizeof(int [2]));

		}

		/* initialise other values */
		str->atom[str->nAtom].seg = 0;

		++ str->nAtom;

		/* increment atoms */
		if (str->nAtom == allocated_atom)
		{
			allocated_atom += 64;
			str->atom = safe_realloc(str->atom, allocated_atom * sizeof(Atom));
			str->sequence.res = safe_realloc(str->sequence.res, allocated_atom * sizeof(char));
		}
	}
	str->sequence.res[str->nResidue] = '\0';
	str->sequence.length = str->nResidue;
}

