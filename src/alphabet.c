/*==============================================================================
alphabet.c : topology alphabet
Copyright (C) 2006-2018 Jens Kleinjung
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alphabet.h"

/*___________________________________________________________________________*/
/* returns topology state */
/* 3 sec.str. states, 7 angle states and 2 contact states */
char topo_state(Str *str, int seg)
{
	int state = 0;
	int ss1, ss2;

	ss1 = str->seg[seg-1][2];
	ss2 = str->seg[seg][2];	

	/* 3 sec.str. states: 1. beta/beta, 2. beta/alpha or alpha/beta, 3. alpha/alpha */
	/* 1. beta/beta [A:C] */
	if (ss1 == 0 && ss2 == 0) {
		state = 65;
	}
	/* 2. beta/alpha [D:F] */
	if (ss1 == 0 && ss2 > 0) {
		state = 68;
	}
	/* 2. alpha/beta [G:I] */
	if (ss1 > 0 && ss2 == 0) {
		state = 71;
	}
	/* 3. alpha/alpha [J:L] */
	if (ss1 > 0 && ss2 > 0) {
		state = 74;
	}
	/* 3 angle states: [0:60[, [60:120[, [120:180[ */
	state += floor((str->phit[seg][0]) / 60);
	if ((str->phit[seg][0] / 60) == 3) {
		-- state;
	}
	/* 2 contact states: yes or no */
	if (str->seg[seg][3] == 0) {
		state += 32; /* convert to lower case */
	}
	/*fprintf(stderr, "%s:%d: ss1 %d, ss2 %d state %d\n",
							__FILE__, __LINE__, ss1, ss2, state);*/

	assert((state >= 65) && (state <= 122));

	return (char)state;
}

/*___________________________________________________________________________*/
/* generate topology sequence */
void topo_sequence(Str *str, char *topseq, char *pdbFileName)
{
	unsigned int seg = 0;
	unsigned int nChain = 0;
	char state;
	char outFileName[256];
	FILE *outFile;

	if (str->nseg >= 2) {
		/* for all other segments */
		for (seg = 1; seg < str->nseg; ++ seg) {
			/* determine topology state */
			state = topo_state(str, seg);
			str->phit[seg][1] = (float)state;
			/* add residue state to topology string */
			if (seg == 1) {
				topseq[0] = tolower(state);
			}
			/* insert '|'character into topology string
				if chain changes between segments */
			if (strcmp(str->atom[str->seg[seg][0]].chain_id, str->atom[str->seg[seg - 1][0]].chain_id) == 0) {
				topseq[seg + nChain] = state;
			} else {
				topseq[seg + nChain] = '|';
				++ nChain;
				topseq[seg + nChain] = state;
			}
		}
	}
	topseq[seg] = '\0';	

	sprintf(&(outFileName[0]), "%s%s", pdbFileName, ".fastt"); 
	outFile = safe_open(outFileName, "w");
	if (strlen(topseq) > 1)
		fprintf(outFile, ">%s\n%s\n", pdbFileName, topseq);
	else
		fprintf(stderr, "Topology file of PDB structure %s is void!\n",
							pdbFileName);
	fclose(outFile);
}

/*___________________________________________________________________________*/
/* generate angle array */
void angle_array(Str *str, char *angleFileName)
{
	unsigned int seg = 0;
	FILE *angleFile; 

	angleFile = safe_open(angleFileName, "w");

	for (seg = 0; seg < str->nseg; ++ seg) {
		fprintf(angleFile, "%3.0f", str->phit[seg][0]);
		if (seg < str->nseg - 1) {
			fprintf(angleFile, "\t");
		} else {
			fprintf(angleFile, "\n");
		}		
	}
}

