/*==============================================================================
 $Id: alphabet.c,v 1.13 2007/04/03 15:50:20 jkleinj Exp $ 
 alphabet.c : topology alphabet
 Copyright (C) 2006 Jens Kleinjung
 GNU GPL License applies
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
	if (ss1 == 0 && ss2 == 0)
		state = 65;
	/* 2. beta/alpha [D:F] */
	if (ss1 == 0 && ss2 > 0)
		state = 68;
	/* 2. alpha/beta [G:I] */
	if (ss1 > 0 && ss2 == 0)
		state = 71;
	/* 3. alpha/alpha [J:L] */
	if (ss1 > 0 && ss2 > 0)
		state = 74;

	/* 3 angle states: [0:60[, [60:120[, [120:180[ */
	state += floor((str->phit[seg][0]) / 60);
	if ((str->phit[seg][0] / 60) == 3) 
		-- state;

	/* 2 contact states: yes or no */
	if (str->seg[seg][3] == 0)
		state += 32; /* convert to lower case */

	/*fprintf(stderr, "ss1 %d, ss2 %d state %d\n", ss1, ss2, state);*/

	assert((state >= 65) && (state <= 122));

	return (char)state;
}

/*___________________________________________________________________________*/
/* generate topology sequence */
void topo_sequence(Str *str, char *topseq, char *pdbfilename)
{
	unsigned int seg = 0;
	char state;

	if (str->nseg >= 2)
	{
		for (seg = 1; seg < str->nseg; ++ seg) /* for all other segments */
		{
			/* determine topology state */
			state = topo_state(str, seg);
			str->phit[seg][1] = (float)state;
			/* add residue state to topology string */
			if (seg == 1)
				topseq[0] = tolower(state);
			topseq[seg] = state;
		}
	}
	topseq[seg] = '\0';	
	if (strlen(topseq) > 1)
		fprintf(stdout, ">%s\n%s\n", pdbfilename, topseq);
	else
		fprintf(stdout, "?\n");
}
