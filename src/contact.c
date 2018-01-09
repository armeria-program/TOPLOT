/*==============================================================================
contact.c : Contacts between sec.str. segments
Copyright (C) 2006-2018 Jens Kleinjung
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structure.h"
#include "contact.h"

#define min(a, b)	(((a) > (b)) ? (b) : (a))

/*____________________________________________________________________________*/
/* rmsd of atom pairs */
__inline__ static float rmsd(Atom *atom0, Atom *atom1)
{
    return sqrt(pow(atom0->x - atom1->x, 2)
              + pow(atom0->y - atom1->y, 2)
              + pow(atom0->z - atom1->z, 2));
}

/*____________________________________________________________________________*/
/* contact between segments */
void contact(Str *str)
{
	unsigned int i, j, k;
	/* distance threshold are expected distance +0.5A to account for distortions */
	float distThreshold_aa = 10.5;
	float distThreshold_ab = 8.0;
	float distThreshold_bb = 5.5;
	int contCount;
	int contThreshold = 3;
	Atom atom0, atom1;
	int ss0, ss1;

	for (i = 1; i < str->nseg; ++ i) {
		contCount = 0;
		for (j = str->seg[i][0]; j < (str->seg[i][0] + str->seg[i][1]); ++ j) {	
			atom0 = str->atom[j];
			ss0 = str->seg[i][2];
			for (k = str->seg[i - 1][0]; k < (str->seg[i - 1][0] + str->seg[i - 1][1]); ++ k) {
				if (j < k + 3) /* minimal sequence distance required */
					continue;
				atom1 = str->atom[k];
				ss1 = str->seg[i-1][2];

				/*fprintf(stderr, "%d %d, dist %f\n", 
					str->seg[i-1][2], str->seg[i][2], rmsd(&atom0, &atom1));*/

				/* minimise operations */
				if (rmsd(&atom0, &atom1) >= distThreshold_aa)
					continue;

				if ((rmsd(&atom0, &atom1) < distThreshold_aa) && ((ss0 > 0) && (ss1 > 0)))
					++ contCount;
				if ((rmsd(&atom0, &atom1) < distThreshold_ab) && (((ss0 > 0) && (ss1 == 0)) || ((ss0 == 0) && (ss1 > 0))))
					++ contCount;
				if ((rmsd(&atom0, &atom1) < distThreshold_bb) && ((ss0 == 0) && (ss1 == 0)))
					++ contCount;
			}					

			if (contCount >= contThreshold) /* shortcut */
				break;
		}	

		/* set contact flag: 0=NO, 1=YES */
		if (contCount >= contThreshold)
			str->seg[i][3] = 1;
		else
			str->seg[i][3] = 0;
	
		if (DEBUG) {
			fprintf(stderr, "%s:%d: seg %d, contCount %d\n",
				__FILE__, __LINE__, i, contCount);
		}
	}
}
