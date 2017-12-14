/*==============================================================================
$Id: parse_args.h,v 1.2 2006/10/17 16:18:14 jkleinj Exp $
helix.c : helix axis
==============================================================================*/

#include "helix.h"

/*____________________________________________________________________________*/
/* helix axis */
/* Using the method: Kahn P. (1989) Defining the axis of a helix. 
	Computers $\&$ Chemistry 13(3): 185-189 (1989).
	The Calpha atoms of a helix turn C1, C2, C3, C4 form two angles
	(C1,C2,C3) and (C2,C3,C4) whose bisection vectors pass through the
	helix axis (although at different height h). The cross product of the 
	bisection vectors yields the vector pointing along the axis. The projection
	of C2 and C3 onto the helix axis defines a distance 'dAxis', which is 
	needed to calculate the helix radius. */
void helix_axis(Str *str, int seg)
{
	unsigned int i, j, k, l;
	int nAxisNormal;
	int nSubAxis;
	Vec *axisNormal;
	Vec *subAxis;
	Vec tmpa1, tmpa2, tmpa3; /* temporary postition vectors */
	Vec dtmpa1, dtmpa2; /* temporary difference vectors */
	Vec stmpa; /* temporary sum vectors */
	Vec tmpb2, tmpb3, tmpb4;
	Vec dtmpb3, dtmpb4;
	Vec tmpc1, tmpc2, tmpc3, tmpc4, tmpc5, tmpc6;
	float *dAxis; /* distance between Calpha projections onto helix axis */
	float *rAxis; /* radius of helix axis */
	float avrAxis = 0; /* average rAxis */
	float t1, t2, t3; /* terms in helix radius formula */

	/*-----------------------------------------------------------------------*/
	/* construct all axis normals, which are angle bisection vectors */
	/* bisection vectors are averaged angle vectors */
	nAxisNormal = str->seg[seg][1] - 1;
	axisNormal = safe_malloc(nAxisNormal * sizeof(Vec));

	for (i = 1; i < nAxisNormal; ++ i)
	{
		j = str->seg[seg][0] + i - 1;
		k = str->seg[seg][0] + i;
		l = str->seg[seg][0] + i + 1;
		/* vectors to the three angle-constituting atoms */
		atom_to_vec(&tmpa1, &str->atom[j]);
		atom_to_vec(&tmpa2, &str->atom[k]);
		atom_to_vec(&tmpa3, &str->atom[l]);
		/* two vectors constituting the angle arms */
		dtmpa1 = diff_vec(&tmpa1, &tmpa2);
		dtmpa2 = diff_vec(&tmpa3, &tmpa2);
		/* angle-disecting vector (going through helix axis) */
		stmpa = sum_vec(&dtmpa1, &dtmpa2);
		/* normalise to length 1A */
		/* this vector is called V1 or V2 in Kahn's paper */
		axisNormal[i] = norm_vec(&stmpa);
	}

	/*-----------------------------------------------------------------------*/
	/* To calculate the central axis, subaxes are constructed from 
		vector cross products of axis normals. */
	/* number of subaxes */
	nSubAxis = nAxisNormal - 1;
	subAxis = safe_malloc(nSubAxis * sizeof(Vec));
	dAxis = safe_malloc(nSubAxis * sizeof(Vec));
	rAxis = safe_malloc(nSubAxis * sizeof(Vec));

	/* vector cross product: this is called 'V1 x V2' or 'H' in Kahn's paper */
	for (i = 1; i < nSubAxis; ++ i)
	{
		/* helix axis direction 'H' */
		subAxis[i] = vec_cro_pro(&axisNormal[i], &axisNormal[i+1]);

		/* Calpha distance projected onto helix axis */
		/* the two central alpha atoms of V1 and V2 */
		k = str->seg[seg][0] + i;
		l = str->seg[seg][0] + i + 1;
		/* their position vectors */
		atom_to_vec(&tmpb2, &str->atom[k]);
		atom_to_vec(&tmpb3, &str->atom[l]);
		/* the difference vector between them */
		dtmpb3 = diff_vec(&tmpb3, &tmpb2);
		/* their distance when projected onto the helix axis */
		/* d = (P2 - P1) H */
		dAxis[i] = vec_dot_pro(&dtmpb3, &subAxis[i]);

		/* helix radius */
		/* r = |d*H|^2 - |(P2 - P1)|^2 / (2 |(P1 - P2) V2|) */
		/* first term */
		tmpb4 = scale_vec(&subAxis[i], dAxis[i]);
		t1 = vec_len_sq(&tmpb4);
		/* second term */
		t2 = vec_len_sq(&dtmpb3);
		/* third term */
		dtmpb4 = diff_vec(&tmpb2, &tmpb3);
		t3 = vec_dot_pro(&dtmpb4, &axisNormal[i+1]);
		/* radius */
		rAxis[i] = (t2 - t1) / (2 * t3);

		avrAxis += rAxis[i];
	}

	/*-----------------------------------------------------------------------*/
	/* Subaxes are averaged to yield the helix axis. */
	print_vec(&str->axis[seg]);
	exit(1);
	str->axis[seg] = average_vec(subAxis, nSubAxis);
	/* average helix radius */
	avrAxis /= (nSubAxis - 1); 
	fprintf(stderr, "radius %f\n", avrAxis);
	copy_vec(&str->axispoint[seg][0], &str->axis[seg]);


	/*-----------------------------------------------------------------------*/
	/* N-terminus */
	/* position vector to atom (= 1) of first axisNormal */
	atom_to_vec(&tmpc1, &str->atom[str->seg[seg][0]] + 1);
	print_vec(&tmpc1);
	tmpc2 = scale_vec(&axisNormal[0], avrAxis);
	fprintf(stderr, "seg %d\n", seg);
	str->axispoint[seg][0] = sum_vec(&tmpc1, &tmpc2);
	print_vec(&str->axispoint[seg][2]);
	/* C-terminus */
	/* position vector to last atom of segment */
	atom_to_vec(&tmpc3, &str->atom[str->seg[seg][0] + nSubAxis - 1]);
	tmpc4 = scale_vec(&axisNormal[str->seg[seg][1] + nSubAxis - 1], avrAxis);
	str->axispoint[seg][1] = sum_vec(&tmpc3, &tmpc4);
	/* midpoint */
	/* half-distance between N-terminus and C-terminus */
	tmpc5 = diff_vec(&str->axispoint[seg][1], &str->axispoint[seg][0]);
	tmpc6 = scale_vec(&tmpc5, 0.5);
	str->axispoint[seg][2] = sum_vec(&str->axispoint[seg][0], &tmpc6);

	/*-----------------------------------------------------------------------*/
	free(axisNormal);
	free(subAxis);
	free(dAxis);
	free(rAxis);
}


