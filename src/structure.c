/*==============================================================================
 structure.c : Routines for structure operations
 Copyright (C) 2004-2006 Jens Kleinjung
 GNU GPL License applies
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structure.h"

/*____________________________________________________________________________*/
/* define */
#define PI (4.0 * atan(1.0))
#define AN (180 / PI)

/*____________________________________________________________________________*/
/* vector length |v| */
float vec_len(Vec *v)
{
    return sqrt(pow(v->x, 2) + pow(v->y, 2) + pow(v->z, 2));
}

/*____________________________________________________________________________*/
/* vector square length |v|^2 */
float vec_len_sq(Vec *v)
{
    return (pow(v->x, 2) + pow(v->y, 2) + pow(v->z, 2));
}

/*____________________________________________________________________________*/
/* normalise vector */
Vec norm_vec(Vec *v1)
{
	Vec v2;
	float len;

	len = sqrt(pow(v1->x, 2) + pow(v1->y, 2) + pow(v1->z, 2));

	v2.x = v1->x / len;
	v2.y = v1->y / len;
	v2.z = v1->z / len;

	return v2;
}

/*____________________________________________________________________________*/
/* scale vector: scale length of vector with given factor */
Vec scale_vec(Vec *v1, float factor)
{
	Vec v2;

    v2.x = v1->x * factor;
    v2.y = v1->y * factor;
    v2.z = v1->z * factor;

	return v2;
}

/*____________________________________________________________________________*/
/* vector dot product v1 * v2 */
float vec_dot_pro(Vec *v1, Vec *v2)
{
	return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
}

/*____________________________________________________________________________*/
/* vector cross product v1 x v2 */
Vec vec_cro_pro(Vec *v1, Vec *v2)
{
	Vec v3;

	v3.x = v1->y * v2->z - v1->z * v2->y;
	v3.y = v1->z * v2->x - v1->x * v2->z;
	v3.z = v1->x * v2->y - v1->y * v2->x;

	return v3;
}

/*____________________________________________________________________________*/
/* difference vector */
Vec diff_vec(Vec *v1, Vec *v2)
{
	Vec v3;

    v3.x = v1->x - v2->x;
    v3.y = v1->y - v2->y;
    v3.z = v1->z - v2->z;

	return v3;
}

/*____________________________________________________________________________*/
/* sum vector */
Vec sum_vec(Vec *v1, Vec *v2)
{
	Vec v3;

	v3.x = v1->x + v2->x;
	v3.y = v1->y + v2->y;
	v3.z = v1->z + v2->z;

	return v3;
}

/*____________________________________________________________________________*/
/* copy vector */
void copy_vec(Vec *v1, Vec *v2)
{
    v1->x = v2->x;
    v1->y = v2->y;
    v1->z = v2->z;
}

/*____________________________________________________________________________*/
/* print vector */
void print_vec(Vec *v)
{
	fprintf(stderr, "x %f, y %f, z %f\n", v->x, v->y, v->z);
}

/*____________________________________________________________________________*/
/* atom to vector */
void atom_to_vec(Vec *v, Atom *a)
{
    v->x = a->x;
    v->y = a->y;
    v->z = a->z;
}

/*____________________________________________________________________________*/
/* vector to atom */
 void vec_to_atom(Atom *a, Vec *v)
{
    a->x = v->x;
    a->y = v->y;
    a->z = v->z;
}

/*____________________________________________________________________________*/
/* vector angle */
/* cosinus of vector angle from dot product */
/* cos angle = v1 * v2 / ||v1|| * ||v2|| */
/* angle = acos (cos angle) */
float vec_ang(Vec *v1, Vec *v2)
{
    return (AN * acos(vec_dot_pro(v1, v2) / (vec_len(v1) * vec_len(v2))));
}

/*____________________________________________________________________________*/
/* average vector : average vector derived from an array of vectors */
Vec average_vec(Vec *v1, int nVec)
{
	unsigned int i;
	Vec v2;
	v2.x = 0.; v2.y = 0.; v2.z = 0.;

	for (i = 0; i < nVec; ++ i)
	{
		v2.x += v1[i].x;
		v2.y += v1[i].y;
		v2.z += v1[i].z;
	}

	v2.x /= nVec;
	v2.y /= nVec;
	v2.z /= nVec;

	return v2;
}

/*____________________________________________________________________________*/
/* calculate the dihedral torsion angle between 4 atoms (Calpha) */
 float calc_diheder_atom(Str *str, int atom1, int atom2, int atom3, int atom4)
{
	Vec d12, d23, d34; /* atom difference vectors */
	Vec p123, p234; /* plane normal vectors */
	float l_p123, l_p234; /* length of plane normal vectors */
	float angle_pp; /* angle between two plane normals */
	float sign;

	/* calculate atom difference vectors  */
	/* atoms 1,2 */
	d12.x = str->atom[atom1].x - str->atom[atom2].x;
	d12.y = str->atom[atom1].y - str->atom[atom2].y;
	d12.z = str->atom[atom1].z - str->atom[atom2].z;
	/* atoms 2,3 */
	d23.x = str->atom[atom3].x - str->atom[atom2].x;
	d23.y = str->atom[atom3].y - str->atom[atom2].y;
	d23.z = str->atom[atom3].z - str->atom[atom2].z;
	/* atoms 3,4 */
	d34.x = str->atom[atom3].x - str->atom[atom4].x;
	d34.y = str->atom[atom3].y - str->atom[atom4].y;
	d34.z = str->atom[atom3].z - str->atom[atom4].z;

	/* plane normals (cross prodduct): planes spanned by two difference vectors */
	/* plane 1,2,3 */
	p123.x = d12.y * d23.z - d12.z * d23.y;
	p123.y = d12.z * d23.x - d12.x * d23.z;
	p123.z = d12.x * d23.y - d12.y * d23.x;
	/* plane 2,3,4 */
	p234.x = d23.z * d34.y - d23.y * d34.z;
	p234.y = d23.x * d34.z - d23.z * d34.x;
	p234.z = d23.y * d34.x - d23.x * d34.y;

	/* length of plane normals */
	l_p123 = sqrt(pow(p123.x, 2) + pow(p123.y, 2) + pow(p123.z, 2));
	l_p234 = sqrt(pow(p234.x, 2) + pow(p234.y, 2) + pow(p234.z, 2));
	if (l_p123 == 0)
		l_p123 = 1e-6;
	if (l_p234 == 0)
		l_p234 = 1e-6;

	/* cosinus of angle between normals (dot product) */
	angle_pp = (p123.x * p234.x + p123.y * p234.y + p123.z * p234.z) / (l_p123 * l_p234);

	/* sign of angle -> rotation sense */
	sign =	d23.x * (p123.z * p234.y - p123.y * p234.z) + \
			d23.y * (p123.x * p234.z - p123.z * p234.x) + \
			d23.z * (p123.y * p234.x - p123.x * p234.y);

	/* account for rotation sense and convert from radian to angle */
	if (sign > 0)
		return (AN * (PI - acos(angle_pp)));
    else
        return (-1 * AN * (PI + (-1 * acos(angle_pp))));
}

/*____________________________________________________________________________*/
/* calculate the dihedral torsion angle between 3 vectors */
float calc_diheder_vector(Vec *d12, Vec *d23, Vec *d34)
{
	/*Vec d12, d23, d34: atom difference vectors */
	Vec p123, p234; /* plane normal vectors */
	float l_p123, l_p234; /* length of plane normal vectors */
	float angle_pp; /* angle between two plane normals */
	float sign;

	/* plane normals : planes spanned by two difference vectors */
	/* plane 1,2,3 */
	p123.x = d12->y * d23->z - d12->z * d23->y;
	p123.y = d12->z * d23->x - d12->x * d23->z;
	p123.z = d12->x * d23->y - d12->y * d23->x;
	/* plane 2,3,4 */
	p234.x = d23->z * d34->y - d23->y * d34->z;
	p234.y = d23->x * d34->z - d23->z * d34->x;
	p234.z = d23->y * d34->x - d23->x * d34->y;

	/* length of plane normals */
	l_p123 = sqrt(pow(p123.x, 2) + pow(p123.y, 2) + pow(p123.z, 2));
	l_p234 = sqrt(pow(p234.x, 2) + pow(p234.y, 2) + pow(p234.z, 2));
	/* avoid failure if vector length too small (at 90 deg) */
	if (l_p123 == 0)
		l_p123 = 1e-6;
	if (l_p234 == 0)
		l_p234 = 1e-6;

	/* normalised angle between normals */
	angle_pp = (p123.x * p234.x + p123.y * p234.y + p123.z * p234.z) / (l_p123 * l_p234);

	/* sign of angle -> rotation sense */
	sign =	d23->x * (p123.z * p234.y - p123.y * p234.z) + \
			d23->y * (p123.x * p234.z - p123.z * p234.x) + \
			d23->z * (p123.y * p234.x - p123.x * p234.y);

	/* account for rotation sense and convert from radian to angle */
	if (sign > 0)
		return (float)(AN * (PI - acos(angle_pp)));
    else
        return (float)(-1 * AN * (PI + (-1 * acos(angle_pp))));
}

/*___________________________________________________________________________*/
/* calculate the PHI/PSI angles per residue and define sec.str. elements */
/* str->ss[res][0] is element number, str->ss[res][1] is sec.str. type */
/*
	float (*phi)[6] : array of backbone angles PHI;
						0: angle, 2-5:4 atom numbers, 6: sec.str.
	float (*psi)[6] : array of backbone angles PSI;
						0: angle, 2-5: 4 atom numbers, 6: sec.str.
	int (*ss)[2] : array of sec.str. elements;
						0: element number, 1: type of sec.str. per residue
*/
void phi_psi(Str *str, int res)
{
	/* PHI */
	if (res > 0) {
		str->phi[res][0] = calc_diheder_atom(str, (int)str->phi[res][1], \
													(int)str->phi[res][2], \
													(int)str->phi[res][3], \
													(int)str->phi[res][4]);
	}

	/* PSI */
	if (res < (str->nResidue - 1)) {
		str->psi[res][0] = calc_diheder_atom(str, (int)str->psi[res][1], \
													(int)str->psi[res][2], \
													(int)str->psi[res][3], \
													(int)str->psi[res][4]);
	}

	/* assign sec.str. from Phi/Psi angles */
	/* beware of logic! : a central range has '&&', a boundary range has '||' */
    /* sec.str. definitions taken from: Kleywegt and Jones (1996) Current Biology 4, 1395-1400 (95% Ramachandran lines) */
	/* sheet, type 0 */
	/* We limit the assignment to residues between (including) 1 and res-1,
		because PHI is not defined for residue 0 and PSI not for res-1. */
	if ((res > 0) && (res < (str->nResidue - 1))) {
		if ((str->phi[res][0] >= -180 && str->phi[res][0] <= -60) && \
			(str->psi[res][0] >=  60  && str->psi[res][0] <= 180)) {
			str->ss[res][1] = 0; /* sec.str. type of this element */
		/* right-handed helix, type 1 */
		} else if ((str->phi[res][0] >= -120 && str->phi[res][0] <= -30) && \
				 (str->psi[res][0] >=  -75 && str->psi[res][0] <=  30)) {
			str->ss[res][1] = 1;
		/* left-handed helix, type 2 */
		} else if ((str->phi[res][0] >=  30 && str->phi[res][0] <= 80) && \
				 (str->psi[res][0] >= -25 && str->psi[res][0] <= 70)) {
			str->ss[res][1] = 2;
		/* anything else (= coil), type 3 */
		} else {
			str->ss[res][1] = 3;
		}
	}

	/* sec.str. elements : these are not! the sec.str. type numbers but */
	/*   numbered elements of contiguous sec.str. types */
	/* We limit the assignment to residues between (including) 1 and res-1,
		because secondary structures are not defined for the terminal residues. */
	if (res == 1) {
		str->ss[res][0] = 0; /* start with element number 0 at residue 1 */
		str->ss[res][1] = 3; /* initalise as ss type 'other' */
	} else if ((res > 1) && (res < (str->nResidue - 1))) {
		if (str->ss[res][1] == str->ss[res - 1][1]) {
			/* continue element if sec.str. types equal */
			str->ss[res][0] = str->ss[res - 1][0];
		} else {
			/* or generate new element */
			str->ss[res][0] = str->ss[res - 1][0] + 1; 
		}
	} else {
		str->ss[res][0] = -1; /* NA */
		str->ss[res][1] = -1; /* NA */
	}

	if (res == 0) {
		fprintf (stderr, "%s:%d: res %d, psi %f, ss_element %d, type %d\n",
			__FILE__, __LINE__, res, str->psi[res][0], str->ss[res][0], str->ss[res][1]);
	} else if (res >= (str->nResidue - 1)) {
		fprintf (stderr, "%s:%d: res %d, phi %f, ss_element %d, type %d\n",
			__FILE__, __LINE__, res, str->phi[res][0], str->ss[res][0], str->ss[res][1]);
	} else {
		fprintf (stderr, "%s:%d: res %d, phi %f, psi %f, ss_element %d, type %d\n",
			__FILE__, __LINE__,
			res, str->phi[res][0], str->psi[res][0], str->ss[res][0], str->ss[res][1]);

	}
}

/*____________________________________________________________________________*/
/* define sec.str. segments based on the sec.str. elements defined above (str->ss) */
/* Two versions coded:
	Version 1 : short segments are included into preceeding long segments
	Version 2 : short segments are excluded from segments */
/*
	int (*seg)[4] : array of sec.str. segment lengths; 
						0: start atom, 1: length, 2: sec.str. type, 3: contact
	float (*phit)[2] : array of segment angles PHIt: 0: angle, 1: topology
	float (*psit)[2] :  array of segment angles PSIt: 0: angle, 1: topology
	int nseg : number of sec.str. segments
*	ss types : 0: sheet, 1: rh helix, 2: lh helix, 3: other
*/
void ss_segments(Str *str)
{
	int i;
	int allocated = 64;

	str->nseg = 0;
    str->seg = safe_malloc(allocated * sizeof(int [4]));
    str->phit = safe_malloc(allocated * sizeof(float [2]));
    str->psit = safe_malloc(allocated * sizeof(float [2]));

	/* 'seg' is function of segment number */
	/* first atom (number 0) */
	str->seg[str->nseg][0] = 0; /* start atom of segment */
	str->seg[str->nseg][1] = 1; /* segment length */
	str->seg[str->nseg][2] = -1; /* ss type unassigned */
	str->seg[str->nseg][3] = -1; /* contact unassigned */

	for (i = 1; i < str->sequence.length; ++ i) {
		if (str->ss[i][0] == str->ss[i-1][0]) { /* if same sec.str. element */
			++ str->seg[str->nseg][1]; /* extend segment */
			str->seg[str->nseg][2] = str->ss[i][1];
		}
		/* segment assignment:
			if too short : don't assign segment; helices at least 5 long,
				sheet at least 4 long, 'others' don't count as segment */
        else if (((str->seg[str->nseg][1] < 5) &&
				 ((str->ss[str->seg[str->nseg][0]][1] == 1) ||
				  (str->ss[str->seg[str->nseg][0]][1] == 2))) ||
                 ((str->seg[str->nseg][1] < 4) &&
				  (str->ss[str->seg[str->nseg][0]][1] == 0)) ||
				  (str->ss[str->seg[str->nseg][0]][1] == -1) ||
				  (str->ss[str->seg[str->nseg][0]][1] == 3))
		{
			str->seg[str->nseg][0] = i; /* segment starts here */
			str->seg[str->nseg][1] = 1; /* ... and is 1 long */
			str->seg[str->nseg][2] = -1; /* ss type unassigned */
			str->seg[str->nseg][3] = -1; /* contact unassigned */
		} else { /* else new segment */
			++ str->nseg;

			/* memory allocation */
			if (str->nseg == allocated)
			{
				allocated += 64;
				str->seg = safe_realloc(str->seg, allocated * sizeof(int [4]));
				str->phit = safe_realloc(str->phit, allocated * sizeof(float [2]));
				str->psit = safe_realloc(str->psit, allocated * sizeof(float [2]));
			}

			str->seg[str->nseg][0] = i; /* start atom of new segment */
			str->seg[str->nseg][1] = 1; /* initialise segment length */
			str->seg[str->nseg][2] = -1; /* ss type unassigned */
			str->seg[str->nseg][3] = -1; /* contact unassigned */
		}
		str->atom[i].seg = str->nseg; /* assign segment number to this atom */

		fprintf (stderr, "%s:%d: res = %d, seg = %d, seg_start = %d, seg_length = %d\n",
					__FILE__, __LINE__,
					i, str->nseg, str->seg[str->nseg][0], str->seg[str->nseg][1]);
	}
	/* check last segment */
    if (((str->seg[str->nseg][1] < 5) &&
		((str->ss[str->seg[str->nseg][0]][1] == 1) ||
         (str->ss[str->seg[str->nseg][0]][1] == 2))) ||
        ((str->seg[str->nseg][1] < 4) &&
		 (str->ss[str->seg[str->nseg][0]][1] == 0)) ||
		 (str->ss[str->seg[str->nseg][0]][1] == 3))
	{
		-- str->nseg;
	}
}

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
	unsigned int i;
	Vec *subAxis;
	/* axis normals */
	Vec *axisNormal;
	Vec tmpa1, tmpa2, tmpa3; /* temporary postition vectors */
	Vec dtmpa1, dtmpa2; /* temporary difference vectors */
	Vec stmpa; /* temporary sum vectors */
	/* subaxes */
	float *dAxis; /* distance between Calpha projections onto helix axis */
	float *rAxis; /* radius of helix axis */
	float avrAxis = 0; /* average rAxis */
	Vec tmpb2, tmpb3;
	Vec dtmpb3, dtmpb4;
	float t1, t2, t3; /* terms in helix radius formula */
	/* axis points */
	Vec tmpc1, tmpc2, tmpc3, tmpc4, tmpc5;

	/*_______________________________________________________________________*/
	/* construct all axis normals, which are angle bisection vectors */
	/* bisection vectors are averaged angle vectors */
	axisNormal = safe_malloc((str->seg[seg][1] - 1) * sizeof(Vec));

	for (i = 1; i < str->seg[seg][1] - 1; ++ i)
	{
		/* vectors to the three angle-constituting atoms */
		atom_to_vec(&tmpa1, &(str->atom[str->seg[seg][0] + i - 1]));
		atom_to_vec(&tmpa2, &(str->atom[str->seg[seg][0] + i]));
		atom_to_vec(&tmpa3, &(str->atom[str->seg[seg][0] + i + 1]));
		/* two vectors constituting the angle arms */
		dtmpa1 = diff_vec(&tmpa1, &tmpa2);
		dtmpa2 = diff_vec(&tmpa3, &tmpa2);
		/* angle-disecting vector (going through helix axis) */
		stmpa = sum_vec(&dtmpa1, &dtmpa2);
		/* normalise to length 1A */
		/* this vector is called V1 or V2 in Kahn's paper */
		axisNormal[i] = norm_vec(&stmpa);
	}

	/*_______________________________________________________________________*/
	/* To calculate the central axis, subaxes are constructed from 
		vector cross products of axis normals. */
	/* number of subaxes */
	subAxis = safe_malloc((str->seg[seg][1] - 2) * sizeof(Vec));
	dAxis = safe_malloc((str->seg[seg][1] - 2) * sizeof(Vec));
	rAxis = safe_malloc((str->seg[seg][1] - 2) * sizeof(Vec));

	/* vector cross product: this is called 'V1 x V2' or 'H' in Kahn's paper */
	for (i = 1; i < str->seg[seg][1] - 2; ++ i)
	{
		/* helix axis direction 'H' */
		subAxis[i] = vec_cro_pro(&(axisNormal[i]), &(axisNormal[i+1]));

		/* Calpha distance projected onto helix axis */
		/* the two central alpha atoms of V1 and V2 */
		/* their position vectors */
		atom_to_vec(&tmpb2, &(str->atom[str->seg[seg][0] + i]));
		atom_to_vec(&tmpb3, &(str->atom[str->seg[seg][0] + i + 1]));
		/* the difference vector between them */
		dtmpb3 = diff_vec(&tmpb3, &tmpb2);
		/* their distance when projected onto the helix axis */
		/* d = (P2 - P1) H */
		dAxis[i] = vec_dot_pro(&dtmpb3, &(subAxis[i]));

		/* helix radius */
		/* r = |d*H|^2 - |(P2 - P1)|^2 / (2 |(P1 - P2) V2|) */
		/* first term */
		t1 = pow(dAxis[i], 2);
		/* second term */
		t2 = vec_len_sq(&dtmpb3);
		/* third term */
		dtmpb4 = diff_vec(&tmpb2, &tmpb3);
		t3 = vec_dot_pro(&dtmpb4, &(axisNormal[i+1]));
		/* radius */
		rAxis[i] = (t2 - t1) / (2 * t3);

		avrAxis += rAxis[i];
	}

	/* average axis radius */
	avrAxis /= (str->seg[seg][1] - 3); 
	/*fprintf(stderr, "radius %f\n", avrAxis);*/

	/*_______________________________________________________________________*/
	/* N-terminus */
	/* position vector to atom (= 1) of first axisNormal */
	atom_to_vec(&tmpc1, &str->atom[str->seg[seg][0] + 1]);
	tmpc2 = scale_vec(&axisNormal[1], rAxis[1]);
	str->axispoint[seg][0] = sum_vec(&tmpc1, &tmpc2);

	/* C-terminus */
	/* position vector to atom of last axisNormal */
	/* note that atom numbers of C-termini are 'seg[seg][0] + seg[seg][1]' */
	atom_to_vec(&tmpc3, &str->atom[str->seg[seg][0] + str->seg[seg][1] - 3]);
	tmpc4 = scale_vec(&axisNormal[str->seg[seg][1] - 3], rAxis[str->seg[seg][1] - 3]);
	str->axispoint[seg][1] = sum_vec(&tmpc3, &tmpc4);

	/* midpoint */
	/* half-distance between N-terminus and C-terminus */
	tmpc5 = diff_vec(&str->axispoint[seg][1], &str->axispoint[seg][0]);
	tmpc5 = scale_vec(&tmpc5, 0.5);
	str->axispoint[seg][2] = sum_vec(&str->axispoint[seg][0], &tmpc5);

	fprintf(stderr, "%s:%d: seg %d, Nt %f, Ct %f, mi %f\n",
				__FILE__, __LINE__, seg,
				str->axispoint[seg][0].x, str->axispoint[seg][1].x,
				str->axispoint[seg][2].x);

	/*_______________________________________________________________________*/
    free(axisNormal);
    free(subAxis);
    free(dAxis);
    free(rAxis);
}

/*____________________________________________________________________________*/
/* strand axis */
/* The strand axis is the average vector of (all) pairs of connection
	vectors (subAxis) between the midpoints of the difference vectors
	C1,C2 and C2,C3. */
/* defines axis points: 0: N-terminus, 1: C-terminus, 2: midpoint */
 void strand_axis(Str *str, int seg)
{
	unsigned int i;
	Vec *subAxis;
	Vec tmp1, tmp2, tmp3, tmp4;
	Vec dtmp1, dtmp2;

	/* subaxes */
	subAxis = safe_malloc((str->seg[seg][1] - 1) * sizeof(Vec));

	for (i = 1; i < str->seg[seg][1] - 1; ++ i) {
		/* vectors to angle forming atoms */
		atom_to_vec(&tmp1, &(str->atom[str->seg[seg][0] + i - 1]));
		atom_to_vec(&tmp2, &(str->atom[str->seg[seg][0] + i]));
		atom_to_vec(&tmp3, &(str->atom[str->seg[seg][0] + i + 1]));
		/* half-length connecting vectors between atoms are on axis */
		dtmp1 = diff_vec(&tmp2, &tmp1);
		dtmp1 = scale_vec(&dtmp1, 0.5);
		dtmp2 = diff_vec(&tmp2, &tmp3);
		dtmp2 = scale_vec(&dtmp2, 0.5);
		/* strand axis: connecting vector is in direction of axis */
		subAxis[i] = diff_vec(&dtmp2, &dtmp1);
		/* axis point N-terminus: position vector to first atom of segment */
		if (i == 1)
			str->axispoint[seg][0] = sum_vec(&tmp1, &dtmp1);
		/* axis point C-terminus: position vector to last atom of segment */
		if (i == str->seg[seg][1] - 2)
			str->axispoint[seg][1] = sum_vec(&tmp1, &dtmp1);
	}
	/* axis point midpoint: half-distance between N-terminus and C-terminus */
	tmp4 = diff_vec(&str->axispoint[seg][1], &str->axispoint[seg][0]);
	tmp4 = scale_vec(&tmp4, 0.5);
	str->axispoint[seg][2] = sum_vec(&str->axispoint[seg][0], &tmp4);

	fprintf(stderr, "%s:%d: seg %d, Nt %f, Ct %f, mi %f\n",
				__FILE__, __LINE__, seg,
				str->axispoint[seg][0].x, str->axispoint[seg][1].x,
				str->axispoint[seg][2].x);

	free(subAxis);
}

/*____________________________________________________________________________*/
/* secondary structures' segment axis and pairwise angles */
/* Define the three points that are needed per residue for diheder calculation:
	0: axis N-terminus, 1: axis C-terminus, 2: axis midpoint. */
void segment_angle(Str *str)
{
	unsigned int seg;
	/*int atom1, atom2, atom3, atom4;*/
	Vec d12, d23, d34;

	/* define segment axes and points on axes */
	for (seg = 0; seg < str->nseg; ++ seg) {
		/* strand axis */
		if (str->ss[str->seg[seg][0]][1] == 0) {
			strand_axis(str, seg);
		/* helix (lh, rh) axis */
		} else if (str->ss[str->seg[seg][0]][1] == 1 ||
			str->ss[str->seg[seg][0]][1] == 2) {
			helix_axis(str, seg);
		} else {
			fprintf(stderr, "%s%d: seg %d, ss %d: This should not happen!\n",
					__FILE__, __LINE__, seg, str->ss[str->seg[seg][0]][1]);
			exit(1);
		}

		if (seg == 0) {
			str->phit[seg][0] = 0;
		} else {
			d12 = diff_vec(&(str->axispoint[seg - 1][0]), &(str->axispoint[seg - 1][2]));
			d23 = diff_vec(&(str->axispoint[seg][2]), &(str->axispoint[seg - 1][2]));
			d34 = diff_vec(&(str->axispoint[seg][2]), &(str->axispoint[seg][1]));

			fprintf(stderr, "%s:%d: axes distance vectors; d12 %f, d23 %f, d34 %f\n",
							__FILE__, __LINE__, d12.x, d23.x, d34.x);

			str->phit[seg][0] = abs(calc_diheder_vector(&d12, &d23, &d34));
		}
	}
}

