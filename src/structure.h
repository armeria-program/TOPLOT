/*===============================================================================
 structure.h : Routines for structure operations
 Copyright (C) 2004 Jens Kleinjung
 GNU GPL License applies
================================================================================*/

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "config.h"
#include "toplot.h"

#define PI (4.0 * atan(1.0))
#define AN (180 / PI)

/*____________________________________________________________________________*/
/* structures */

/* vector is array of floats */
typedef struct
{
   float x, y, z;
} Vec;

/* sequence */
typedef struct  
{
    char *name; /* sequence name */
    char *res; /* array of residues = sequence */
    int length; /* length of sequence */
	int N; /* N-terminal position */
	int C; /* C-terminal position */
} Seq;

/* atom : definition of PDB atom format, numbers indicate columns */
/* only the entries used in this program are activated so save resources */
typedef struct
{
	/* PDB data */
	char record_ne[8]; /* Record type; 1 -  6*/
	int atom_nr; /* Atom serial number;  7 - 11 */
	char atom_ne[8]; /* Atom name; 13 - 16 */
	char alt_loc[2]; /* Alternate location indicator; 17 */
	char res_ne[4]; /* Residue name; 18 - 20 */
	char chain_id[2]; /* Chain identifier; 22 */
	int res_nr; /* Residue sequence number; 23 - 26 */
	/*char icode[2];*/ /* Code for insertion of residues; 27 */
	float x; /* Orthogonal coordinates for X in Angstroms; 31 - 38 */
	float y; /* Orthogonal coordinates for Y in Angstroms; 39 - 46 */
	float z; /* Orthogonal coordinates for Z in Angstroms; 47 - 54 */
	/*float occupancy;*/ /* Occupancy; 55 - 60 */
	/*float temp_f;*/ /* Temperature factor; 61 - 66 */
	/*char seg_id[5];*/ /* Segment identifier; 73 - 76 */
	/*char element[3];*/ /* Element symbol; 77 - 78 */
	/*char charge[3];*/ /* Charge on the atom; 79 - 80 */
	char descrip[32]; /* everything before coordinates */
	float tx; /* transformed x : tx = U * x + t */
	float ty; /* transformed y : ty = U * y + t */
	float tz; /* transformed z : tz = U * z + t */
	/* derived data */
	int seg; /* number of sec. str. segment this atom belongs to */
} Atom;

/* structure */
typedef struct
{
	Atom *atom; /* array of selected (CA) atoms constituting structure */
	int nAtom; /* number of selected (CA) atoms */
	int nResidue; /* number of selected (CA) atoms */
	int nChain; /* number of chains */

	Seq sequence; /* sequence of structure */
	float (*phi)[6]; /* array of backbone angles PHI:
						0: angle, 2-5: 4 atom numbers, 6: sec.str. */
	float (*psi)[6]; /* array of backbone angles PSI:
						0: angle, 2-5: 4 atom numbers, 6: sec.str. */
	int (*ss)[2]; /* array of sec.str. elements: 0:
						0: element number, 1: type of sec.str. per residue */
	int (*seg)[4]; /* array of sec.str. segment lengths; 
						0: start atom, 1: length, 2: sec.str. type, 3: contact */
	float (*phit)[2]; /* array of segment angles PHIt: 0: angle, 1: topology */
	float (*psit)[2]; /* array of segment angles PSIt: 0: angle, 1: topology */
	int nseg; /* number of sec.str. segments */
	Vec *axis; /* vector for segment axis definition */
	Vec (*axispoint)[3]; /* vector for points on axis:
						0: N-terminus, 1: C-terminus, 2: midpoint */
} Str;

/* prototypes */
float vec_len(Vec *v);
float vec_len_sq(Vec *v);
Vec norm_vec(Vec *v1);
Vec scale_vec(Vec *v1, float factor);
Vec vec_cro_pro(Vec *v1, Vec *v2);
float vec_dot_pro(Vec *v1, Vec *v2);
Vec sum_vec(Vec *v1, Vec *v2);
void copy_vec(Vec *v1, Vec *v2);
void print_vec(Vec *v);
void atom_to_vec(Vec *v, Atom *a);
Vec diff_vec(Vec *v1, Vec *v2);
float vec_ang(Vec *v1, Vec *v2);
Vec average_vec(Vec *v1, int nVec);
void phi_psi(Str *str, int res);
void ss_segments(Str *str);
void segment_angle(Str *str);

#endif
