/* HEADER.h
 *
 *Copyright 1995-96 Roger P. Woods, M.D.
 *Modified 1/6/96
 *
 *This file defines the image header structure to be used by AIR.
 *
 *The current default file format defined here for AIR is intended to be compatible
 *with the ANALYZE format created by Richard A. Robb's group at the Mayo Clinic. I only have
 *direct access to ANALYZE version 2.0 which runs on a Sun workstation, and  I therefore cannot
 *guarantee compatibility of the struct below with newer versions of ANALYZE or other platforms.
 *Consequently, the ANALYZE struct names for the various variables used by AIR are provided at the 
 *end of this file to facilitate reconciliation of any problems.
 *
 *You should be able to adapt AIR to read and write image files in the format of your choice by
 *modifying this file and possibly the functions read_header(), write_header(),
 *read_image(), quickread(), save(), and map_value().
 *
 *Note that AIR only uses part of the information in the header and does not provide a 
 *general mechanism for propagating variable information that is not intrinsic to AIR.
 *
 *Note: ANALYZE is a registered trademark of the Mayo Foundation. A description of ANALYZE can be
 *found in: Robb, RA, Hanson, DP. A software system for interactive and quantitative visualization
 *of multidimensional biomedical images. Australasian Physical & Engineering Sciences in Medicine
 *1991;14:9-30.
 */
 

struct hdr{
	int		sizeof_hdr;	/* For ANALYZE compatibility only */
	char		pad1[28];
	int		extents;	/* For ANALYZE compatibility only */
	char		pad2[2];
	char		regular;	/* For ANALYZE compatibility only */
	char		pad3;
	short int	dims;		/* For ANALYZE compatibility only */
	short int	x_dim;		/* AIR */
	short int	y_dim;		/* AIR */
	short int	z_dim;		/* AIR */
	short int	t_dim;		/* For ANALYZE compatibility only */
	char		pad4[20];
	short int	datatype;	/* For ANALYZE compatibility only */
	short int	bits;		/* AIR */
	char		pad5[6];
	float		x_size;		/* AIR */
	float		y_size;		/* AIR */
	float		z_size;		/* AIR */
	char		pad6[48];
	int		glmax;		/* AIR */
	int		glmin;		/* AIR */
	char		descrip[80];	/* AIR (non-essential) */
	char		pad7[120];
};


/* Translation table:

	AIR header struct		ANALYZE header struct

 	hdr.sizeof_hdr			hdr.hk.sizeof_hdr
 	hdr.extents			hdr.hk.extents
 	hdr.regular			hdr.hk.regular
	hdr.dims			hdr.dime.dim[0]
	hdr.x_dim			hdr.dime.dim[1]
	hdr.y_dim			hdr.dime.dim[2]
	hdr.z_dim			hdr.dime.dim[3]
	hdr.t_dim			hdr.dime.dim[4]
	hdr.datatype			hdr.dime.datatype
	hdr.bits			hdr.dime.bitpix
	hdr.x_size			hdr.dime.pixdim[1]
	hdr.y_size			hdr.dime.pixdim[2]
	hdr.z_size			hdr.dime.pixdim[3]
	hdr.glmax			hdr.dime.glmax
	hdr.glmin			hdr.dime.glmin
	hdr.descrip			hdr.hist.descrip
*/
	
