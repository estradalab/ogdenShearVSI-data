/* AIR.h --header file for Automated Image Registration Subroutines AIR 3.00 */
/* Copyright 1995-96 Roger P. Woods, M.D. */
/* Modified 10/26/96 */

/* This header file utilizes a 4 x 4 matrix rather than the older 4 x 3 format*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* The following information is site specific */

#define PLATF 1			/* Used to properly define machine dependent constants*/
				/* Currently supported:					*/
				/* 1 IEEE big-endian (Motorola 68000's, Sun 3, SPARC)	*/
				/* 2 IEEE little-endian (doesn't work for DEC alpha) */
				/* 3 Platform independent solution for C compilers with float.h header file */
				/* (DEC alpha includes float.h, so use 3 for alphas) */

#define VERBOSITY 0             /*0 turns off nonessential screen printing*/

#define PIX_SIZE_ERR .0001	/* Voxel sizes that differ by less than this value are assumed identical */

#define PIXEL_MAX_SIZE 10	/* warning issued if voxel size exceeds this*/
#define PIXEL_MIN_SIZE .1	/* warning issued if voxel size less than this*/

#define EDITED 5.0		/* warning issued if MRI contains less than this percentage of zero pixel values*/
				/* suggesting that it has not been edited--value of 0.0 inactivates this feature*/
				/* This feature is supported only in the AIR 1.0 programs alignpettomri */
				/* and alignmritopet. The new program alignlinear does not check to see that */
				/* files have been edited.							*/

#define OLD_ANALYZE 1		/* if 1, warning will be issued whenever*/
				/* a file is loaded that has voxels sizes of */
				/* 1.000 x 1.000 x1.000 suggesting that*/
				/* file was created by a version of ANALYZE	*/
				/* that sometimes fails to preserve voxel sizes		*/
				/* (particularly when a file has been edited*/
				/* set to 0 (zero) to inactivate*/



#define COMPRESS 1		/* if nonzero, image files (but not headers) can be decompressed*/
				/* and recompressed using UNIX compress and uncompress commands */
				/* This should be set to zero if your operating system does not */
				/* support the uncompress and compress commands.                */
				/* Note that only files that are compressed when initially      */
				/* encountered will be recompressed and that not all programs   */
				/* support this feature.                                        */
				/* If your system does not support this feature, you may crash  */
				/* instead of exiting gracefully when you ask for a nonexistant */
				/* file to be loaded. Set COMPRESS to 0 to prevent this.        */
				/* Alternatively, you can modify the subroutines                */
				/* src/decompress.c and src/recompress.c to make system calls   */
				/* appropriate for your systems architecture.                   */
				/* Likewise, if you prefer other compression programs, these    */
				/* subroutines can be modified to call them instead.            */

#define CLOCK 0			/* if nonzero, verbose mode in alignlinear will print elapsed
				 * time with each iteration
				 *
				 * this is a nonessential function and causes trouble with some
				 * compilers. If you get a clock() related error when you compile,
				 * set this variable to zero
				 */

#define OUTBITS 8		/* Options are 8 and 16 */
				/* Controls internal representation of data */
				/* Controls bits/pixel of output data*/
				/* Data with a different number of bits/pixel can be input*/
				/* but will be converted to this number of bits/pixel*/
				/* using bit shifting to increase the number*/
				/* and using bit shifting, possibly combined with rescaling*/
				/* to decrease the number of bits/pixel.*/
				/* (Whether or not to rescale is dictated by the header*/
				/* global maximum					*/


#define REP16 3			/* This variable is only relevant when OUTBITS==16 	*/
				/* However, its value must always be defined		*/
				/* If 1, 16 bit data will be written to disk as		*/
				/* 	unsigned short ints (NOT ANALYZE compatible) 	*/
				/*	header minimum will be 0, maximum 65535		*/
				/* If 2, 16 bit data will be written to disk as		*/
				/* 	short ints w/o using negative values		*/
				/*	this effectively reduces storage to 15 bits	*/
				/* 	and reduction is done via bit shifting		*/
				/*	header minimum will be 0, maximum 32767		*/
				/* If 3, 16 bit data will be written to disk as		*/
				/* 	short ints using negative values		*/
				/*	an actual value of zero will be mapped to	*/
				/* 	-32768 in short int representation		*/
				/*	header minimum will be -32767 			*/
				/*	header maximum will be 32767 			*/

#if(OUTBITS==8)

#define THRESHOLD1 55		/* These are the 8 bit thresholds (see below) */
#define THRESHOLD2 55

#elif(OUTBITS==16)

#define THRESHOLD1 7000		/* These are the 16 bit thresholds (see below) */
#define THRESHOLD2 7000

#endif
				/* THRESHOLD1 and THRESHOLD2 are the default thresholds
				 * that AIR will offer when you register images.
				 *
				 * You probably want to use the same value for THRESHOLD1 and THRESHOLD2.
				 *
				 * You can always override the defaults, but it's nice to compile
				 * in some reasonable values if you can.
				 *
				 * Ideally, you should look at some typical data, and pick a default
				 * value that will generally separate brain from non-brain structures.
				 *
				 * In other words, pick a threshold such that values outside the structure
				 * of interest have pixel values below that threshold.
				 *
				 * An 8 bit value of 55 is good for PET data, but you'll probably 
				 * want something closer to 10 for 8 bit MRI data
				 *
				 * For 16 bit data, it's hard to guess what value will be reasonable
				 * since the full 16 bit range is often not utilized. A value around
				 * 7000 is reasonable for PET data which is often effectively 15 positive bits.
				 * For MRI data, as little as 12 bits is commonly used, in which case a
				 * reasonable default threshold is probably in the 100-200 range. On the other
				 * hand, if all 16 bits are effectively utilized, values in the 1000-2000 range
				 * may be more appropriate.
				 *
				 * Best bet--look at some data before choosing a default and caution users
				 * not to rely on the default values unless your site is very consistent
				 * about data formats and image intensities.
				 */


/* The remaining information should not generally require modification */

#if(OUTBITS==8)

typedef unsigned char my_pixels;
#define MAX_POSS_VALUE 255

#elif(OUTBITS==16)

typedef unsigned short int my_pixels;
#define MAX_POSS_VALUE 65535	/* This is true regardless of the value of REP16 above*/
				/* 16-bit data is always represented internally is*/
				/* unsigned short ints, even when it is written to*/
				/* output files in other formats*/

#endif


#define HDR_SUFF ".hdr"		/* This suffix will be used in trying to open an image header*/
#define IMG_SUFF ".img"		/* This suffix will be used in trying to open an image file*/

#define LOGFILE "\0"		/* Use "\0" if .img files don't have corresponding .log files*/
#define NORM_SUFF ".nrm"	/* Use "\0" if .log files containing ASCII weighting factors are not available*/


/* These are the internal structs*/
/* The external struct (for image headers as stored on disk)*/
/* is stored separately in HEADER.h */

struct key_info{
        int bits;
        int x_dim;
        int y_dim;
        int z_dim;
        double x_size;
        double y_size;
        double z_size;
};

struct air16{
        double  		e[4][4];
        char    		s_file[128];
        struct key_info		s;
        char    		r_file[128];
        struct key_info 	r;
        char    		comment[128];
	unsigned long int	s_hash;
	unsigned long int	r_hash;
	unsigned short		s_volume;	/* Not used in this version of AIR */
	unsigned short		r_volume;	/* Not used in this version of AIR */
	char			reserved[116];
};

struct oldair{					/*Allows AIR 2.0 to read AIR 1.0 .air files*/
        double                  e[4][3];
        char                    s_file[128];
        struct key_info         s;
        char                    r_file[128];
        struct key_info         r;
        char                    comment[128];
        unsigned long int       s_hash;
        unsigned long int       r_hash;
	unsigned short		s_volume;	/* Not used in this version of AIR */
	unsigned short		r_volume;	/* Not used in this version of AIR */
        char                    reserved[116];
};
