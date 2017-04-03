/*                                                                
**  Copyright (C) 2004-2008,2016,2017  Smithsonian Astrophysical Observatory 
*/

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */




#include <dslib.h>
#include <ascdm.h>
#include <stdlib.h>
#include <math.h>
#include <histlib.h>

#include <cxcregion.h>
#include <dsnan.h>

#define FLOOR(x)  ((double)(x))
#define CEIL(x)   ((double)(x))

regRegion *GlobalMaskRegion;


/* Okay, global variables are a bad idea, excpet when doing recursion and
 * we are passing the same data to each level ... more data gets pushed on the
 * stack and causes crashes.  Also using global variables can make it
 * run much faster so we'll bite the bullet here. */
void *GlobalData;               /* i: data array */
float *GlobalDErr;              /* i: error array */
float *GlobalOutData;           /* o: output array */
float *GlobalOutArea;           /* o: output area array */
float *GlobalOutSNR;            /* o: output SNR */
unsigned long *GlobalMask;      /* o: output mask */
long GlobalXLen;                /* i: length of x-axis (full img) */
long GlobalYLen;                /* i: length of y-axis (full img) */
long *GlobalLAxes;            /* X,Y lens togeether */
dmDataType GlobalDataType;
dmDescriptor *GlobalXdesc = NULL;
dmDescriptor *GlobalYdesc = NULL;
short *GlobalPixMask = NULL;


float GlobalSNRThresh;          /* i: SNR threshold */
enum { ZERO_ABOVE = 0, ONE_ABOVE, TWO_ABOVE, THREE_ABOVE, ALL_ABOVE } GlobalSplitCriteria;
double GlobalX0;                // =  4274.5; // 4001.9; //4274.0;
double GlobalY0;                // =  3956.0; // 3850.5; // 3954.0;
double GlobalInnerRadius;       // = 5.0;
double GlobalOuterRadius;       // = 1000.0;
double GlobalStartAngle;        // = 0.0;
double GlobalStopAngle;         // = 360.0;
double GlobalMinRadius;         // = 0.5;
double GlobalMinAngle;          // = 1 ;
double GlobalEllipticity;
int (*GlobalShapeFunction)(regRegion *r, double a_m, double b_m, double a_l, double b_l) = NULL;
void (*GlobalLimitsFunction)( double a_m, double b_m, double a_l, double b_l, double p[4], double q[4]) = NULL;


typedef struct {
    char infile[DS_SZ_FNAME];
    char errimg[DS_SZ_FNAME];
    char outfile[DS_SZ_FNAME];
    char areafile[DS_SZ_FNAME];
    char maskfile[DS_SZ_FNAME];
    char snrfile[DS_SZ_FNAME];
    char shape[100];
    short method;
    short clobber;
} Parameters;


/* Using the dmtools/dmimgio routines removes lots of duplicate code that was
 * originally here.  Also allow us to keep track of NULL/NaN value pixels
 * more easily
 */
#include "dmimgio.h"

#ifdef DEMO_MODE
FILE *global_demo_file;
#endif


/* ------Prototypes ----------------------- */

int load_error_image(char *errimg);
int abin(void);

double get_snr(double a_min, double b_min, double a_len, double b_len,
               float *oval, long *area);
void abin_rec(double a_min, double b_min, double a_len, double b_len);
int convert_coords(dmDescriptor * xdesc, dmDescriptor * ydesc, double xx,
                   double yy, double *xat, double *yat);
int invert_coords(dmDescriptor * xdesc, dmDescriptor * ydesc, double xat,
                  double yat, double *ii, double *jj);
regRegion *make_region(regRegion * inreg, double a_min, double b_min,
                    double a_len, double b_len, long *xs, long *xl,
                    long *ys, long *yl);
void fill_region(double a_min, double b_min, double a_len, double b_len);


int write_single_output(dmBlock *inBlock, char *outfile, void *outdata, 
        short clobber, dmDataType dt, regRegion *outreg);
int write_outputs( dmBlock *inBlock, Parameters *pp);

dmBlock *load_infile( char *infile);
int allocate_memory(void);
int autoname( Parameters *pp );
int map_method( short method ) ;
void free_globals(void);
Parameters *load_parameters(void);
int aux_output(char *outfile, char *auxoutfile, char *suffix, short clobber );
int map_shapes( char *shape );


/*
 *  Note:  In the interface, if using polar grid then:
 * 
 *     a_min is the starting radius, a_len is the length (so from a_min to a_min+a_len).
 *     b_min is the starting angle, b_len is the length (so from b_min to b_min+b_len).
 * 
 * If using a rectangular grid, then
 *     a_min is the *middle* of the rectangle x-axis, a_len is the total lenthg (a_min-a_len/2 to a_min+a_len/2)
 *     b_min is the *middle* of the rectangle y-axis, b_len is the total lenthg (b_min-b_len/2 to b_min+b_len/2)
 *
 */
 
int make_pie( regRegion *reg, double a_min, double b_min, double a_len, double b_len);
int make_epanda( regRegion *reg, double a_min, double b_min, double a_len, double b_len);
int make_bpanda( regRegion *reg, double a_min, double b_min, double a_len, double b_len);
int make_rotbox( regRegion *reg, double a_min, double b_min, double a_len, double b_len);
void polar_limits( double a_min, double b_min, double a_len, double b_len, double pp[4],double qq[4] );
void cartesian_limits( double a_min, double b_min, double a_len, double b_len, double pp[4],double qq[4] );



/* ----------------------------- */


int convert_coords(dmDescriptor *xdesc, 
                  dmDescriptor *ydesc, 
                  double xx,       /*Note: Using double rather than long here! */
                  double yy,   /* */
                  double *xat, 
                  double *yat)
{
    if (xdesc) {
        double lgc[2];
        double phy[2];
        lgc[0] = xx + 1;
        lgc[1] = yy + 1;
        dmCoordCalc_d(xdesc, lgc, phy);
        if (ydesc) {
            dmCoordCalc_d(ydesc, lgc + 1, phy + 1);
        }
        *xat = phy[0];
        *yat = phy[1];
    } else {
        *xat = xx;
        *yat = yy;
    }
    return (0);
}


int invert_coords(dmDescriptor *xdesc, 
                 dmDescriptor *ydesc, 
                 double xx,        /*Note: Using double rather than long here! */
                 double yy,    /* */
                 double *ii, 
                 double *jj)
{
    if (xdesc) {
        double lgc[2];
        double phy[2];
        phy[0] = xx;
        phy[1] = yy;

        dmCoordInvert_d(xdesc, phy, lgc);
        if (ydesc) {
            dmCoordInvert_d(ydesc, phy + 1, lgc + 1);
        }
        *ii = lgc[0] - 1;
        *jj = lgc[1] - 1;
    } else {
        *ii = xx;
        *jj = yy;
    }
    return (0);
}

/* ---------------------------------------------- */

/* Which shapes use which parameters?
 * 
 *  | shape    |rstart |rstop |astart |astop |minrad |minang |ellip |
 *  |----------|-------|------|-------|------|-------|-------|------|
 *  | pie      | +     | +    |  +    | +    |  +    |  +    |  o   |
 *  | epanda   | +     | +    |  +    | +    |  +    |  +    |  +   |
 *  | bpanda   | +     | +    |  +    | +    |  +    |  +    |  +   |
 *  | box      | o     | +    |  +    | o    |  +    |  o    |  +   |
 * 
 */

double FUDGE_FACTOR = 0.00001;


/* technically a pie could be made using the epanda logic, but that runs
 * much slower so we use the dedicated shape instead
 */
int make_pie( regRegion *reg, 
               double a_min, 
               double b_min,
               double a_len, 
               double b_len
               )
{
    double reg_rad[2];
    double reg_ang[2];
    reg_rad[0] = a_min;
    reg_rad[1] = a_min + a_len;
    reg_ang[0] = b_min;
    reg_ang[1] = b_min + b_len;
    regAppendShape(reg, "Pie", 1, 1, &GlobalX0, &GlobalY0, 1, reg_rad,
                   reg_ang, 0, 0);
    return(0);
}


/*
 * Create a series of ellptical pie's (elliptical annulus over a fraction of angles)
 * 
 * ellipse(x0,y0,mjr_outer,mnr_outer,rotangle)*!ellipse(x0,y0,mjr_inner,mnr_inner,rotangle)*sector(x0,y0,start_ang,stop_ang)
 * 
 */
int make_epanda( regRegion *reg, 
               double a_min, 
               double b_min,
               double a_len, 
               double b_len
               )
{

    double reg_ang[2];
    reg_ang[0] = b_min;
    reg_ang[1] = b_min + b_len;

    double r[2];
    r[0] = a_min+a_len;
    r[1] = r[0] * GlobalEllipticity;
    regAppendShape(reg, "ellipse", 1, 1, &GlobalX0, &GlobalY0, 1, r,
                   &GlobalStartAngle, 0, 0);
    r[0] = a_min;
    r[1] = r[0] * GlobalEllipticity;
    regAppendShape(reg, "ellipse", 0, 0, &GlobalX0, &GlobalY0, 1, r,
                   &GlobalStartAngle, 0, 0);

    if ( b_len < 360.0 ) {
        regAppendShape(reg, "sector", 1, 0, &GlobalX0, &GlobalY0, 1, NULL,
                    reg_ang, 0, 0);
    }
    
    return(0);
}


/*
 * 
 * similar to above, but now using rotbox's instead of ellipses
 * 
 */
int make_bpanda( regRegion *reg, 
               double a_min, 
               double b_min,
               double a_len, 
               double b_len
               )
{
    double reg_ang[2];
    reg_ang[0] = b_min;
    reg_ang[1] = b_min + b_len;

    double r[2];
    r[0] = a_min+a_len;
    r[1] = r[0] * GlobalEllipticity;
    regAppendShape(reg, "rotbox", 1, 1, &GlobalX0, &GlobalY0, 1, r,
                   &GlobalStartAngle, 0, 0);
    r[0] = a_min;
    r[1] = r[0] * GlobalEllipticity;
    regAppendShape(reg, "rotbox", 0, 0, &GlobalX0, &GlobalY0, 1, r,
                   &GlobalStartAngle, 0, 0);

    if (b_len < 360.0) {
        regAppendShape(reg, "sector", 1, 0, &GlobalX0, &GlobalY0, 1, NULL,
                    reg_ang, 0, 0);
    }
    return(0);
}


/*
 *  For the cartesian grid we only use rotbox's
 */
int make_rotbox( regRegion *reg, 
               double a_min, 
               double b_min,
               double a_len, 
               double b_len
               )
{
    
    /*
     * Ugly ... but what I've seen is that due to numerical precision,
     * we can end up with regions where the 4 corners meet that land "exactly"
     * in the middle of a pixel then that pixel may never get included.
     * The limited precision in the location or lengths leaves the "gap".
     * 
     * Basically I need to increase the size of the regions by a tiny fraction
     * so that regions overlap ever so slightly.  That way they fully cover
     * the image.
     * 
     * Now, this does introduce some ambiguity when the pixel center lands
     * within this overlap region -- but we already have that (technically
     * the edge of the region is included, so abutting edges overlap). 
     * So I'm okay with this.
     */

    double ll[2] = { a_len+FUDGE_FACTOR, b_len+FUDGE_FACTOR };
    regAppendShape(reg, "Rotbox", 1, 1, &a_min, &b_min, 1, ll,
                   &GlobalStartAngle, 0, 0);
    return(0);
}


/*
 *  Create the region and return a bounding box around it.
 */
regRegion *make_region(regRegion *inreg, 
                   double a_min, 
                   double b_min,
                   double a_len, 
                   double b_len,  
                   long *xs, 
                   long *xl,
                   long *ys, 
                   long *yl)
{
    /*  */

    regRegion *reg;
    reg = inreg ? inreg : regCreateEmptyRegion();

    GlobalShapeFunction( reg, a_min, b_min, a_len, b_len );

    if (xs && ys && xl && yl) {
        double xrange[2], yrange[2];
        double hack_range[2] = { -1.0e60, 1.0e60 };       // The image better be fit.

        regExtent(reg, hack_range, hack_range, xrange, yrange);
        double ilo, jlo, ihi, jhi;
        invert_coords(GlobalXdesc, GlobalYdesc, xrange[0], yrange[0], &ilo,
                      &jlo);
        invert_coords(GlobalXdesc, GlobalYdesc, xrange[1], yrange[1], &ihi,
                      &jhi);

        /* Add some padding to the bounding box limits */
        /* INT/Truncate/etc check!!!!!!! */
        *xs = ilo - 2;
        *ys = jlo - 2;
        *xl = ihi - ilo + 4;
        *yl = jhi - jlo + 4;
    }



#ifdef DEMO_MODE    

    if (inreg) {
        // "epanda(4182.4908,3882.4997,44.993503,89.993503,1,100,75,400,300,1,19.993503)"

        fprintf( global_demo_file, "physical; epanda(%g,%g,%g,%g,1,%g,%g,%g,%g,1,%g)\n", 
          GlobalX0,GlobalY0, b_min-GlobalStartAngle, b_min+b_len-GlobalStartAngle, 
          a_min, a_min* GlobalEllipticity, 
          (a_min+a_len), (a_min+a_len)*GlobalEllipticity, GlobalStartAngle );
    }


#endif

    return reg;
}



/* Compute the signal to noise ratio in the sub-image.  Also returns the
 * sum of the pixel values and the area (number of non-null pixels) */
double get_snr(double a_min, 
               double b_min, 
               double a_len, 
               double b_len, 
               float *oval,     /* o: sum of pixel values */
               long *area       /* o: number of pixels */
    )
{
    long xs;                    /* i: start of x-axis (sub img) */
    long ys;                    /* i: start of y-axis (sub img) */
    long xl;                    /* i: length of x-axis (sub img) */
    long yl;                    /* i: length of y-axis (sub img) */

    double px, py;

    float val;
    float noise;
    float locsnr;
    long ii, jj;
    double pixval;

    val = 0.0;
    noise = 0.0;
    *area = 0;

    regRegion *reg = make_region(NULL, a_min, b_min, a_len, b_len, &xs, &xl, &ys, &yl);


    /* Determine SNR for current sub-image */
    for (ii = xs-1; ii <= (xl + xs)+1; ii++) {
        if (ii < 0 || ii >= GlobalXLen) {
            continue;
        }
        for (jj = ys-1; jj <= (yl + ys)+1; jj++) {
            long pix;
            if (jj < 0 || jj >= GlobalYLen) {
                continue;
            }

            convert_coords(GlobalXdesc, GlobalYdesc, ii, jj, &px, &py);
            if (0 == regInsideRegion(reg, px, py)) {
                continue;
            }

            pix = ii + (jj * GlobalXLen);

            if (0 != GlobalMask[pix]) {
                continue;
            }

            pixval = get_image_value(GlobalData, GlobalDataType, ii, jj,
                                     GlobalLAxes, GlobalPixMask);
            if (ds_dNAN(pixval)) {
                continue;
            }
            val += pixval;
            noise += (GlobalDErr[pix] * GlobalDErr[pix]);
            *area += 1;


        } // end for jj
    }  // end for ii
    
    regFree(reg);

    locsnr = val / sqrt(noise);
    *oval = val;

    return locsnr;
}



void fill_region(double a_min, 
                double b_min, 
                double a_len, 
                double b_len)
{
    long xs;                    /* i: start of x-axis (sub img) */
    long ys;                    /* i: start of y-axis (sub img) */
    long xl;                    /* i: length of x-axis (sub img) */
    long yl;                    /* i: length of y-axis (sub img) */

    static unsigned long mask_no;       /* keep as static to avoid yet another
                                           param to function; gets updated for
                                           each recursive call.  Could make this a global value too. */

    float locsnr;
    long ii, jj;
    float val;
    long area;
    double pixval;
    double px, py;

    regRegion *reg = make_region(NULL, a_min, b_min, a_len, b_len, &xs, &xl, &ys, &yl);

    locsnr = get_snr(a_min, b_min, a_len, b_len, &val, &area);

    if (0 == area)
        return;

    val /= area;

    mask_no += 1;               /* statically increases per bin */

    /* store output values */
    for (ii = xs-1; ii <= xs + xl+1; ii++) {
        if (ii < 0 || ii >= GlobalXLen) {
            continue;
        }
        for (jj = ys-1; jj <= ys + yl+1; jj++) {
            long pix;

            if (jj < 0 || jj >= GlobalYLen) {
                continue;
            }

            convert_coords(GlobalXdesc, GlobalYdesc, ii, jj, &px, &py);
            if (0 == regInsideRegion(reg, px, py)) {
                continue;
            }


            pix = ii + (jj * GlobalXLen);

            if (0 != GlobalMask[pix]) {
                continue;
            }

            pixval = get_image_value(GlobalData, GlobalDataType, ii, jj,
                                     GlobalLAxes, GlobalPixMask);

            if (ds_dNAN(pixval)) {
                GlobalOutData[pix] = pixval; // ie NaN
                GlobalOutArea[pix] = pixval;
                GlobalMask[pix] = 0;
                GlobalOutSNR[pix] = pixval;
            } else {
                GlobalOutData[pix] = val;
                GlobalOutArea[pix] = area;
                GlobalMask[pix] = mask_no;
                GlobalOutSNR[pix] = locsnr;
            }

        }                       // end for jj
    }                           // end for ii
    regFree(reg);


    /* Add to the global region definition */
    make_region(GlobalMaskRegion, a_min, b_min, a_len, b_len, NULL, NULL, NULL, NULL);


    return;

}


void polar_limits( double a_min, 
              double b_min,  
              double a_len, 
              double b_len,
              double pp[4],
              double qq[4] )
{
    double dpp;
    double dqq;
    dpp = a_len/2.0;
    dqq = b_len/2.0;
    pp[0] = a_min;     qq[0] = b_min;
    pp[1] = a_min+dpp; qq[1] = b_min;
    pp[2] = a_min;     qq[2] = b_min+dqq;
    pp[3] = a_min+dpp; qq[3] = b_min+dqq;        
}


void cartesian_limits( double a_min, 
              double b_min,  
              double a_len, 
              double b_len,
              double pp[4],
              double qq[4] )
{
    double dpp;
    double dqq;
    dpp = a_len/4.0;
    dqq = b_len/4.0;

    double rot_x, rot_y;
    double rad_ang = GlobalStartAngle*DEG2RAD;
    double c_a = cos(rad_ang);
    double s_a = sin(rad_ang);


    // rotate a_min,b_min backwards: note cos(a) = cos(-a) and 
    // -sin(a)=sin(-a)    
    double r0 =  (a_min-GlobalX0)*c_a + (b_min-GlobalY0)*s_a;
    double a0 = -(a_min-GlobalX0)*s_a + (b_min-GlobalY0)*c_a;

    // Lower left -- rotate back
    rot_x = r0-dpp;
    rot_y = a0-dqq;
    pp[0] = rot_x*c_a - rot_y*s_a + GlobalX0;
    qq[0] = rot_x*s_a + rot_y*c_a + GlobalY0;
    
    // Lower right -- rotate back
    rot_x = r0+dpp;
    rot_y = a0-dqq;
    pp[1] = rot_x*c_a - rot_y*s_a + GlobalX0;
    qq[1] = rot_x*s_a + rot_y*c_a + GlobalY0;

    // upper left -- rotate back
    rot_x = r0-dpp;
    rot_y = a0+dqq;
    pp[2] = rot_x*c_a - rot_y*s_a + GlobalX0;
    qq[2] = rot_x*s_a + rot_y*c_a + GlobalY0;

    // upper right -- rotate back
    rot_x = r0+dpp;
    rot_y = a0+dqq;
    pp[3] = rot_x*c_a - rot_y*s_a + GlobalX0;
    qq[3] = rot_x*s_a + rot_y*c_a + GlobalY0;


}


/* Recursive binning routine */
void abin_rec(double a_min, 
              double b_min,  
              double a_len, 
              double b_len)
{

    double pp[4], qq[4], dpp, dqq;

    GlobalLimitsFunction( a_min, b_min, a_len, b_len, pp, qq );

    dpp = a_len/2.0;
    dqq = b_len/2.0;


    short check = 0;
    if (ZERO_ABOVE == GlobalSplitCriteria) {
        /* This is the original method -- if the current block is above SNR
         * then split it. */
        float at, oval;
        long npix;
        at = get_snr(a_min, b_min, a_len, b_len, &oval, &npix);
        check = (at > GlobalSNRThresh);

    } else {
        /* Determine SNR for current sub-image */

        /* This new method will only split the 2x2 if when the sub images
         * are created, some/all of them will remain above the SNR limit.
         * 
         */

        float ll, lr, ul, ur;
        float oval_ll, oval_lr, oval_ul, oval_ur;
        long npix_ll, npix_lr, npix_ul, npix_ur;

        short ill, ilr, iul, iur;

        ll = get_snr(pp[0], qq[0], dpp, dqq, &oval_ll, &npix_ll); /* low-left */
        lr = get_snr(pp[1], qq[1], dpp, dqq, &oval_lr, &npix_lr); /* low-rite */
        ul = get_snr(pp[2], qq[2], dpp, dqq, &oval_ul, &npix_ul); /* up-left */
        ur = get_snr(pp[3], qq[3], dpp, dqq, &oval_ur, &npix_ur); /* up-rite */

        /*
         * It is OK to split if sub-cell has no valid pixel; but not all of them.
         */
        ill = (ll >= GlobalSNRThresh) || (npix_ll == 0);
        ilr = (lr >= GlobalSNRThresh) || (npix_lr == 0);
        iul = (ul >= GlobalSNRThresh) || (npix_ul == 0);
        iur = (ur >= GlobalSNRThresh) || (npix_ur == 0);

        /* If there are no pixels, no reason to recurse */
        if ((npix_ll + npix_lr + npix_ul + npix_ur) == 0) {
            check = 0;

        } else if (ONE_ABOVE == GlobalSplitCriteria) {
            /* If any one of the sub images is above snr, then split */
            check = (ill + ilr + iul + iur);

        } else if (TWO_ABOVE == GlobalSplitCriteria) {
            /* If two, then they have to be side-by side, not diagonal */
            check = ((ill && ilr) || (ilr && iur) || (iur && iul) || (iul && ill));

        } else if (THREE_ABOVE == GlobalSplitCriteria) {
            check = (ill + ilr + iul + iur) >= 3 ? 1 : 0;

        } else if (ALL_ABOVE == GlobalSplitCriteria) {
            check = (ill + ilr + iul + iur) == 4 ? 1 : 0;

        } else {
            err_msg("This should not have happened, something's amiss");
            return;
        }

    }                           // end else 


    if ((check) && (a_len > GlobalMinRadius) && (b_len > GlobalMinAngle)) {
        /* Enter recursion */
        
        abin_rec(pp[0],qq[0],dpp,dqq); /* low-left */
        abin_rec(pp[1],qq[1],dpp,dqq); /* low-rite */
        abin_rec(pp[2],qq[2],dpp,dqq); /* up-left */
        abin_rec(pp[3],qq[3],dpp,dqq); /* up-rite */
        return;
    }

    fill_region(a_min, b_min, a_len, b_len);

}



int load_error_image(char *errimg)
{

    /* Read Error Image */
    unsigned long npix = GlobalXLen * GlobalYLen;

    if ((strlen(errimg) == 0) || (ds_strcmp_cis(errimg, "none") == 0)) {

        double pixval;
        long xx, yy, jj;

        for (yy = 0; yy < GlobalYLen; yy++) {
            for (xx = 0; xx < GlobalXLen; xx++) {
                jj = xx + yy * GlobalXLen;
                pixval = get_image_value(GlobalData, GlobalDataType, xx, yy,
                                    GlobalLAxes, GlobalPixMask);

                if (ds_dNAN(pixval)) {
                    GlobalDErr[jj] = 0;
                } else {
                    GlobalDErr[jj] = sqrt(pixval);      // assumes Gaussian stats
                }
            }                   // end xx
        }                       // end yy

    } else {
        long enAxes;
        long *elAxes;
        dmDescriptor *errDs;
        dmBlock *erBlock;

        erBlock = dmImageOpen(errimg);
        if (erBlock == NULL) {
            err_msg("ERROR: Could not open file '%s'\n", errimg);
            return (-1);
        }
        errDs = dmImageGetDataDescriptor(erBlock);
        enAxes = dmGetArrayDimensions(errDs, &elAxes);
        if ((enAxes != 2) ||
            (elAxes[0] != GlobalXLen) || (elAxes[1] != GlobalYLen)) {
            err_msg
                ("ERROR: Error image must be 2D image with non-zero axes\n");
            return (-1);
        }
        // get_data( errDs, npix, GlobalDErr );
        dmGetArray_f(errDs, GlobalDErr, npix);
        dmImageClose(erBlock);
    }

    return 0;
}


int write_single_output(dmBlock *inBlock, char *outfile, void *outdata, short clobber, dmDataType dt, 
    regRegion *outreg)
{
    dmBlock *outBlock;
    dmDescriptor *outDes;

    if ((strlen(outfile) == 0) || (ds_strcmp_cis(outfile, "none") == 0)) {
        return(0);
    }

    outBlock = dmImageCreate(outfile, dt, GlobalLAxes, 2);
    if (outBlock == NULL) {
        err_msg("ERROR: Could not create output '%s'\n", outfile);
    }

    outDes = dmImageGetDataDescriptor(outBlock);
    dmBlockCopy(inBlock, outBlock, "HEADER");
    ds_copy_full_header(inBlock, outBlock, "dmradar", 0);
    put_param_hist_info(outBlock, "dmradar", NULL, 0);
    dmBlockCopyWCS(inBlock, outBlock);

    if ( dmFLOAT == dt )  {
        dmSetArray_f(outDes, (float*)outdata, GlobalLAxes[0]*GlobalLAxes[1]);
    } else if (dmULONG == dt ) {
        dmSetArray_ul(outDes, (unsigned long*)outdata, GlobalLAxes[0]*GlobalLAxes[1]);
    } else {
        return(-1);
    }

    if ( outreg ) {
        dmBlockClose(dmTableWriteRegion(dmBlockGetDataset(outBlock),
                                        "REGION", NULL, outreg));
    }

    dmImageClose(outBlock);

    return(0);
    
}



dmBlock *load_infile( char *infile)
{
    /* Read the data */
    dmBlock *inBlock;
 
    inBlock = dmImageOpen(infile);
    if (!inBlock) {
        err_msg("ERROR: Could not open infile='%s'\n", infile);
        return (NULL);
    }

    regRegion *dss = NULL;
    long null;
    short has_null;
    GlobalDataType = get_image_data(inBlock, &GlobalData, &GlobalLAxes, &dss, &null,
                       &has_null);
    get_image_wcs(inBlock, &GlobalXdesc, &GlobalYdesc);
    GlobalPixMask = get_image_mask(inBlock, GlobalData, GlobalDataType, GlobalLAxes, dss,
                       null, has_null, GlobalXdesc, GlobalYdesc);
    GlobalXLen = GlobalLAxes[0];
    GlobalYLen = GlobalLAxes[1];

    return(inBlock);    
}


int allocate_memory(void)
{
    long npix;
    npix = (GlobalLAxes[0] * GlobalLAxes[1]);
    if (npix == 0) {
        err_msg("ERROR: Image is empty (one axis is 0 length)\n");
        return (-1);
    }

    /* Allocate memory for the products */
    if (( NULL == (GlobalDErr = (float *) calloc(npix, sizeof(float)))) ||
        ( NULL == (GlobalOutData = (float *) calloc(npix, sizeof(float)))) ||
        ( NULL == (GlobalOutArea = (float *) calloc(npix, sizeof(float)))) ||
        ( NULL == (GlobalOutSNR = (float *) calloc(npix, sizeof(float)))) ||
        ( NULL == (GlobalMask = (unsigned long *) calloc(npix, sizeof(unsigned long)))) ||
        ( NULL == (GlobalMaskRegion = regCreateEmptyRegion())) ) {
        err_msg("ERROR: Problem allocating memory");
        return(-34);            
    }
    
    return(0);
}



/* help for autoname f() below for the aux files */
int aux_output(char *outfile, char *auxoutfile, char *suffix, short clobber )
{
    ds_autoname(outfile, auxoutfile, suffix, DS_SZ_FNAME);
    if (strlen(auxoutfile)>0 && 0!=ds_strcmp_cis(auxoutfile,"none")) {
        if (ds_clobber(auxoutfile, clobber, NULL) != 0) {
            return(1);
        }        
    }
    return(0);
}


/* Setup the output file names and clobber if they exist */
int autoname( Parameters *pp)
{
    long maxlen = DS_SZ_FNAME;

    if ( 0 == strlen(pp->infile)) {
        err_msg("ERROR: infile cannot be blank");
        return(-1);
    }

    /* Go ahead and take care of the autonaming */
    ds_autoname(pp->infile, pp->outfile, "abinimg", maxlen);

    if ( (0 == strlen(pp->outfile)) || ( 0 == ds_strcmp_cis(pp->outfile, "none"))) {
        err_msg("ERROR: outfile cannot be blank");
        return(-2);
    }

    if (ds_clobber(pp->outfile, pp->clobber, NULL) != 0) {
        return(1);
    }

    if (( 0 != aux_output( pp->outfile, pp->maskfile, "maskimg", pp->clobber)) ||
        ( 0 != aux_output( pp->outfile, pp->snrfile, "snrimg", pp->clobber)) ||
        ( 0 != aux_output( pp->outfile, pp->areafile, "areaimg", pp->clobber))) {
        return(1);
    }


    return(0);
}


/* Determine which algorithm to use based on the number of sub-images
 * that need to be *above* SNR limit in the final image.
 * 
 * 0 is the orignal dmnautilus algorithm, where it keeps trying to split
 * until pixels are below SNR
 * 
 * '4' requires all 4 sub-images to be above SNR limit.
 */
int map_method( short method ) 
{

    switch (method) {
    case 0:  GlobalSplitCriteria = ZERO_ABOVE;   break;
    case 1:  GlobalSplitCriteria = ONE_ABOVE;    break;
    case 2:  GlobalSplitCriteria = TWO_ABOVE;    break;
    case 3:  GlobalSplitCriteria = THREE_ABOVE;  break;
    case 4:  GlobalSplitCriteria = ALL_ABOVE;    break;
    default:
        err_msg("Invalid method parameter value");
        return (-1);
        break;
    };
    return(0);
}


void free_globals(void)
{
    /* make valgrind happy */
    free(GlobalData);
    free(GlobalDErr);
    free(GlobalOutData);
    free(GlobalOutArea);
    free(GlobalOutSNR);
    free(GlobalMask);

}


int write_outputs( dmBlock *inBlock, Parameters *pp)
{
    /* Write out files -- NB: mask file has different datatypes and different extensions */
    if ( ( 0 != write_single_output( inBlock, pp->outfile, GlobalOutData, pp->clobber, dmFLOAT, NULL )) ||
         ( 0 != write_single_output( inBlock, pp->areafile, GlobalOutArea, pp->clobber, dmFLOAT, NULL )) ||    
         ( 0 != write_single_output( inBlock, pp->snrfile, GlobalOutSNR, pp->clobber, dmFLOAT, NULL )) ||
         ( 0 != write_single_output( inBlock, pp->maskfile, GlobalMask, pp->clobber, dmULONG, GlobalMaskRegion))) {
        return(-1);
     }

    return(0);
}

Parameters *load_parameters(void)
{
    Parameters *par = (Parameters*)calloc(1,sizeof(Parameters));
    if ( NULL == par ) {
        return(NULL);
    }
    
    /* Read in all the data */
    clgetstr("infile", par->infile, DS_SZ_FNAME);
    clgetstr("outfile", par->outfile, DS_SZ_FNAME);
    GlobalSNRThresh = clgetd("snr");
    GlobalX0 = clgetd("xcenter");
    GlobalY0 = clgetd("ycenter");
    par->method = clgeti("method");
    clgetstr( "shape", par->shape, 99 );

    GlobalInnerRadius = clgetd("rstart");
    GlobalOuterRadius = clgetd("rstop");
    GlobalStartAngle = clgetd("astart");
    GlobalStopAngle = clgetd("astop");

    GlobalEllipticity = clgetd("ellipticity");

    GlobalMinRadius = clgetd("minradius");
    GlobalMinAngle = clgetd("minangle");

    clgetstr("inerrfile", par->errimg, DS_SZ_FNAME);
    clgetstr("outmaskfile", par->maskfile, DS_SZ_FNAME);
    clgetstr("outsnrfile", par->snrfile, DS_SZ_FNAME);
    clgetstr("outareafile", par->areafile, DS_SZ_FNAME);
    par->clobber = clgetb("clobber");

    return(par);
}


int map_shapes( char *shape )
{
    if ( 0 ==  strcmp(shape, "pie") ) {
        GlobalShapeFunction = make_pie;
        GlobalLimitsFunction = polar_limits;        
    } else if ( 0 == strcmp(shape, "epanda")) {
        GlobalShapeFunction = make_epanda;
        GlobalLimitsFunction = polar_limits;                
    } else if ( 0 == strcmp(shape, "bpanda")) {
        GlobalShapeFunction = make_bpanda;
        GlobalLimitsFunction = polar_limits;                
    } else if ( 0 == strcmp(shape, "box")) {
        GlobalShapeFunction = make_rotbox;
        GlobalLimitsFunction = cartesian_limits;
    } else {
        err_msg("ERROR: Unknown shape='%s'", shape);
        return(-1);
    }

    return(0);
    
}



/* Main routine; does all the work of a quad-tree adaptive binning routine*/
int abin(void)
{

    Parameters *pp;
    if ( NULL == (pp = load_parameters())) {
        return(-1);
    }

    if ( 0 != map_method(pp->method) ) {
        return(-1);
    }   

    if ( 0 != map_shapes(pp->shape)) {
        return(-1);
    }

    if ( 0 != autoname(pp)) {
        return(-1);
    }

    dmBlock *inBlock;
    if ( NULL == (inBlock = load_infile( pp->infile))) {
        return(-1);
    }

    if ( 0 != allocate_memory() ) {
        return(-1);
    }

    if (0 != load_error_image(pp->errimg)) {
        return -1;
    }


#ifdef DEMO_MODE
    global_demo_file = fopen("demo.reg", "w");
#endif



    /* Start Algorithm */
    if ( polar_limits == GlobalLimitsFunction ) {
        if (GlobalInnerRadius > 0) {
            fill_region(0, GlobalStartAngle, GlobalInnerRadius, GlobalStopAngle);
        }
        abin_rec(GlobalInnerRadius, GlobalStartAngle, GlobalOuterRadius,
                 GlobalStopAngle);
    } else {
        abin_rec( GlobalX0, GlobalY0, GlobalOuterRadius, GlobalOuterRadius*GlobalEllipticity);        
    }

    if ( 0 != write_outputs( inBlock, pp)) {
        return(-1);
    }


#ifdef DEMO_MODE
    fclose( global_demo_file);
#endif


    /* Must keep open until now to do all the wcs/hdr copies */
    dmImageClose(inBlock);
    free_globals();
    return (0);

}
