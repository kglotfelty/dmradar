/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
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

regRegion *maskRegion;


/* Okay, global variables are a bad idea, excpet when doing recursion and
 * we are passing the same data to each level ... more data gets pushed on the
 * stack and causes crashes.  Also using global variables can make it
 * run much faster so we'll bite the bullet here. */
void  *GlobalData;     /* i: data array */
float *GlobalDErr;     /* i: error array */
float *GlobalOutData;  /* o: output array */
float *GlobalOutArea;  /* o: output area array */
float *GlobalOutSNR;   /* o: output SNR */
unsigned long *GlobalMask; /* o: output mask */
long GlobalXLen;        /* i: length of x-axis (full img) */
long GlobalYLen;        /* i: length of y-axis (full img) */
long GlobalLAxes[2];    /* X,Y lens togeether */
float GlobalSNRThresh;   /* i: SNR threshold */
enum { ZERO_ABOVE=0, ONE_ABOVE, TWO_ABOVE, THREE_ABOVE, ALL_ABOVE} GlobalSplitCriteria;
dmDataType GlobalDataType;
dmDescriptor *GlobalXdesc=NULL;
dmDescriptor *GlobalYdesc=NULL;
short *GlobalPixMask=NULL;


double GlobalX0; // =  4274.5; // 4001.9; //4274.0;
double GlobalY0; // =  3956.0; // 3850.5; // 3954.0;
double GlobalInnerRadius; // = 5.0;
double GlobalOuterRadius; // = 1000.0;
double GlobalStartAngle; // = 0.0;
double GlobalStopAngle; // = 360.0;
double GlobalMinRadius; // = 0.5;
double GlobalMinAngle; // = 1 ;





/* Using the dmtools/dmimgio routines removes lots of duplicate code that was
 * originally here.  Also allow us to keep track of NULL/NaN value pixels
 * more easily
 */
#include "dmimgio.h"



/* ------Prototypes ----------------------- */

int load_error_image( char *errimg );
int abin(void);

double get_snr(double r_min, double a_min, double r_len, double a_len, float *oval, long *area);
void abin_rec (double r_min, double a_min, double r_len, double a_len);   
int convert_coords( dmDescriptor *xdesc, dmDescriptor *ydesc, double xx, double yy, double *xat, double *yat);
int invert_coords( dmDescriptor *xdesc, dmDescriptor *ydesc, double xat, double yat, double *ii, double *jj);
regRegion *make_pie( regRegion *inreg, double r_min, double r_len, double a_min, double a_len, long *xs, long *xl, long *ys, long *yl);
void fill_region( double r_min,double a_min,double r_len,double a_len) ;

/* ----------------------------- */


int convert_coords( dmDescriptor *xdesc,
                    dmDescriptor *ydesc,
                    double xx,  /*Note: Using double rather than long here! */
                    double yy, /* */
                    double *xat,
                    double *yat
                    )
{
  if ( xdesc ) {
    double lgc[2];
    double phy[2];
    lgc[0] = xx+1;
    lgc[1] = yy+1;
    dmCoordCalc_d( xdesc, lgc, phy );
    if ( ydesc ) {
      dmCoordCalc_d( ydesc, lgc+1, phy+1 );
    }
    *xat = phy[0];
    *yat = phy[1];
  } else {
    *xat = xx;
    *yat = yy;
  }
  return(0);
}


int invert_coords( dmDescriptor *xdesc,
                    dmDescriptor *ydesc,
                    double xx,  /*Note: Using double rather than long here! */
                    double yy, /* */
                    double *ii,
                    double *jj
                    )
{
  if ( xdesc ) {
    double lgc[2];
    double phy[2];
    phy[0] = xx;
    phy[1] = yy;

    dmCoordInvert_d( xdesc, phy, lgc);
    if ( ydesc ) {
      dmCoordInvert_d( ydesc, phy+1, lgc+1);
    }
    *ii = lgc[0]-1;
    *jj = lgc[1]-1;
  } else {
    *ii = xx;
    *jj = yy;
  }
  return(0);
}





regRegion *make_pie( regRegion *inreg, double r_min, double r_len, double a_min, double a_len, long *xs, long *xl, long *ys, long *yl)
{
      /*  */

  regRegion *reg;
  reg = inreg ? inreg : regCreateEmptyRegion();

  if ( 1 ) {
      double reg_rad[2];
      double reg_ang[2];  
      reg_rad[0] = r_min;
      reg_rad[1] = r_min + r_len;
      reg_ang[0] = a_min;
      reg_ang[1] = a_min + a_len;
      regAppendShape( reg, "Pie", 1, 1, &GlobalX0, &GlobalY0, 1, reg_rad, reg_ang, 0, 0 );
  } else {
      double xx = r_min+r_len/2.0;
      double yy = a_min+a_len/2.0;
      double ll[2] = { r_len, a_len };
      regAppendShape( reg, "Rotbox", 1, 1, &xx, &yy, 1, ll, &GlobalStartAngle, 0, 0 );
  }  


  if ( xs && ys && xl && yl ) {
      double xrange[2], yrange[2];
      double hack_range[2] = { -1.0e6, 1.0e6}; // The image better be fit.

      regExtent( reg, hack_range, hack_range, xrange, yrange );
      double ilo, jlo, ihi, jhi;
      invert_coords( GlobalXdesc, GlobalYdesc, xrange[0], yrange[0], &ilo, &jlo );
      invert_coords( GlobalXdesc, GlobalYdesc, xrange[1], yrange[1], &ihi, &jhi );

      /* Add some padding to the bounding box limits */
      /* INT/Truncate/etc check!!!!!!! */
      *xs = ilo - 2;  
      *ys = jlo - 2;
      *xl = ihi-ilo + 4;
      *yl = jhi-jlo + 4; 
  }
  
  return reg;
}



/* Compute the signal to noise ratio in the sub-image.  Also returns the
 * sum of the pixel values and the area (number of non-null pixels) */
double get_snr( 
        double r_min,
        double a_min,
        double r_len,
        double a_len,
        float  *oval,   /* o: sum of pixel values */
        long   *area    /* o: number of pixels */
  )
{
  long   xs;     /* i: start of x-axis (sub img) */
  long   ys;     /* i: start of y-axis (sub img) */
  long   xl;     /* i: length of x-axis (sub img) */
  long   yl;      /* i: length of y-axis (sub img) */

  double px, py;
  
  float val;
  float noise;
  float locsnr;
  long ii,jj;
  double pixval;

  val = 0.0;
  noise = 0.0;
  *area = 0;

  regRegion *reg = make_pie( NULL, r_min, r_len, a_min, a_len, &xs, &xl, &ys, &yl );


  /* Determine SNR for current sub-image */
  for (ii=xs; ii<(xl+xs); ii++ ) {
      if ( ii< 0 || ii >= GlobalXLen ) {
        continue; 
      }
    for (jj=ys; jj<(yl+ys); jj++) {
      long pix;
      if ( jj<0 || jj >= GlobalYLen ) {
        continue; 
      }
            
     convert_coords( GlobalXdesc,GlobalYdesc, ii, jj, &px, &py );
     if ( 0 ==  regInsideRegion( reg, px, py ) ) {
        continue;
     }
     

      pix = ii+(jj*GlobalXLen);

      if ( 0 != GlobalMask[pix] ) {
        continue;
      }

      pixval = get_image_value( GlobalData, GlobalDataType, ii, jj,
          GlobalLAxes, GlobalPixMask);
      if ( ds_dNAN(pixval) ) {          
          continue;
      }
      val += pixval;
      noise += ( GlobalDErr[pix] * GlobalDErr[pix] );
      *area += 1;


    }
  }


  locsnr = val / sqrt(noise);
  *oval = val;

  return locsnr;
}



void fill_region( 
        double r_min,
        double a_min,
        double r_len,
        double a_len
        )   
{
  long   xs;     /* i: start of x-axis (sub img) */
  long   ys;     /* i: start of y-axis (sub img) */
  long   xl;     /* i: length of x-axis (sub img) */
  long   yl;      /* i: length of y-axis (sub img) */

  regRegion *reg = make_pie( NULL, r_min, r_len, a_min, a_len, &xs, &xl, &ys, &yl );

    static unsigned long mask_no; /* keep as static to avoid yet another
                   param to function; gets updated for
                   each recursive call.  Could make this a global value too.*/

    float locsnr;
    long ii,jj;
    float val;
    long area;
    double pixval;
    double px, py;

    locsnr = get_snr( r_min, a_min, r_len, a_len, &val, &area );

    if ( 0 == area )
      return;
    
    val /= area;

    mask_no += 1; /* statically increases per bin */

    /* store output values */
    for (ii=xs; ii<xs+xl; ii++ ) {
      if ( ii<0 || ii >= GlobalXLen ) { 
      continue; 
      }
      for (jj=ys; jj<ys+yl; jj++) {
      long pix;

      if ( jj<0 || jj >= GlobalYLen ) { 
        continue; 
      }

     convert_coords( GlobalXdesc,GlobalYdesc, ii, jj, &px, &py );
     if ( 0 ==  regInsideRegion( reg, px, py ) ) {
        continue;
     }


    pix = ii+(jj*GlobalXLen);

      if ( 0 != GlobalMask[pix] ) {
        continue;
      }

      pixval = get_image_value( GlobalData, GlobalDataType, ii, jj,
          GlobalLAxes, GlobalPixMask);

      if ( ds_dNAN(pixval) ) {
        GlobalOutData[pix] = pixval;
        GlobalOutArea[pix] = pixval;
        GlobalMask[pix] = 0;
        GlobalOutSNR[pix] = pixval;
      } else {
        GlobalOutData[pix] = val;
        GlobalOutArea[pix] = area;
        GlobalMask[pix] = mask_no;
        GlobalOutSNR[pix] = locsnr;
      }
      
     } // end for jj
    } // end for ii
    regFree(reg);


    /* Add to the global region definition */   
    make_pie( maskRegion, r_min, r_len, a_min, a_len, NULL, NULL, NULL, NULL );


    return;

}



/* Recursive binning routine */
void abin_rec ( 
        double r_min,
        double a_min,
        double r_len,
        double a_len
        )   
{


  long   xs;     /* i: start of x-axis (sub img) */
  long   ys;     /* i: start of y-axis (sub img) */
  long   xl;     /* i: length of x-axis (sub img) */
  long   yl;      /* i: length of y-axis (sub img) */

  regRegion *reg = make_pie( NULL, r_min, r_len, a_min, a_len, &xs, &xl, &ys, &yl );
  regFree(reg);  // We just need the bounds 
  

  short check = 0;
  if ( ZERO_ABOVE == GlobalSplitCriteria ) {
      /* This is the original method -- if the current block is above SNR
       * then split it. */
      float at, oval;
      long npix;  
      at = get_snr( r_min, a_min, r_len, a_len , &oval, &npix );
      check = ( at > GlobalSNRThresh );

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

      ll = get_snr( r_min, a_min, FLOOR(r_len/2.0), FLOOR(a_len/2.0), &oval_ll, &npix_ll ); /* low-left */
      lr = get_snr( r_min+FLOOR(r_len/2.0), a_min, CEIL(r_len/2.0), FLOOR(a_len/2.0), &oval_lr, &npix_lr ); /* low-rite*/
      ul = get_snr( r_min, a_min+FLOOR(a_len/2.0), FLOOR(r_len/2.0), CEIL(a_len/2.0), &oval_ul, &npix_ul); /* up-left */
      ur = get_snr( r_min+FLOOR(r_len/2.0), a_min+FLOOR(a_len/2.0), CEIL(r_len/2.0), CEIL(a_len/2.0), &oval_ur, &npix_ur ); /* up-rite */

      /*
       * It is OK to split if sub-cell has no valid pixel; but not all of them.
       */       
      ill = ( ll >= GlobalSNRThresh) || ( npix_ll == 0);  
      ilr = ( lr >= GlobalSNRThresh) || ( npix_lr == 0);
      iul = ( ul >= GlobalSNRThresh) || ( npix_ul == 0);
      iur = ( ur >= GlobalSNRThresh) || ( npix_ur == 0);

      /* If there are no pixels, no reason to recurse */
      if ((npix_ll+npix_lr+npix_ul+npix_ur) ==0 ) {
          check =0; 

      } else if ( ONE_ABOVE == GlobalSplitCriteria ) {
          /* If any one of the sub images is above snr, then split */
          check = ( ill + ilr + iul + iur  );            

      } else if ( TWO_ABOVE == GlobalSplitCriteria ) {
          /* If two, then they have to be side-by side, not diagonal */
         check = ( (ill && ilr ) || 
                   (ilr && iur ) ||
                   (iur && iul ) ||
                   (iul && ill )
                  ); 

      } else if ( THREE_ABOVE == GlobalSplitCriteria ) {
          check = ( ill + ilr + iul + iur ) >= 3 ? 1 : 0;

      } else if ( ALL_ABOVE == GlobalSplitCriteria ) {
          check = ( ill + ilr + iul + iur ) == 4 ? 1 : 0;

      } else {
          err_msg("This should not have happened, something's amiss");
          return;
      }

  } // end else 

  
  if (( check )&&(r_len>GlobalMinRadius)&&(a_len>GlobalMinAngle)) {
    /* Enter recursion */

    /* need to use floor() and ceil() because input image may
       not be square, or 2^n.  This will bias the left-upper image w/ 1 more 
       pixel per bin; but that's a limit of not using square images.
       The alternative is only use square smallest sub-image or pad image to 2**N.*/

        abin_rec( r_min, a_min, FLOOR(r_len/2.0), FLOOR(a_len/2.0) ); /* low-left */
        abin_rec( r_min+FLOOR(r_len/2.0), a_min, CEIL(r_len/2.0), FLOOR(a_len/2.0) ); /* low-rite*/
        abin_rec( r_min, a_min+FLOOR(a_len/2.0), FLOOR(r_len/2.0), CEIL(a_len/2.0) ); /* up-left */
        abin_rec( r_min+FLOOR(r_len/2.0), a_min+FLOOR(a_len/2.0), CEIL(r_len/2.0), CEIL(a_len/2.0)); /* up-rite */
        return; 
    }
      
    fill_region( r_min, a_min, r_len, a_len );
  
}



int load_error_image( char *errimg ) {
    
  /* Read Error Image */
  unsigned long npix = GlobalXLen*GlobalYLen;

  if ( ( strlen(errimg) == 0 ) ||
       ( ds_strcmp_cis(errimg,"none" ) == 0 ) ) {

     double pixval;
     long xx,yy,jj;    

     for (yy=0; yy<GlobalYLen; yy++) {
        for ( xx=0;xx<GlobalXLen;xx++) {

            jj = xx + yy*GlobalXLen;
           pixval = get_image_value( GlobalData, GlobalDataType, xx,yy,
                    GlobalLAxes, GlobalPixMask);

            if (ds_dNAN(pixval) ) {
                GlobalDErr[jj] = 0;
            } else {
                GlobalDErr[jj] = sqrt(pixval);  // assumes Gaussian stats
            }
        
       }  // end xx
     } // end yy

  } else {
    long enAxes;
    long *elAxes;
    dmDescriptor *errDs;
    dmBlock *erBlock;

    erBlock = dmImageOpen( errimg );
    if ( erBlock == NULL ) {
      err_msg("ERROR: Could not open file '%s'\n", errimg);
      return(-1);
    }
    errDs = dmImageGetDataDescriptor(erBlock );
    enAxes = dmGetArrayDimensions( errDs, &elAxes );
    if ( (enAxes != 2 ) || 
     ( elAxes[0] != GlobalXLen ) ||
     ( elAxes[1] != GlobalYLen )    ) {
      err_msg("ERROR: Error image must be 2D image with non-zero axes\n");
      return(-1);
    }
    // get_data( errDs, npix, GlobalDErr );

    dmGetArray_f( errDs, GlobalDErr, npix );

    dmImageClose( erBlock );
  }

  return 0;
}



/* Main routine; does all the work of a quad-tree adaptive binning routine*/
int abin (void)
{

  char infile[DS_SZ_FNAME];
  char errimg[DS_SZ_FNAME];

  char outfile[DS_SZ_FNAME];
  char areafile[DS_SZ_FNAME];
  char maskfile[DS_SZ_FNAME];
  char snrfile[DS_SZ_FNAME];
  short method;
  short clobber;

  long npix;

  dmBlock *inBlock;

  /* Read in all the data */
  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  GlobalSNRThresh = clgetd( "snr" );
  GlobalX0 = clgetd("xcenter");
  GlobalY0 = clgetd("ycenter");
  method = clgeti("method");

  GlobalInnerRadius = clgetd("rstart");
  GlobalOuterRadius = clgetd("rstop");
  GlobalStartAngle = clgetd("astart");
  GlobalStopAngle = clgetd("astop");
  GlobalMinRadius = clgetd("minradius");
  GlobalMinAngle = clgetd("minangle");


  clgetstr( "inerrfile",   errimg, DS_SZ_FNAME );
  clgetstr( "outmaskfile", maskfile, DS_SZ_FNAME );
  clgetstr( "outsnrfile",  snrfile,  DS_SZ_FNAME );
  clgetstr( "outareafile", areafile, DS_SZ_FNAME );
  clobber = clgetb( "clobber" );

  switch (method) 
  {
    case 0: GlobalSplitCriteria = ZERO_ABOVE; break;
    case 1: GlobalSplitCriteria = ONE_ABOVE; break;
    case 2: GlobalSplitCriteria = TWO_ABOVE; break;
    case 3: GlobalSplitCriteria = THREE_ABOVE; break;
    case 4: GlobalSplitCriteria = ALL_ABOVE; break;
    default:
      err_msg("Invalid method parameter value");
      return(-1);
      break;
  };


  /* Go ahead and take care of the autonaming */
  ds_autoname( infile, outfile, "abinimg", DS_SZ_FNAME );
  ds_autoname( outfile, maskfile, "maskimg", DS_SZ_FNAME );
  ds_autoname( outfile, snrfile, "snrimg", DS_SZ_FNAME );
  ds_autoname( outfile, areafile, "areaimg", DS_SZ_FNAME );


  /* Read the data */
  /* TODO: Replace with dmimgio routines */

  inBlock = dmImageOpen( infile );
  if ( !inBlock ) {
    err_msg("ERROR: Could not open infile='%s'\n", infile );
    return(-1);
  }



  long *lAxes=NULL;
  regRegion *dss=NULL;
  long null;
  short has_null;

  GlobalDataType = get_image_data( inBlock, &GlobalData, &lAxes, &dss, &null, &has_null );
  get_image_wcs( inBlock, &GlobalXdesc, &GlobalYdesc );
  GlobalPixMask = get_image_mask( inBlock, GlobalData, GlobalDataType, lAxes, dss, null, has_null, 
                         GlobalXdesc, GlobalYdesc );


  npix = ( lAxes[0]*lAxes[1]);
  if (npix ==0 ) {
    err_msg("ERROR: Image is empty (one axis is 0 length)\n");
    return(-1);
  }


  GlobalLAxes[0] = GlobalXLen = lAxes[0];
  GlobalLAxes[1] = GlobalYLen = lAxes[1];


  /* Allocate memory for the products */
  //GlobalData = (float*)calloc(npix,sizeof(float));
  GlobalDErr = (float*)calloc(npix,sizeof(float));
  GlobalOutData = (float*)calloc(npix,sizeof(float));
  GlobalOutArea = (float*)calloc(npix,sizeof(float));
  GlobalOutSNR = (float*)calloc(npix,sizeof(float));
  GlobalMask = (unsigned long*)calloc(npix,sizeof(unsigned long));
  maskRegion = regCreateEmptyRegion();

  if ( 0 != load_error_image( errimg ) ) {
        return -1;
  }
  
  /* Start Algorithm */


  if ( GlobalInnerRadius > 0 ) {
    fill_region( 0, GlobalStartAngle, GlobalInnerRadius, GlobalStopAngle);
  }
  abin_rec( GlobalInnerRadius, GlobalStartAngle, GlobalOuterRadius, GlobalStopAngle);


  /* Write out files -- NB: mask file has different datatypes and different extensions */
  dmBlock *outBlock;
  dmDescriptor *outDes;

  if ( ds_clobber( outfile, clobber, NULL) == 0 ) {
    outBlock = dmImageCreate(outfile, dmFLOAT, lAxes, 2 );
    if ( outBlock == NULL ) {
      err_msg("ERROR: Could not create output '%s'\n", outfile);
    }
    outDes = dmImageGetDataDescriptor( outBlock );
    dmBlockCopy( inBlock, outBlock, "HEADER"); 
    ds_copy_full_header( inBlock, outBlock, "dmradar", 0 );
    put_param_hist_info( outBlock, "dmradar", NULL, 0 );
    dmBlockCopyWCS( inBlock, outBlock);
    dmSetArray_f( outDes, GlobalOutData, npix );
    dmImageClose( outBlock );
  } else {
    return(-1);
  }

  if ( (strlen(areafile)>0) && (ds_strcmp_cis(areafile,"none")!=0) ) {
    if ( ds_clobber( areafile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(areafile, dmFLOAT, lAxes, 2 );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", areafile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmradar", 0 );
      put_param_hist_info( outBlock, "dmradar", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_f( outDes, GlobalOutArea, npix );
      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }

  if ( (strlen(snrfile)>0) && (ds_strcmp_cis(snrfile,"none")!=0) ) {
    if ( ds_clobber( snrfile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(snrfile, dmFLOAT, lAxes, 2 );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", snrfile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmradar", 0 );
      put_param_hist_info( outBlock, "dmradar", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_f( outDes, GlobalOutSNR, npix );
      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }


  if ( (strlen(maskfile)>0) && (ds_strcmp_cis(maskfile,"none")!=0) ) {
    if ( ds_clobber( maskfile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(maskfile, dmULONG, lAxes, 2 );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", maskfile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmradar", 0 );
      put_param_hist_info( outBlock, "dmradar", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_ul( outDes, GlobalMask, npix );

      dmBlockClose( dmTableWriteRegion( dmBlockGetDataset( outBlock ),
              "REGION", NULL, maskRegion ));

      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }
  
  /* Must keep open until now to do all the wcs/hdr copies */
  dmImageClose( inBlock ); 
  
  /* make valgrind happy */
  free(GlobalData);
  free(GlobalDErr);
  free(GlobalOutData);
  free(GlobalOutArea);
  free(GlobalOutSNR);
  free(GlobalMask);


  return(0);

}
