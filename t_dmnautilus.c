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

/***********************************************************************/
/*     Copyrights:                                                      */
/*                                                                      */
/*     Copyright (c) 2001 Smithsonian Astrophysical Observatory         */
/*                                                                      */
/*     Permission to use, copy, modify, distribute, and  sell  this     */
/*     software  and  its  documentation  for any purpose is hereby     */
/*     granted without  fee,  provided  that  the  above  copyright     */
/*     notice  appear  in  all  copies and that both that copyright     */
/*     notice and this permission notice appear in supporting docu-     */
/*     mentation,  and  that  the  name  of the  Smithsonian Astro-     */
/*     physical Observatory not be used in advertising or publicity     */
/*     pertaining to distribution of the software without specific,     */
/*     written  prior  permission.   The Smithsonian  Astrophysical     */
/*     Observatory makes no representations about  the  suitability     */
/*     of  this  software for any purpose.  It is provided  "as is"     */
/*     without express or implied warranty.                             */
/*     THE  SMITHSONIAN  ASTROPHYSICAL  OBSERVATORY  DISCLAIMS  ALL     */
/*     WARRANTIES  WITH  REGARD  TO  THIS  SOFTWARE,  INCLUDING ALL     */
/*     IMPLIED  WARRANTIES  OF  MERCHANTABILITY AND FITNESS, IN  NO     */
/*     EVENT SHALL THE  SMITHSONIAN  ASTROPHYSICAL  OBSERVATORY  BE     */
/*     LIABLE FOR  ANY SPECIAL, INDIRECT  OR  CONSEQUENTIAL DAMAGES     */
/*     OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS  OF USE,  DATA     */
/*     OR  PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR     */
/*     OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION  WITH     */
/*     THE  USE OR PERFORMANCE OF THIS SOFTWARE.                        */
/*                                                                      */
/************************************************************************/
 
/* H*****************************************************************
 *
 * FILE NAME: t_abin.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * IRAF wrapper for abin tool
 *  
 *
 *
H***************************************************************** */

#include "dslib.h"


void main(int argc, char** argv)
{
    int fail_status = 0; 

    dsErrInitLib(dsPTGRPERR, "dmradar");

    /* INIT THE IRAF ENVIRONMENT */
    /* IRAFHOST_INIT;		*/

    /* OPEN THE PARAMETER FILE */
    if(clinit(argv, argc, "rw") == NULL)
    {
       err_msg( "Problem opening parameter file %s.par\n", argv[0]);
       err_msg( "Parameter library error: %s.\n", paramerrstr());
       fail_status = -1;
    }
    else
    {    
       /* EXECUTE OUR PROGRAM */ 
       fail_status = abin();
    
       /* CLOSE PARAMETER FILE AND RETURN TO THE OS */
       clclose();
    } 
    
    exit(fail_status); 
}
