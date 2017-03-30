/******************************************************************************
 *
 * File:        allhs.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_ALLHS_H)

typedef void (*H_fn) (dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, void* psrc, ENSOBSTYPE dst[]);

H_fn getH(int issurface, char H_tag[]);

#define _ALLHS_H
#endif
