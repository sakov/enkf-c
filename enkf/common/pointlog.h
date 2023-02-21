/******************************************************************************
 *
 * File:        pointlog.h       
 *
 * Created:     7/10/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_POINTLOG_H)

struct pointlog {
    int id;
    double lon, lat;
    /*
     * Grid id of the "native grid". If < 0 (default) -- then log all grids.
     */
    int gridid;
    /*
     * grid identifier used to get gridid at a later stage
     */
    char* gridname;
    /*
     * arrays of fractional grid indices [ngrid][3] corresponding to
     * model.grids[]
     */
    double** fij;
};

void plogs_add(int* nplog, pointlog** plogs, double lon, double lat, char* gridname);
void plogs_destroy(int nplog, pointlog plogs[]);

#define _POINTLOG_H
#endif
