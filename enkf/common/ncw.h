/******************************************************************************
 *
 * File:           ncw.h
 *  
 * Created         19/10/2000
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *  
 * Purpose:        Simple wrappers to netcdf library procedures for
 *                 better error messaging.
 *                 Some straightforward extensions to netcdf library 
 *                 procedures.
 *
 * Revisions:      PS 19/08/2003 Changed library name from "nc" to "ncw". 
 *                    Changed procedure prefixes from "NC_" to "ncw_".
 *
 *****************************************************************************/

#if !defined(_NCW_H)

#include <netcdf.h>

extern const char ncw_version[];

/* These procedures are straightforward wrappers of the corresponding
 * procedures in netcdf library.
 */
void ncw_create(const char fname[], int mode, int* ncid);
void ncw_open(const char fname[], int mode, int* ncid);
void ncw_redef(const char fname[], int ncid);
void ncw_enddef(const char fname[], int ncid);
void ncw_sync(const char fname[], int ncid);
void ncw_close(const char fname[], int ncid);
void ncw_inq(const char fname[], int ncid, int* ndims, int* nvars, int* natts, int* unlimdimid);
void ncw_inq_ndims(const char fname[], int ncid, int* ndims);
void ncw_inq_nvars(const char fname[], int ncid, int* nvars);
void ncw_inq_natts(const char fname[], int ncid, int* natts);
void ncw_inq_unlimdimid(const char fname[], int ncid, int* unlimdimid);
void ncw_def_dim(const char fname[], int ncid, const char dimname[], size_t len, int* dimid);
void ncw_inq_dimid(const char fname[], int ncid, const char dimname[], int* dimid);
void ncw_inq_dim(const char fname[], int ncid, int dimid, char dimname[], size_t* len);
void ncw_inq_dimname(const char fname[], int ncid, int dimid, char dimname[]);
void ncw_inq_dimlen(const char fname[], int ncid, int dimid, size_t* len);
void ncw_rename_dim(const char fname[], int ncid, const char oldname[], const char newname[]);
void ncw_def_var(const char fname[], int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid);
void ncw_inq_varid(const char fname[], int ncid, const char varname[], int* varid);
void ncw_inq_var(const char fname[], int ncid, int varid, char varname[], nc_type* xtype, int* ndims, int dimids[], int* natts);
void ncw_inq_varname(const char fname[], int ncid, int varid, char varname[]);
void ncw_inq_vartype(const char fname[], int ncid, int varid, nc_type* xtype);
void ncw_inq_varndims(const char fname[], int ncid, int varid, int* ndims);
void ncw_inq_vardimid(const char fname[], int ncid, int varid, int dimids[]);
void ncw_inq_varnatts(const char fname[], int ncid, int varid, int* natts);
void ncw_rename_var(const char fname[], int ncid, const char oldname[], const char newname[]);
void ncw_put_var_text(const char fname[], int ncid, int varid, const char v[]);
void ncw_put_var_short(const char fname[], int ncid, int varid, const short int v[]);
void ncw_put_var_int(const char fname[], int ncid, int varid, const int v[]);
void ncw_put_var_float(const char fname[], int ncid, int varid, const float v[]);
void ncw_put_var_double(const char fname[], int ncid, int varid, const double v[]);
void ncw_get_var_text(const char fname[], int ncid, int varid, char v[]);
void ncw_get_var_short(const char fname[], int ncid, int varid, short int v[]);
void ncw_get_var_int(const char fname[], int ncid, int varid, int v[]);
void ncw_get_var_float(const char fname[], int ncid, int varid, float v[]);
void ncw_get_var_double(const char fname[], int ncid, int varid, double v[]);
void ncw_get_var1_double(const char fname[], int ncid, int varid, const size_t len[], double* in);
void ncw_put_vara_text(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const char v[]);
void ncw_put_vara_short(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const short int v[]);
void ncw_put_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const int v[]);
void ncw_put_vara_float(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const float v[]);
void ncw_put_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const double v[]);
void ncw_get_vara_text(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], char v[]);
void ncw_get_vara_short(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], short int v[]);
void ncw_get_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], int v[]);
void ncw_get_vara_float(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], float v[]);
void ncw_get_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], double v[]);
void ncw_put_att_text(const char fname[], int ncid, int varid, const char attname[], const char v[]);
void ncw_put_att_int(const char fname[], int ncid, int varid, const char attname[], size_t len, const int v[]);
void ncw_put_att_float(const char fname[], int ncid, int varid, const char attname[], size_t len, const float v[]);
void ncw_put_att_double(const char fname[], int ncid, int varid, const char attname[], size_t len, const double v[]);
void ncw_inq_attname(const char fname[], int ncid, int varid, int attrid, char attname[]);
void ncw_inq_att(const char fname[], int ncid, int varid, const char attname[], nc_type* xtype, size_t* len);
void ncw_inq_attlen(const char fname[], int ncid, int varid, const char attname[], size_t* len);
void ncw_copy_att(const char fname_src[], int ncid_src, int varid_src, const char attname[], const char fname_dst[], int ncid_dst, int varid_dst);
void ncw_rename_att(const char fname[], int ncid, const char varname[], const char oldname[], const char newname[]);
void ncw_del_att(const char fname[], int ncid, int varid, const char name[]);
void ncw_get_att_text(const char fname[], int ncid, int varid, const char attname[], char v[]);
void ncw_get_att_int(const char fname[], int ncid, int varid, const char attname[], int v[]);
void ncw_get_att_float(const char fname[], int ncid, int varid, const char attname[], float v[]);
void ncw_get_att_double(const char fname[], int ncid, int varid, const char attname[], double v[]);

/* These procedures do not have direct analogues in the netcdf library.
 */
int ncw_inq_nrecords(const char fname[], int ncid);
const char* ncw_nctype2str(nc_type type);
size_t ncw_sizeof(nc_type type);
void ncw_copy_dims(const char* fname_src, int ncid_src, const char* fname_dst, int ncid_dst);
int ncw_copy_vardef(const char* fname_src, int ncid_src, int varid_src, const char* fname_dst, int ncid_dst);
void ncw_copy_vardata(const char* fname_src, int ncid_src, int varid_src, const char* fname_dst, int ncid_dst);
void ncw_copy_var(const char* fname_src, int ncid_src, const char varname[], const char* fname_dst, int ncid_dst);
void ncw_inq_dimid2(const char fname[], int ncid, const char dimname1[], const char dimname2[], int* dimid);
void ncw_get_att_int2(const char fname[], int ncid, int varid, const char attname1[], const char attname2[], int v[]);
void ncw_find_vars(const char fname[], int ncid, int ndims, const int dims[], const char attr[], const void* attval, int* nvars, int** vids);
void ncw_find_timevarid(const char fname[], int ncid, int* varid);
int ncw_att_exists(int ncid, int varid, const char attname[]);
int ncw_var_exists(int ncid, const char varname[]);
int ncw_dim_exists(int ncid, const char dimname[]);
void ncw_copy_atts(const char* fname_src, int ncid_src, int varid_src, const char* fname_dst, int ncid_dst, int varid_dst);
void ncw_def_var_as(const char fname[], int ncid, const char oldvarname[], const char newvarname[]);
void ncw_get_var_double_record(const char fname[], int ncid, int varid, int r, double v[]);
void ncw_get_var_float_record(const char fname[], int ncid, int varid, int r, float v[]);
void ncw_put_var_double_record(const char fname[], int ncid, int varid, int r, double v[]);
void ncw_put_var_float_record(const char fname[], int ncid, int varid, int r, float v[]);

typedef void (*ncw_quit_fn) (char* format, ...);
void ncw_set_quitfn(ncw_quit_fn quit_fn);

#define _NCW_H
#endif
