/******************************************************************************
 *
 * File:           ncw.h
 *  
 * Created         19/10/2000
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *                 Bureau of Meteorology
 *  
 * Purpose:        Simple wrappers to netcdf library procedures for
 *                 better error messaging.
 *                 Some straightforward extensions to netcdf library 
 *                 procedures.
 *
 * Revisions:      PS 19/08/2003 Changed library name from "nc" to "ncw". 
 *                    Changed procedure prefixes from "NC_" to "ncw_".
 *                 PS 24/08/2016 Version 2.00: eliminated file names from the
 *                    functions' arguments, using nc_inq_path() instead. Most
 *                    of the wrappers have now the same format as their siblings
 *                    in the NetCDF library.
 *
 *****************************************************************************/

#if !defined(_NCW_H)

#include <netcdf.h>

extern const char ncw_version[];
extern int ncw_chunkbylayers;

/* It is possible to set the quit procedure. By default, the internal procedure
 * is used.
 */
typedef void (*ncw_quit_fn) (char* format, ...);
void ncw_set_quitfn(ncw_quit_fn quit_fn);

/* These procedures are the straightforward wrappers of the corresponding
 * procedures in netcdf library.
 */
char* ncw_get_path(int ncid);
void ncw_create(const char fname[], int mode, int* ncid);
void ncw_open(const char fname[], int mode, int* ncid);
void ncw_redef(int ncid);
void ncw_enddef(int ncid);
void ncw_sync(int ncid);
void ncw_close(int ncid);
void ncw_inq(int ncid, int* ndims, int* nvars, int* natts, int* unlimdimid);
void ncw_inq_ndims(int ncid, int* ndims);
void ncw_inq_nvars(int ncid, int* nvars);
void ncw_inq_natts(int ncid, int* natts);
void ncw_inq_unlimdim(int ncid, int* unlimdimid);
void ncw_inq_format(int ncid, int* format);
void ncw_def_dim(int ncid, const char dimname[], size_t len, int* dimid);
void ncw_inq_dimid(int ncid, const char dimname[], int* dimid);
void ncw_inq_dim(int ncid, int dimid, char dimname[], size_t* len);
void ncw_inq_dimname(int ncid, int dimid, char dimname[]);
void ncw_inq_dimlen(int ncid, int dimid, size_t* len);
void ncw_rename_dim(int ncid, const char oldname[], const char newname[]);
void ncw_def_var(int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid);
void ncw_def_var_deflate(int ncid, int varid, int shuffle, int deflate, int deflate_level);
void ncw_def_var_fill(int ncid, int varid, int nofill, void* fillvalue);
void ncw_inq_varid(int ncid, const char varname[], int* varid);
void ncw_inq_var(int ncid, int varid, char varname[], nc_type* xtype, int* ndims, int dimids[], int* natts);
void ncw_inq_varname(int ncid, int varid, char varname[]);
void ncw_inq_vartype(int ncid, int varid, nc_type* xtype);
void ncw_inq_varndims(int ncid, int varid, int* ndims);
void ncw_inq_vardimid(int ncid, int varid, int dimids[]);
void ncw_inq_varnatts(int ncid, int varid, int* natts);
void ncw_inq_varsize(int ncid, int varid, size_t* size);
void ncw_inq_var_deflate(int ncid, int varid, int* shuffle, int* deflate, int* deflate_level);
void ncw_inq_var_fill(int ncid, int varid, int* nofill, void* fillvalue);
void ncw_rename_var(int ncid, const char oldname[], const char newname[]);
void ncw_put_var(int ncid, int varid, const void* v);
void ncw_put_var_text(int ncid, int varid, const char v[]);
void ncw_put_var_schar(int ncid, int varid, const signed char v[]);
void ncw_put_var_uchar(int ncid, int varid, const unsigned char v[]);
void ncw_put_var_short(int ncid, int varid, const short int v[]);
void ncw_put_var_ushort(int ncid, int varid, const unsigned short int v[]);
void ncw_put_var_int(int ncid, int varid, const int v[]);
void ncw_put_var_uint(int ncid, int varid, const unsigned int v[]);
void ncw_put_var_long(int ncid, int varid, const long v[]);
void ncw_put_var_float(int ncid, int varid, const float v[]);
void ncw_put_var_double(int ncid, int varid, const double v[]);
void ncw_get_var(int ncid, int varid, void* v);
void ncw_get_var_text(int ncid, int varid, char v[]);
void ncw_get_var_schar(int ncid, int varid, signed char v[]);
void ncw_get_var_uchar(int ncid, int varid, unsigned char v[]);
void ncw_get_var_short(int ncid, int varid, short int v[]);
void ncw_get_var_ushort(int ncid, int varid, unsigned short int v[]);
void ncw_get_var_int(int ncid, int varid, int v[]);
void ncw_get_var_uint(int ncid, int varid, unsigned int v[]);
void ncw_get_var_float(int ncid, int varid, float v[]);
void ncw_get_var_float_fixerange(int ncid, int varid, float v[]);
void ncw_get_var_double(int ncid, int varid, double v[]);
void ncw_get_var1_double(int ncid, int varid, const size_t len[], double* in);
void ncw_put_vara(int ncid, int varid, const size_t start[], const size_t count[], void* v);
void ncw_put_vara_text(int ncid, int varid, const size_t start[], const size_t count[], const char v[]);
void ncw_put_vara_short(int ncid, int varid, const size_t start[], const size_t count[], const short int v[]);
void ncw_put_vara_ushort(int ncid, int varid, const size_t start[], const size_t count[], const unsigned short int v[]);
void ncw_put_vara_int(int ncid, int varid, const size_t start[], const size_t count[], const int v[]);
void ncw_put_vara_float(int ncid, int varid, const size_t start[], const size_t count[], const float v[]);
void ncw_get_vara_float_fixerange(int ncid, int varid, const size_t start[], const size_t count[], float v[]);
void ncw_put_vara_double(int ncid, int varid, const size_t start[], const size_t count[], const double v[]);
void ncw_get_vara(int ncid, int varid, const size_t start[], const size_t count[], void* v);
void ncw_get_vara_text(int ncid, int varid, const size_t start[], const size_t count[], char v[]);
void ncw_get_vara_short(int ncid, int varid, const size_t start[], const size_t count[], short int v[]);
void ncw_get_vara_int(int ncid, int varid, const size_t start[], const size_t count[], int v[]);
void ncw_get_vara_float(int ncid, int varid, const size_t start[], const size_t count[], float v[]);
void ncw_get_vara_double(int ncid, int varid, const size_t start[], const size_t count[], double v[]);
void ncw_put_att_text(int ncid, int varid, const char attname[], const char v[]);
void ncw_put_att_uchar(int ncid, int varid, const char attname[], size_t len, const unsigned char v[]);
void ncw_put_att_short(int ncid, int varid, const char attname[], size_t len, const short int v[]);
void ncw_put_att_ushort(int ncid, int varid, const char attname[], size_t len, const unsigned short int v[]);
void ncw_put_att_int(int ncid, int varid, const char attname[], size_t len, const int v[]);
void ncw_put_att_long(int ncid, int varid, const char attname[], size_t len, const long int v[]);
void ncw_put_att_float(int ncid, int varid, const char attname[], size_t len, const float v[]);
void ncw_put_att_double(int ncid, int varid, const char attname[], size_t len, const double v[]);
void ncw_inq_attname(int ncid, int varid, int attrid, char attname[]);
void ncw_inq_att(int ncid, int varid, const char attname[], nc_type* xtype, size_t* len);
void ncw_inq_atttype(int ncid, int varid, const char attname[], nc_type* xtype);
void ncw_inq_attlen(int ncid, int varid, const char attname[], size_t* len);
void ncw_copy_att(int ncid_src, int varid_src, const char attname[], int ncid_dst, int varid_dst);
void ncw_rename_att(int ncid, int varid, const char oldname[], const char newname[]);
void ncw_del_att(int ncid, int varid, const char name[]);
void ncw_get_att(int ncid, int varid, const char attname[], void* v);
void ncw_get_att_text(int ncid, int varid, const char attname[], char v[]);
void ncw_get_att_schar(int ncid, int varid, const char attname[], signed char v[]);
void ncw_get_att_short(int ncid, int varid, const char attname[], short int v[]);
void ncw_get_att_int(int ncid, int varid, const char attname[], int v[]);
void ncw_get_att_uint(int ncid, int varid, const char attname[], unsigned int v[]);
void ncw_get_att_float(int ncid, int varid, const char attname[], float v[]);
void ncw_get_att_double(int ncid, int varid, const char attname[], double v[]);

/* These procedures do not have direct analogues in the netcdf library.
 */
const char* ncw_nctype2str(nc_type type);
size_t ncw_sizeof(nc_type type);

int ncw_inq_nrecords(int ncid);
void ncw_inq_vardims(int ncid, int varid, int maxndims, int* ndims, size_t dimlen[]);

void ncw_copy_dims(int ncid_src, int ncid_dst);
void ncw_copy_dim(int ncid_src, const char dimname[], int ncid_dst);
int ncw_copy_vardef(int ncid_src, int varid_src, int ncid_dst);
void ncw_copy_atts(int ncid_src, int varid_src, int ncid_dst, int varid_dst);
void ncw_copy_vardata(int ncid_src, int varid_src, int ncid_dst);
void ncw_copy_var(int ncid_src, const char varname[], int ncid_dst);
void ncw_inq_dimid2(int ncid, const char dimname1[], const char dimname2[], int* dimid);
void ncw_get_att_int2(int ncid, int varid, const char attname1[], const char attname2[], int v[]);

void ncw_def_deflate(int ncid, int shuffle, int deflate, int deflate_level);

void ncw_find_timevarid(int ncid, int* varid);
void ncw_find_vars(int ncid, int ndims, const int dims[], const char attr[], const void* attval, int* nvars, int** vids);

int ncw_file_opens(const char fname[], int mode);
int ncw_att_exists(int ncid, int varid, const char attname[]);
int ncw_var_exists(int ncid, const char varname[]);
int ncw_dim_exists(int ncid, const char dimname[]);
int ncw_att_exists2(int ncid, int varid, const char attname[]);

void ncw_def_var_as(int ncid, const char oldvarname[], const char newvarname[]);
size_t ncw_get_varsize(int ncid, int varid);
void ncw_get_var_double_record(int ncid, int varid, int r, double v[]);
void ncw_get_var_float_record(int ncid, int varid, int r, float v[]);
void ncw_get_var_int_record(int ncid, int varid, int r, int v[]);
void ncw_put_var_double_record(int ncid, int varid, int r, double v[]);
void ncw_put_var_float_record(int ncid, int varid, int r, float v[]);
void ncw_put_var_int_record(int ncid, int varid, int r, int v[]);

void ncw_check_attlen(int ncid, int varid, const char attname[], size_t len);
void ncw_check_dimlen(int ncid, const char dimname[], size_t len);
void ncw_check_varndims(int ncid, int varid, int ndims);
void ncw_check_vardims(int ncid, int varid, int ndims, size_t dimlen[]);
void ncw_check_varsize(int ncid, int varid, size_t size);

int ncw_var_hasunlimdim(int ncid, int varid);

#define _NCW_H
#endif
