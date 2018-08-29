/******************************************************************************
 *
 * File:           stringtable.h
 *
 * Created:        27/09/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        String table
 *
 * Description:    An expandable array of strings.
 *
 * Revisions:      12/2012 PS: Changed prefixes "stringtable_" to "st_"
 *
 *****************************************************************************/

#if !defined(_STRINGTABLE_H)

struct stringtable;
typedef struct stringtable stringtable;
typedef void (*st_quit_fn) (char* format, ...);

stringtable* st_create(char* name);
void st_destroy(stringtable* st);
stringtable* st_copy(stringtable* parent);
int st_getsize(stringtable* st);
int st_add(stringtable* st, char* s, int index);
int st_add_ifabsent(stringtable* st, char* s, int index);
int st_findindexbystring(stringtable* st, char* s);
char* st_findstringbyindex(stringtable* st, int i);
void st_reset(stringtable* st);
void st_sort(stringtable* st);
void st_print(stringtable* st);
void st_printentries(stringtable* st, char* sep);
void st_printentry(stringtable* st, int i);
void st_set_quitfn(st_quit_fn quit_fn);

#define _STRINGTABLE_H
#endif
