/******************************************************************************
 *
 * File:           stringtable.h
 *
 * Created:        27/09/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        String table -- header
 *
 * Description:    Stringtable is basically an expandable array of strings
 *                 and associated indices.
 *
 * Revisions:      12/2012 PS: Changed prefixes "stringtable_" to "st_"
 *
 *****************************************************************************/

#if !defined(_STRINGTABLE_H)

typedef struct {
    char* s;
    int index;                  /* An index. Is set during the string entry.
                                 * If "-1" entered, is set to index of this
                                 * entry in the stringtable. */
    int naccess;                /* Number of times looked for. */
} stringentry;

#if !defined(STRINGTABLE_STRUCT)
struct stringtable;
typedef struct stringtable stringtable;

#define STRINGTABLE_STRUCT
#endif

struct stringtable {
    char* name;
    int unique;                 /* flag: whether all entries must be unique;
                                 * 1 by default */
    int n;
    int sorted;                 /* flag */
    stringentry** se;
};

/** Stringtable constructor.
 * @param name Stringtable name
 * @return Stringtable
 */
stringtable* st_create(char* name);

/** Stringtable constructor.
 * @param parent Parent stringtable
 * @return Stringtable
 */
stringtable* st_copy(stringtable* parent);

/** Stringtable destructor.
 * @param st Stringtable
 */
void st_destroy(stringtable* st);

/** Adds an entry to a stringtable.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 * @return String index
 */
int st_add(stringtable* st, char* s, int index);

/** Adds an entry to a stringtable if it is not already there.
 * @param st Stringtable
 * @param s String to be added
 * @param len Length of the string
 * @param index External index of the added string to be stored
 * @return String index
 */
int st_add_ifabsent(stringtable* st, char* s, int index);

/** Finds index associated with a string in a stringtable. Uses 
 * binary search for a sorted table; cycles trhough all entries otherwise.
 * @param st Stringtable
 * @param s String
 * @return Index of the string if found; -1 otherwise
 */
int st_findindexbystring(stringtable* st, char* s);

/** Finds string associated with an index in a stringtable. Uses 
 * linear search.
 * @param st Stringtable
 * @param i Index
 * @return String if found; NULL otherwise
 */
char* st_findstringbyindex(stringtable* st, int i);

/** Resets contents of a stringtable.
 * @param st Stringtable
 */
void st_reset(stringtable* st);

/** Sorts a stringtable using qsort().
 * @param st Stringtable
 */
void st_sort(stringtable* st);

/** Prints contents of a stringtable to standard error.
 * @param st Stringtable
 */
void st_print(stringtable* st);

/** Prints stringtable entries with specified separator to standard error.
 * @param st Stringtable
 * @param sep Separator
 */
void st_printentries(stringtable* st, char* sep);

/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 */
void st_printentry(stringtable* st, int i);

#define _STRINGTABLE_H
#endif
