/******************************************************************************
 *
 * File:           stringtable.c
 *
 * Created:        27/09/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        String table -- implementation
 *
 * Description:    An expandable array of strings.
 *
 * Revisions:      12/2012 PS: Changed prefixes "stringtable_" to "st_"
 *
 *****************************************************************************/

#define STRINGTABLE_NINC 50

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include "stringtable.h"

extern int verbose;

typedef struct {
    char* s;
    int index;                  /* An index. Is set during the string entry.
                                 * If "-1" entered, is set to index of this
                                 * entry in the stringtable. */
    int naccess;                /* Number of times looked for. */
} stringentry;

struct stringtable {
    char* name;
    int unique;                 /* flag: whether all entries must be unique;
                                 * 1 by default */
    int n;
    int sorted;                 /* flag */
    stringentry** se;
};

static void quit_def(char* format, ...);
static st_quit_fn quit = quit_def;

static void quit_def(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  error: stringtable: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    exit(1);
}

/** Compare procedure for qsort(). Non-case-sensitive.
 * @param se1 Pointer to stringentry
 * @param se2 Pointer to stringentry
 * @return 0 if se1 and se2 are identical, -1 or 1 otherwise.
 */
static int cmpse(const void* se1, const void* se2)
{
    return strcasecmp((*(stringentry **) se1)->s, (*(stringentry **) se2)->s);
}

/** Stringtable constructor.
 * @param name Stringtable name
 * @return Stringtable
 */
stringtable* st_create(char* name)
{
    stringtable* st = malloc(sizeof(stringtable));

    if (name != NULL)
        st->name = strdup(name);
    else
        st->name = strdup("stringtable");
    st->unique = 1;
    st->n = 0;
    st->sorted = 0;
    st->se = NULL;

    return st;
}

/** Stringtable destructor.
 * @param st Stringtable
 */
void st_destroy(stringtable* st)
{
    int i;

    if (st == NULL)
        return;

    for (i = 0; i < st->n; ++i) {
        free(st->se[i]->s);
        free(st->se[i]);
    }
    free(st->name);
    free(st->se);
    free(st);
}

/** Stringtable constructor.
 * @param parent Parent stringtable
 * @return Stringtable
 */
stringtable* st_copy(stringtable* parent)
{
    stringtable* st;
    int i;

    if (parent == NULL)
        return NULL;

    st = malloc(sizeof(stringtable));
    st->name = strdup(parent->name);
    st->unique = parent->unique;
    st->n = parent->n;
    st->sorted = parent->sorted;
    st->se = malloc((st->n / STRINGTABLE_NINC + 1) * STRINGTABLE_NINC * sizeof(void*));
    for (i = 0; i < st->n; ++i) {
        stringentry* se0 = parent->se[i];
        stringentry* se = malloc(sizeof(stringentry));

        se->s = strdup(se0->s);
        se->index = se0->index;
        se->naccess = 0;

        st->se[i] = se;
    }

    return st;
}

/** Returns the size of the table.
 * @param st Stringtable
 */
int st_getsize(stringtable* st)
{
    return st->n;
}

/** Resets contents of a stringtable.
 * @param st Stringtable
 */
void st_reset(stringtable* st)
{
    st->n = 0;
}

/** Adds an entry to a stringtable.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 * @return String index
 */
int st_add(stringtable* st, char* s, int index)
{
    stringentry* se = NULL;

    if (st->unique) {
        int previndex = st_findindexbystring(st, s);

        if (previndex >= 0) {
            if (previndex != index)
                quit("\"%s\": entry \"%s\" duplicated\n", st->name, s);
            else
                return index;
        }
    }

    if (st->n % STRINGTABLE_NINC == 0)
        st->se = realloc(st->se, (st->n + STRINGTABLE_NINC) * sizeof(void*));

    st->se[st->n] = malloc(sizeof(stringentry));
    se = st->se[st->n];
    se->s = strdup(s);
    se->index = (index >= 0) ? index : st->n;
    se->naccess = 0;
    st->n++;
    st->sorted = 0;

    return se->index;
}

/** Adds an entry to a stringtable if it is not already there.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 * @return String index
 */
int st_add_ifabsent(stringtable* st, char* s, int index)
{
    stringentry* se = NULL;
    int previndex = st_findindexbystring(st, s);

    if (previndex >= 0)
        return previndex;

    if (st->n % STRINGTABLE_NINC == 0)
        st->se = realloc(st->se, (st->n + STRINGTABLE_NINC) * sizeof(void*));

    st->se[st->n] = malloc(sizeof(stringentry));
    se = st->se[st->n];
    se->s = strdup(s);
    se->index = (index >= 0) ? index : st->n;
    se->naccess = 0;
    st->n++;
    st->sorted = 0;

    return se->index;
}

/** Finds index associated with a string in a stringtable. Uses 
 * binary search for a sorted table; cycles trhough all entries otherwise.
 * @param st Stringtable
 * @param s String
 * @return Index of the string if found; -1 otherwise
 */
int st_findindexbystring(stringtable* st, char* s)
{
    int index = -1;

    if (st->sorted) {
        stringentry tmp;

        tmp.s = s;
        tmp.index = -1;

        {
            stringentry* se = &tmp;
            stringentry** ans = bsearch(&se, st->se, st->n, sizeof(void*), cmpse);

            if (ans != NULL) {
                index = (*ans)->index;
                (*ans)->naccess++;
            }
        }
    } else {
        int n = st->n;
        int i;

        for (i = 0; i < n; ++i)
            if (strcasecmp(st->se[i]->s, s) == 0) {
                index = st->se[i]->index;
                st->se[i]->naccess++;
                break;
            }
    }

    return index;
}

/** Finds string associated with an index in a stringtable. Uses 
 * linear search.
 * @param st Stringtable
 * @param index Index
 * @return String if found; NULL otherwise
 */
char* st_findstringbyindex(stringtable* st, int index)
{
    int i;

    for (i = 0; i < st->n; ++i)
        if (st->se[i]->index == index)
            return st->se[i]->s;

    return NULL;
}

/** Sorts a stringtable using qsort().
 * @param st pointer to stringtable
 */
void st_sort(stringtable* st)
{
    if (st->sorted)
        return;
    qsort(st->se, st->n, sizeof(void*), cmpse);
    st->sorted = 1;
}

/** Prints contents of a stringtable to standard error.
 * @param st Stringtable
 */
void st_print(stringtable* st)
{
    int i;

    printf("%s:\n", st->name);
    for (i = 0; i < st->n; ++i)
        printf("  %s(%d)\n", st->se[i]->s, st->se[i]->index);
    fflush(stdout);
}

/** Prints stringtable entries with specified separator to standard output.
 * @param st Stringtable
 * @param sep Separator
 */
void st_printentries(stringtable* st, char* sep)
{
    int n = st->n;
    int i;

    assert(sep != NULL);

    if (n < 1)
        return;

    printf("%s", st->se[0]->s);
    for (i = 1; i < n; ++i)
        printf("%s%s", sep, st->se[i]->s);
    fflush(stdout);
}

/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 */
void st_printentry(stringtable* st, int i)
{
    fprintf(stdout, "  %s: %s(%d)\n", st->name, st->se[i]->s, st->se[i]->index);
    fflush(stdout);
}
