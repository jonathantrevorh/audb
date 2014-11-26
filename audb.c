/*
 * modified 2014-11-25 to implement vp tree
 *
 */

#include "sqlite3ext.h"
SQLITE_EXTENSION_INIT1
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>

#ifndef SQLITE_OMIT_VIRTUALTABLE

long MAX_FEATURE_DIST = (~(unsigned long) 0) >> 1;

/*
** Forward declaration of objects used by this implementation
*/
typedef struct audb_vtab audb_vtab;
typedef struct audb_cursor audb_cursor;
typedef struct vp_node vp_node;
typedef struct feature_t feature_t;

struct feature_t {
    int* feature_start;
    int key_id;
};

struct vp_node {
    sqlite3_int64 id;
    feature_t vp;
    long radius;
    char* path;
    vp_node* left;
    vp_node* right;
};

static vp_node* vpNodeNext(vp_node* node) {
    return NULL;
}

long feature_distance(int* f1, int* f2) {
    long dist = 0;
    int i;
    int len = sizeof(f2)/sizeof(int);
    for (i = 0 ; i < len ; i++) {
        unsigned long t = f1[i] ^ f2[i];
        while (t != 0) {
            dist += t & 1;
            t >>= 1;
        }
    }
    return dist;
}


static void internal_search(vp_node* root, int* sample, long* feature_dist, feature_t* match, int isRightCheck) {
    if (root == NULL) {
        *feature_dist = MAX_FEATURE_DIST;
        return;
    }

    long cDist = feature_distance(sample, root->vp.feature_start);
    if (isRightCheck && cDist >= root->radius) {
        return;
    }

    if (cDist < root->radius) {
        // within bounding radius, so search left
        long tDist = cDist;
        feature_t found;
        internal_search(root->left, sample, &tDist, &found, 0);

        if (tDist < root->radius - cDist) {
            if (tDist < cDist) {
                *feature_dist = tDist;
                *match = found;
            } else {
                *feature_dist = cDist;
                *match = root->vp;
            }
        } else {
            long rDist = tDist;
            feature_t rFound;
            internal_search(root->right, sample, &rDist, &rFound, 1);
            if (rDist < tDist) {
                *feature_dist = rDist;
                *match = rFound;
            } else {
                if (tDist < cDist) {
                    *feature_dist = tDist;
                    *match = found;
                } else {
                    *feature_dist = cDist;
                    *match = root->vp;
                }
            }
        }
    } else {
        long tDist = cDist;
        feature_t found;
        internal_search(root->right, sample, &tDist, &found, 0);
        if (tDist < cDist) {
            *feature_dist = tDist;
            *match = found;
        } else {
            *feature_dist = cDist;
            *match = root->vp;
        }
    }
}

/* Insert a new node pNew. Return NULL on success. If the key is not
** unique, then do not perform the insert but instead leave pNew unchanged
** and return a pointer to an existing node with the same key.
*/
static void vpNodeInsert(vp_node *pHead, vp_node *pNew) {
  vp_node *p = ppHead;
  // TODO: insert just this node
}

/* Walk the tree can call xDestroy on each node
*/
static void audbDestroy(vp_node *p, void (*xDestroy)(vp_node*)){
  if( p ){
    audbDestroy(p->left, xDestroy);
    audbDestroy(p->right, xDestroy);
    xDestroy(p);
  }
}
/*
** End of the AVL Tree implementation
******************************************************************************/

/*
** A closure virtual-table object
*/
struct audb_vtab {
  sqlite3_vtab base;         /* Base class - must be first */
  char *zDb;                 /* Name of database.  (ex: "main") */
  char *zSelf;               /* Name of this virtual table */
  sqlite3 *db;               /* The database connection */
  int nCursor;               /* Number of pending cursors */
  vp_node *root;
};

/* A closure cursor object */
struct audb_cursor {
  sqlite3_vtab_cursor base;  /* Base class - must be first */
  audb_vtab *pVtab;       /* The virtual table this cursor belongs to */
  vp_node *pCurrent;         /* Current element of output */
};

/*
** This function converts an SQL quoted string into an unquoted string
** and returns a pointer to a buffer allocated using sqlite3_malloc()
** containing the result. The caller should eventually free this buffer
** using sqlite3_free.
**
** Examples:
**
**     "abc"   becomes   abc
**     'xyz'   becomes   xyz
**     [pqr]   becomes   pqr
**     `mno`   becomes   mno
*/
static char *audbDequote(const char *zIn){
  int nIn;                        /* Size of input string, in bytes */
  char *zOut;                     /* Output (dequoted) string */

  nIn = (int)strlen(zIn);
  zOut = sqlite3_malloc(nIn+1);
  if( zOut ){
    char q = zIn[0];              /* Quote character (if any ) */

    if( q!='[' && q!= '\'' && q!='"' && q!='`' ){
      memcpy(zOut, zIn, nIn+1);
    }else{
      int iOut = 0;               /* Index of next byte to write to output */
      int iIn;                    /* Index of next byte to read from input */

      if( q=='[' ) q = ']';
      for(iIn=1; iIn<nIn; iIn++){
        if( zIn[iIn]==q ) iIn++;
        zOut[iOut++] = zIn[iIn];
      }
    }
    assert( (int)strlen(zOut)<=nIn );
  }
  return zOut;
}

/*
** Deallocate an audb_vtab object
*/
static void audbFree(audb_vtab *p){
  if( p ){
    sqlite3_free(p->zDb);
    sqlite3_free(p->zSelf);
    memset(p, 0, sizeof(*p));
    sqlite3_free(p);
  }
}

/*
** xDisconnect/xDestroy method for the audb module.
*/
static int audbDisconnect(sqlite3_vtab *pVtab){
  audb_vtab *p = (audb_vtab*)pVtab;
  assert( p->nCursor==0 );
  audbFree(p);
  return SQLITE_OK;
}

/*
** Check to see if the argument is of the form:
**
**       KEY = VALUE
**
** If it is, return a pointer to the first character of VALUE.
** If not, return NULL.  Spaces around the = are ignored.
*/
static const char *audbValueOfKey(const char *zKey, const char *zStr){
  int nKey = (int)strlen(zKey);
  int nStr = (int)strlen(zStr);
  int i;
  if( nStr<nKey+1 ) return 0;
  if( memcmp(zStr, zKey, nKey)!=0 ) return 0;
  for(i=nKey; isspace(zStr[i]); i++){}
  if( zStr[i]!='=' ) return 0;
  i++;
  while( isspace(zStr[i]) ){ i++; }
  return zStr+i;
}

/*
** xConnect/xCreate method for the audb module. Arguments are:
**
**   argv[0]    -> module name  ("audb")
**   argv[1]    -> database name
**   argv[2]    -> table name
**   argv[3...] -> arguments
*/
static int audbConnect(sqlite3 *db, void *pAux, int argc, const char *const*argv, sqlite3_vtab **ppVtab, char **pzErr) {
  int rc = SQLITE_OK;              /* Return code */
  audb_vtab *pNew = 0;          /* New virtual table */
  const char *zDb = argv[1];
  const char *zVal;
  int i;

  (void)pAux;
  *ppVtab = 0;
  pNew = sqlite3_malloc( sizeof(*pNew) );
  if( pNew==0 ) return SQLITE_NOMEM;
  rc = SQLITE_NOMEM;
  memset(pNew, 0, sizeof(*pNew));
  pNew->db = db;
  pNew->zDb = sqlite3_mprintf("%s", zDb);
  if( pNew->zDb==0 ) goto audbConnectError;
  pNew->zSelf = sqlite3_mprintf("%s", argv[2]);
  if( pNew->zSelf==0 ) goto audbConnectError;
  for(i=3; i<argc; i++){
    continue;
    *pzErr = sqlite3_mprintf("unrecognized argument: [%s]\n", argv[i]);
    audbFree(pNew);
    *ppVtab = 0;
    return SQLITE_ERROR;
  }
  rc = sqlite3_declare_vtab(db,
         "CREATE TABLE x(id,path,root HIDDEN,tablename HIDDEN,"
                        "idcolumn HIDDEN,parentcolumn HIDDEN)"
       );
#define AUDB_COL_ID              0
#define AUDB_COL_PATH            1
#define AUDB_COL_ROOT            2
#define AUDB_COL_TABLENAME       3
#define AUDB_COL_IDCOLUMN        4
#define AUDB_COL_PARENTCOLUMN    5
  if( rc!=SQLITE_OK ){
    audbFree(pNew);
  }
  *ppVtab = &pNew->base;
  return rc;

audbConnectError:
  audbFree(pNew);
  return rc;
}

/*
** Open a new audb cursor.
*/
static int audbOpen(sqlite3_vtab *pVTab, sqlite3_vtab_cursor **ppCursor){
  printf("audbopen\n");
  audb_vtab *p = (audb_vtab*)pVTab;
  audb_cursor *pCur;
  pCur = sqlite3_malloc( sizeof(*pCur) );
  if( pCur==0 ) return SQLITE_NOMEM;
  memset(pCur, 0, sizeof(*pCur));
  pCur->pVtab = p;
  *ppCursor = &pCur->base;
  p->nCursor++;
  return SQLITE_OK;
}

/*
** Free up all the memory allocated by a cursor.  Set it rLimit to 0
** to indicate that it is at EOF.
*/
static void audbClearCursor(audb_cursor *pCur){
  audbDestroy(pCur->pVtab->root, (void(*)(vp_node*))sqlite3_free);
  pCur->pCurrent = 0;
  pCur->pVtab->root = 0;
}

/*
** Close a closure cursor.
*/
static int audbClose(sqlite3_vtab_cursor *cur){
  audb_cursor *pCur = (audb_cursor *)cur;
  audbClearCursor(pCur);
  pCur->pVtab->nCursor--;
  sqlite3_free(pCur);
  return SQLITE_OK;
}

/*
** Advance a cursor to its next row of output
*/
static int audbNext(sqlite3_vtab_cursor *cur){
  audb_cursor *pCur = (audb_cursor*)cur;
  pCur->pCurrent = vpNodeNext(pCur->pCurrent);
  return SQLITE_OK;
}

/*
** Allocate a node
*/
static vp_node *vpNodeCreate(
  sqlite3_int64 id,       /* The node ID */
  feature_t *vp        /* The generation number for this node */
){
  vp_node *node = sqlite3_malloc( sizeof(*node) );
  if( node==0 ) return;
  memset(node, 0, sizeof(*node));
  node->id = id;
  node->left = NULL;
  node->right = NULL;
  node->vp = *vp;
  node->radius = 0;
  return node;
}

/*
 * start the search
*/
static int audbFilter(sqlite3_vtab_cursor *pVtabCursor, int idxNum, const char *idxStr, int argc, sqlite3_value **argv) {
    audb_cursor *pCur = (audb_cursor *)pVtabCursor;
    audb_vtab *pVtab = pCur->pVtab;
    (void)idxStr;  /* Unused parameter */
    (void)argc;    /* Unused parameter */
    audbClearCursor(pCur);
    pCur->pCurrent = pCur->pVtab->root;
    int rc = SQLITE_OK;

    return rc;
}

/*
 * given a cursor, and a column index (0 through 5), fill ctz with a sqlite3_result_{type} with value of the column's value for the row at this cursor
 * null return is acceptable
*/
static int audbColumn(sqlite3_vtab_cursor *cur, sqlite3_context *ctx, int i){
  audb_cursor *pCur = (audb_cursor*)cur;
  switch( i ){
    case AUDB_COL_ID: {
      sqlite3_result_int64(ctx, pCur->pCurrent->id);
      break;
    }
    case AUDB_COL_PATH: {
      sqlite3_result_text(ctx, pCur->pCurrent->path, -1, SQLITE_TRANSIENT);
      break;
    }
    case AUDB_COL_ROOT: {
      sqlite3_result_null(ctx);
      break;
    }
  }
  return SQLITE_OK;
}

/*
** The rowid.  For the closure table, this is the same as the "id" column.
*/
static int audbRowid(sqlite3_vtab_cursor *cur, sqlite_int64 *pRowid){
  audb_cursor *pCur = (audb_cursor*)cur;
  *pRowid = pCur->pCurrent->id;
  return SQLITE_OK;
}

static int audbUpdate(sqlite3_vtab *pTab, int argc, sqlite3_value **argv, sqlite_int64 *pRowid) {
    audb_vtab *pVTab = (audb_vtab*)pTab;
    int isDelete = argc == 1;
    if (isDelete) {
        // TODO: handle deletion of things
        const char errorMsg[] = "DELETE command not supported";
        printf(errorMsg);
        pTab->zErrMsg = sqlite3_mprintf(errorMsg);
        return SQLITE_ERROR;
    }

    int zerothArgType = sqlite3_value_type(argv[0]);
    int isInsert = zerothArgType == SQLITE_NULL;
    if (isInsert) {
        const unsigned char* id = sqlite3_value_text(argv[2]);
        const unsigned char* path = sqlite3_value_text(argv[3]);
        printf("insert id: %s, path: %s\n", id, path);
        // TODO: read file, generate feature, create vp_node, insert into tree
        feature_t *feature = sqlite3_malloc( sizeof(*feature) );
        vp_node *node = vpNodeCreate(0, feature);
        vpNodeInsert(pVTab->root, node);
        return SQLITE_OK;
    }

    int updatedRowId = argv[0] == argv[1];
    if (updatedRowId) {
        const char errorMsg[] = "updating row id not supported";
        printf(errorMsg);
        pTab->zErrMsg = sqlite3_mprintf(errorMsg);
        return SQLITE_ERROR;
    }

    // update values
    const char errorMsg[] = "updating row values";
    printf(errorMsg);
    pTab->zErrMsg = sqlite3_mprintf(errorMsg);
    return SQLITE_ERROR;
}


/*
** EOF indicator
*/
static int audbEof(sqlite3_vtab_cursor *cur){
  audb_cursor *pCur = (audb_cursor*)cur;
  return pCur->pCurrent==0;
}

/*
 * tells sqlite how to search things, using heuristics to find the most optimal attributes or conditions to search on
 * we can ignore this for now. the conditions are not acutally used by sqlite, but rather passed into the filter function
 * so if we don't use them in filter, we don't need to make them here
 *
** Search for terms of these forms:
**
**   (A)    root = $root
**   (B1)   depth < $depth
**   (B2)   depth <= $depth
**   (B3)   depth = $depth
**   (C)    tablename = $tablename
**   (D)    idcolumn = $idcolumn
**   (E)    parentcolumn = $parentcolumn
**
**
**
**   idxNum       meaning
**   ----------   ------------------------------------------------------
**   0x00000001   Term of the form (A) found
**   0x00000002   The term of bit-2 is like (B1)
**   0x000000f0   Index in filter.argv[] of $depth.  0 if not used.
**   0x00000f00   Index in filter.argv[] of $tablename.  0 if not used.
**   0x0000f000   Index in filter.argv[] of $idcolumn.  0 if not used
**   0x000f0000   Index in filter.argv[] of $parentcolumn.  0 if not used.
**
** There must be a term of type (A).  If there is not, then the index type
** is 0 and the query will return an empty set.
*/
static int audbBestIndex(
  sqlite3_vtab *pTab,             /* The virtual table */
  sqlite3_index_info *pIdxInfo    /* Information about the query */
){
  int iPlan = 0;
  int i;
  int idx = 1;
  int seenMatch = 0;
  const struct sqlite3_index_constraint *pConstraint;
  audb_vtab *pVtab = (audb_vtab*)pTab;
  double rCost = 10000000.0;
  pIdxInfo->estimatedCost = rCost;

  return SQLITE_OK;
}

/*
** A virtual table module that implements the "transitive_closure".
*/
static sqlite3_module audbModule = {
  0,                      /* iVersion */
  audbConnect,         /* xCreate */
  audbConnect,         /* xConnect */
  audbBestIndex,       /* xBestIndex */
  audbDisconnect,      /* xDisconnect */
  audbDisconnect,      /* xDestroy */
  audbOpen,            /* xOpen - open a cursor */
  audbClose,           /* xClose - close a cursor */
  audbFilter,          /* xFilter - configure scan constraints */
  audbNext,            /* xNext - advance a cursor */
  audbEof,             /* xEof - check for end of scan */
  audbColumn,          /* xColumn - read data */
  audbRowid,           /* xRowid - read data */
  audbUpdate,          /* xUpdate */
  0,                      /* xBegin */
  0,                      /* xSync */
  0,                      /* xCommit */
  0,                      /* xRollback */
  0,                      /* xFindMethod */
  0,                      /* xRename */
  0,                      /* xSavepoint */
  0,                      /* xRelease */
  0                       /* xRollbackTo */
};

#endif /* SQLITE_OMIT_VIRTUALTABLE */

/*
** Register the closure virtual table
*/
#ifdef _WIN32
__declspec(dllexport)
#endif
int sqlite3_audb_init(
  sqlite3 *db,
  char **pzErrMsg,
  const sqlite3_api_routines *pApi
){
  int rc = SQLITE_OK;
  SQLITE_EXTENSION_INIT2(pApi);
  (void)pzErrMsg;
#ifndef SQLITE_OMIT_VIRTUALTABLE
  rc = sqlite3_create_module(db, "audb_tree", &audbModule, 0);
#endif /* SQLITE_OMIT_VIRTUALTABLE */
  return rc;
}
