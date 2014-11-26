/*
VP_TREE.CC REGION!
*/
// VP_Tree.cc
// Vantage-Point tree implementation
// Jonathon Deriso
// CS4420 AudioDB Project
// 11/14/2014

#include <vector>
#include "stdlib.h"
#include "time.h"
#include <cstdlib>
#include "stdio.h"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include "VP_Tree.h"
// SQLITE Extension
#include "sqlite3ext.h"
SQLITE_EXTENSION_INIT1
#include <assert.h>
#include <ctype.h>
#ifndef SQLITE_OMIT_VIRTUALTABLE

// END SQLITE Extension
using namespace std;


// Note: these typedef are used so that
// any changes to how we represent the feature can
// happen here
int FEATURE_LENGTH = 25;
string FOLDER = "";
vector<string> keys_paths;
vector<string> sample_paths;
vector<string> keys_names;
vector<string> keys_artists;
vector<string> sample_names;
vector<vector<int> > vp_global_cache;
typedef long Feature_Dist;
Feature_Dist MAX_FEATURE_DIST = (~(unsigned long) 0) >> 1;
typedef vector<int> Feature;
struct Feature_t {
	//Feature* feature;
	int* feature_start;
	int Key_ID;
};
typedef vector<Feature_t> FeatureSet;

class VP_Node{
public:
	Feature_t VantagePoint;
	Feature_Dist Radius;
	VP_Node* Left;
	VP_Node* Right;
	VP_Node(VP_Node* node) : VantagePoint(node->VantagePoint), Radius(node->Radius), Left(node->Left), Right(node->Right){}
	VP_Node() {}
};
VP_Node index_root;
void internal_free_vp(VP_Node* root){
	//  int hitCount = 0;
	//cout << ++hitCount << endl;
	if (root == NULL)
		return;
	internal_free_vp(root->Left);
	internal_free_vp(root->Right);
	delete root;
}
void Free_VP(VP_Node root){
	internal_free_vp(root.Left);
	internal_free_vp(root.Right);
}

Feature_Dist Distance(int* f1, Feature f2){
	Feature_Dist dist = 0;
	for (int i = 0; i < f2.size(); i++){
		unsigned long t = f1[i] ^ f2[i];
		while (t != 0){
			dist += t & 1;
			t >>= 1;
		}
	}
	return dist;
}
Feature_Dist Distance(Feature f2, int* f1){
	return Distance(f1, f2);
}
Feature_Dist Distance(int* f1, int*f2){
	Feature_Dist dist = 0;
	for (int i = 0; i < FEATURE_LENGTH; i++){
		unsigned long t = f1[i] ^ f2[i];
		while (t != 0){
			dist += t & 1;
			t >>= 1;
		}
	}
	return dist;
}
Feature_Dist Distance(Feature f1, Feature f2){
	// Hamming Distance
	Feature_Dist dist = 0;
	for (int i = 0; i < f1.size(); i++){
		unsigned long t = f1[i] ^ f2[i];
		while (t != 0){
			dist += t & 1;
			t >>= 1;
		}
	}
	return dist;
}
Feature_Dist GetEnergy(FeatureSet S, int* P, Feature_Dist radius){
	// using variance as energy
	Feature_Dist energy, t;
	energy = 0; t = 0;
	for (int i = 0; i < S.size(); i++){
		//t = Distance(*(S[i].feature), P) - radius;
		t = Distance(S[i].feature_start, P) - radius;
		energy += (t * t);
	}
	return energy / S.size(); // will have truncate problems but close enough
}
Feature_Dist GetEnergy(FeatureSet S, Feature P, Feature_Dist radius){
	// using variance as energy
	Feature_Dist energy, t;
	energy = 0; t = 0;
	for (int i = 0; i < S.size(); i++){
		//t = Distance(*(S[i].feature), P) - radius;
		t = Distance(S[i].feature_start, P) - radius;
		energy += (t * t);
	}
	return energy / S.size(); // will have truncate problems but close enough
}
FeatureSet GetSample(FeatureSet S){
	// return random sample of Features in the set
	FeatureSet P;
	int sz = S.size();
	int sampleSize = sqrt(sz);
	if (sampleSize == 0){
		for (int i = 0; i < sz; i++){
			P[i] = S[i];
		}
		return P;
	}

	for (int i = 0; i < sampleSize; i++)
		P.push_back(S[rand() % sz]);
	return P;
}
Feature_Dist GetMedian(int* vp, FeatureSet S){
	// Return median Distance of vp to features
	// in set S
	vector<Feature_Dist> dist;
	for (int i = 0; i < S.size(); i++){
		//dist.push_back(Distance(vp, *(S[i].feature)));
		dist.push_back(Distance(vp, S[i].feature_start));
	}
	sort(dist.begin(), dist.end());
	return dist[dist.size() / 2];
}
Feature_Dist GetMedian(Feature vp, FeatureSet S){
	// Return median Distance of vp to features
	// in set S
	vector<Feature_Dist> dist;
	for (int i = 0; i < S.size(); i++){
		//dist.push_back(Distance(vp, *(S[i].feature)));
		dist.push_back(Distance(vp, S[i].feature_start));
	}
	sort(dist.begin(), dist.end());
	return dist[dist.size() / 2];
}


void internal_search(VP_Node* root, Feature sample, Feature_Dist* dist, Feature_t* match, bool isRightCheck = false){

	//If the sample feature is within the bounding radius
	//	If there is a feature in the set whose distance to the sample feature is less than the distance of the sample feature to the bounding radius then that feature is the best match
	//	ELSE we need to search the right child and see if we can find a point that is closer but still keep
	//	track of the closest item within the current set
	//	ELSE
	//	Look in the right child
	if (root == NULL){
		*dist = MAX_FEATURE_DIST;
		return;
	}
//	Feature_Dist cDist = Distance(sample, *(root->VantagePoint.feature));
	Feature_Dist cDist = Distance(sample, root->VantagePoint.feature_start);
	if (isRightCheck && cDist >= root->Radius)
		return;
	if (cDist < root->Radius){
		// within bounding radius
		// search left
		Feature_Dist tdist = cDist;
		Feature_t found;
		internal_search(root->Left, sample, &tdist, &found);
		// if we have found a point closer to the sample
		// than the sample is to the boundary of the radius
		// we know there are no points outside of the left child
		// that could be closer
		if (tdist < root->Radius - cDist){
			if (tdist < cDist){
				*dist = tdist;
				*match = found;
			} // the vantage point at this level could still be closest
			else {
				*dist = cDist;
				*match = root->VantagePoint;
			}
		}
		else { // look right, since there might be a closer point outside of this vantage point radius
			Feature_Dist rdist = tdist;
			Feature_t rfound;
			internal_search(root->Right, sample, &rdist, &rfound, true);
			if (rdist < tdist){ // if we found a closer point on the right side
				*dist = rdist;
				*match = rfound;
			}
			else { // else the closest is between the current vantage point and the closest
				// found in the left child
				if (tdist < cDist){
					*dist = tdist;
					*match = found;
				}
				else {
					*dist = cDist;
					*match = root->VantagePoint;
					
				}
			}
		}

	}
	else { // outside of bounding radius ->
		// look right
		Feature_Dist tdist = cDist;
		Feature_t found;
		internal_search(root->Right, sample, &tdist, &found);
		if (tdist < cDist){ // compare best match from right child with current vantage point
			*dist = tdist;
			*match = found;
		}
		else {
			*dist = cDist;
			*match = root->VantagePoint;
			
		}
	}
} //
Feature_t Search(VP_Node* root, Feature sample, Feature_Dist* distance){
	// wrapper function for search
	*distance = MAX_FEATURE_DIST;
	Feature_t rtn;
	if (sample.size() == 0)
		return rtn;
	internal_search(root, sample, distance, &rtn);
	return rtn;
}
Feature_t Select_VP(FeatureSet S){
	// select optimal (estimated) vantage point for
	// set S by randomly selecting a set of candidate vantage points
	// then by comparing the variance between in a set of sampled points

	//randomly sample Select_VP
	FeatureSet P;
	P = GetSample(S);
	Feature_Dist cRad = 0;
	Feature_Dist best_spread = -1;
	Feature_Dist spread = 0;
	Feature_t rtn;
	for (int i = 0; i < P.size(); i++){
		FeatureSet D;
		D = GetSample(S);
		/*cRad = GetMedian(*(P[i].feature), D);
		spread = GetEnergy(D, *(P[i].feature), cRad);*/
		cRad = GetMedian(P[i].feature_start, D);
		spread = GetEnergy(D, P[i].feature_start, cRad);
		// using variance as energy
		// want to maximize energy
		if (spread > best_spread){
			best_spread = spread;
			rtn = P[i];
		}
	}
	return rtn;
}
// Consider implementing iteratively using by emulating
// stack behavior. Stack overflow can occur relatively fast
// with the number of keys we are dealing with
VP_Node Build_VP(FeatureSet S){
	// Build Vantage Point Tree
	VP_Node node;
	if (S.size() == 0)
		return node;
	if (S.size() == 1){
		node.Left = NULL;
		node.Right = NULL;
		node.VantagePoint = S[0];
		node.Radius = 0;
		return node;
	}


	node.VantagePoint = Select_VP(S); // select Vantage Point
	//node.Radius = GetMedian(*(node.VantagePoint.feature), S); // Get radius for vantage point
	node.Radius = GetMedian(node.VantagePoint.feature_start, S); // Get radius for vantage point
	FeatureSet left, right;
	// Partition into two sets based on if there are within radius or not
	Feature_Dist prevD = -1;
	for (int i = 0; i < S.size(); i++){
		// Might need to consider not including
		// checking vantage point with itself
		// (it is still in the set)
		//Feature_Dist d = Distance(*(node.VantagePoint.feature), *(S[i].feature));
		Feature_Dist d = Distance(node.VantagePoint.feature_start, S[i].feature_start);
		if (d < node.Radius || d == 0 || (d < node.Radius && d == prevD) ){
			left.push_back(S[i]);
		}
		else {
			right.push_back(S[i]);
		}
		prevD = d;
	}

	// Need to deallocate these nodes
	if (left.size() > 0)
		node.Left = new VP_Node(Build_VP(left));
	else
		node.Left = NULL;
	if (right.size() > 0)
		node.Right = new VP_Node(Build_VP(right));
	else
		node.Right = NULL;
	return node;
}


/*
END OF VP_TREE.CC REGION
*/

// could potentially add some error handling
vector<int> readfile(const char* file){
	ifstream fs(file);

	string t;

	getline(fs, t);
	getline(fs, t);
	// start at index 9
	int len = atoi(t.substr(9, t.size() - 9).c_str());
	vector<int> out;
	getline(fs, t, ',');
	t = t.substr(12, t.size() - 12);
	out.push_back(atoi(t.c_str()));
	while (getline(fs, t, ',')){
		out.push_back(atoi(t.c_str()));
	}
	return out;
}
bool vp_cache_init_flag = false;
bool vp_cache_reset = false;
// needs to be called before building feature sets
// to initialize cache with enough memory so that
// pointers aren't invalidated when the vector must reallocate
// memory
void FeatureSetCacheInit(int numKeys){
	vp_cache_init_flag = true;
	vp_cache_reset = !vp_cache_reset;
	vp_global_cache.clear();
	vp_global_cache.resize(numKeys);
}
FeatureSet BuildFeatureSet(vector<vector<int> > keys){
	  int index = 0;
	  bool cacheReset = false;
	if (cacheReset != vp_cache_reset){
	  // Cache has been reset (cleared and resized)
	  cacheReset = vp_cache_reset;
          index = 0;
	}
	if (!vp_cache_init_flag){
		cout << "WARNING!! internal cache should be initialized using the function 'FeatureSetCacheInit(int numKeys'." << endl;
	}
	FeatureSet S;
	for (int i = 0; i < keys.size(); i++){
		//for (int k = 0; k <= keys[i].size() - FEATURE_LENGTH; k++){
		//	vector<int> f(FEATURE_LENGTH);
		//	for (int j = 0; j < FEATURE_LENGTH; j++){
		//		f[j] = keys[i][k + j];
		//	}
		//	Feature_t t;
		//	//t.feature = new Feature();
		//	//*(t.feature) = f;
		//	vp_global_cache.push_back(f);
		//	t.feature_start = &(vp_global_cache[vp_global_cache.size() - 1][0]);
		//	t.Key_ID = i;
		//	S.push_back(t);

		//}
		if (index == vp_global_cache.size()){
			cout << "WARNING!! Ran out of internal cache for feature sets, 'FeatureSetCacheInit' was not called with enough space, storage locations will be reused and loss of data will occur" << endl;
		}
		vp_global_cache[index++] = keys[i]; //.push_back(keys[i]);
		for (int k = 0; k <= keys[i].size() - FEATURE_LENGTH; k++){
			Feature_t t;
			t.feature_start = &(vp_global_cache[index - 1][k]);//&(vp_global_cache[vp_global_cache.size() - 1][k]);
			t.Key_ID = i;
			S.push_back(t);
		}
	}
	return S;
}
vector<Feature> GenerateFeatures(vector<string> paths, bool useGlobalFolder = true){
	vector<Feature> features;
	for (int i = 0; i < paths.size(); i++){
		string tm = (useGlobalFolder) ? "fpcalc -raw " + FOLDER + paths[i] + " > out_temp.txt " : "fpcalc -raw " + paths[i] + " > out_temp.txt ";
		char* ar = new char[tm.size() + 1];
		tm.copy(ar, tm.size());
		ar[tm.size()] = 0;
		system(ar);
		features.push_back(readfile("out_temp.txt"));
		delete [] ar;
	}
	return features;
}
void ClearKeys(){
	keys_names.clear();
	keys_paths.clear();
}
void ClearSamples(){
	sample_names.clear();
	sample_paths.clear();
}
void AddtoPaths(string path, vector<string>* paths, vector<string>* names, string name = ""){
	paths->push_back(path);
	if (name == "")
		names->push_back(path);// consider pruning out certain characters first (such as '"')
	else
		names->push_back(name);

}

void AddSample(string path, string name = ""){
	AddtoPaths(path, &sample_paths, &sample_names, name);
}
void Add_Sample(char* path, char* name){
	AddSample(string(path), string(name));
}
void AddKey(string path, string name = ""){
	AddtoPaths(path, &keys_paths, &keys_names, name);
}
void Add_Key(char* path, char* songName, char* artist){
	keys_paths.push_back(string(path));
	keys_names.push_back(string(songName));
	keys_artists.push_back(string(artist));
}
//void Add_Key(char* path, char* name){
//	AddKey(string(path), string(name));
//}
void SetFolder(string folder){
	FOLDER = folder;
}
void Set_Folder(char* folder){
	SetFolder(string(folder));
}
void BuildIndex(){
	vector<Feature> keys = GenerateFeatures(keys_paths);
        FeatureSetCacheInit(keys.size());
	FeatureSet kSet = BuildFeatureSet(keys);
	index_root = Build_VP(kSet);
}
// return type needs to be determined
// based on what we want to return
// currently returns the name of the matched song
int SearchIndexForKey(char* samples_path){
	// the sample_path may contain multiple
	// feature windows, so we will find the one
	// that contains the best match
	string sample_path(samples_path);
	vector<string> dummyVector;
	dummyVector.push_back(sample_path);
	vector<Feature> sampleFeatures = GenerateFeatures(dummyVector, false);
	//FeatureSet sSet = BuildFeatureSet(sampleFeatures);
	Feature_Dist dist = MAX_FEATURE_DIST;
	Feature_t min;
	for (int i = 0; i <= sampleFeatures[0].size() - FEATURE_LENGTH; i++){
		Feature_Dist d;
		Feature_t found;
     		//cout << "Entered search round " << i << endl;
		found = Search(&index_root, Feature(&(sampleFeatures[0][i]), &(sampleFeatures[0][i+ FEATURE_LENGTH])), &d);
		if (d < dist){
			dist = d;
			min = found;
		}
	}
	//cout << "distance found is " << dist << " Key found is " << min.Key_ID << endl;
	return min.Key_ID;
}
// return type needs to be determined
// based on what we want to return
// currently returns the name of the matched song
char* SearchIndex(char* samples_path, char* field){
	// the sample_path may contain multiple
	// feature windows, so we will find the one
	// that contains the best match
	string sample_path(samples_path);
	vector<string> dummyVector;
	dummyVector.push_back(sample_path);
	vector<Feature> sampleFeatures = GenerateFeatures(dummyVector, false);
	FeatureSet sSet = BuildFeatureSet(sampleFeatures);
	Feature_Dist dist = MAX_FEATURE_DIST;
	Feature_t min;
	for (int i = 0; i < sampleFeatures.size(); i++){
		Feature_Dist d;
		Feature_t found;
		found = Search(&index_root, Feature(sSet[i].feature_start, sSet[i].feature_start + FEATURE_LENGTH), &d);
		if (d < dist){
			dist = d;
			min = found;
		}
	}
	string fstr(field);
	string rstr;
	if (fstr.find("path") != string::npos || fstr == "*")
		rstr += keys_paths[min.Key_ID] + "|";
	if (fstr.find("song") != string::npos || fstr == "*")
		rstr += keys_names[min.Key_ID] + "|";
	if (fstr.find("artist") != string::npos || fstr == "*")
		rstr += keys_artists[min.Key_ID] + "|";
	//string rstr = keys_names[min.Key_ID];
	char* rtn = (char*)malloc(rstr.length() + 1);
	strcpy(rtn, rstr.c_str());
	return rtn;//keys_names[min.Key_ID].c_str();
}
int AudioDBMain(){
	for (int key_num = 1; key_num < 22; key_num++){
	cout << "Executing with " << key_num << " songs in the index" << endl;
	//int key_num = 20;//23;
	int sample_num = 10;
	string key_paths [] = { "\"GorillazWav.wav\"", "\"Martin Garrix - Animals.wav\"", "\"Everything at Once.mp3\"", "\"RHCP - Look Around.mp3\"", "\"Aoki - Boneless.mp3\"", "\"AWOLNATION - Sail.mp3\"", "\"Check My Steezo.mp3\"", "\"No Beef.mp3\"", "\"Get Up.mp3\"" , "Rattle.mp3", "Stampede.mp3", "\"Too Turnt Up.mp3\"", "\"Breakn A Sweat.mp3\"", "\"Duck Sauce.mp3\"", "Existence.mp3" , "\"Hold me Close.mp3\"", "Rage.mp3", "Freak.mp3", "\"Express Yourself.mp3\"", "\"Turn Down For What.mp3\"", "Bangarang.mp3", "\"Can't You See.mp3\"", "Turbulence.mp3"};
	string key_names [] = { "Gorillaz", "Animals", "Everything at Once", "Look Around", "Boneless" , "Sail", "Check My Steezo", "No Beef", "Get Up", "Rattle", "Stampede", "TTU", "Breakn a Sweat",  "Duck Sauce", "Existence", "Hold Me Close", "Rage", "Freak", "Express Yourself", "Turn Down For What", "Bangarang", "Can't You See", "Turbulence"};
		string sample_paths2[] = { "\"Martin2.wav\"", "\"Gorillaz2.wav\"", "\"Martin3.wav\"", "\"Gorillaz3.wav\"", "\"Everything2.wav\"", "\"Everything3.wav\"", "\"RHCP2.wav\"", "\"RHCP3.wav\"", "Aoki2.wav", "Aoki3.wav" };
	string sample_names2[] = { "Animals", "Gorillaz", "Animals", "Gorillaz", "Everything at Once", "Everything at Once", "Look Around", "Look Around", "Boneless", "Boneless" };
	ClearKeys();
	ClearSamples();
	for (int i = 0; i < key_num; i++){
		AddKey(key_paths[i], key_names[i]);
	}
	for (int i = 0; i < sample_num; i++){
		AddSample(sample_paths2[i], sample_names2[i]);
	}
	vector<vector<int> > keys = GenerateFeatures(keys_paths);

	vector<vector<int> > samples = GenerateFeatures(sample_paths);

	for (FEATURE_LENGTH = 25; FEATURE_LENGTH <= 30; FEATURE_LENGTH = 31){// FEATURE_LENGTH++){ 
		//cout << "\n\n*** Feature Length " << FEATURE_LENGTH << "***" << endl;
		FeatureSetCacheInit(keys.size() + samples.size());
		FeatureSet kSet = BuildFeatureSet(keys);
		FeatureSet sSet = BuildFeatureSet(samples);
		clock_t t1;
		time_t t2;
		long total = sSet.size();
		long hit = 0;
		long miss = 0;
		t1 = clock();
		// REGION LINEAR SEARCH
		for (int i = 0; i < sSet.size(); i++){
			int minDist = (~((unsigned) 0)) >> 1;
			int ind = -1;
			for (int j = 0; j < kSet.size(); j++){
				int t = Distance(kSet[j].feature_start, sSet[i].feature_start);
				if (t < minDist){
					minDist = t;
					ind = j;
				}
			}
			if (sample_names[sSet[i].Key_ID] == keys_names[kSet[ind].Key_ID])
				hit++;
			else
				miss++;
			//	cout << "Sample: " << sample_names[sSet[i].Key_ID] << " Song: " << keys_names[kSet[ind].Key_ID] << " Distance: " << minDist << endl;
		}
		// REGION END LINEAR SEARCH
		t2 = clock();
		cout << "{Linear Search} Execution Time: " << t2 - t1 << "   Total Samples: " << total << "   Accuracy: " << ((double) hit) / total << endl;
		t1 = clock();
		VP_Node root = Build_VP(kSet);
		t2 = clock();
		cout << "Time to build index with " << kSet.size() << " keys:  " << t2 - t1 << endl;
		hit = 0;
		miss = 0;
		t1 = clock();
		for (int i = 0; i < sSet.size(); i++){
			Feature_Dist d;
			Feature_t found;
			found = Search(&root, Feature(sSet[i].feature_start, sSet[i].feature_start + FEATURE_LENGTH), &d);
			//cout << "Sample: " << sample_names[sSet[i].Key_ID] << " Song: " << keys_names[found.Key_ID] << " Distance: " << d<< endl;
			if (sample_names[sSet[i].Key_ID] == keys_names[found.Key_ID])
				hit++;
			else
				miss++;
		}
		t2 = clock();
		cout << "{Indexed Search} Execution Time: " << t2 - t1 << "   Total Samples: " << total << "   Accuracy: " << ((double) hit) / total << endl;
		Free_VP(root);

	}
 }
	return 0;
}


typedef struct  audb_vtab  audb_vtab;
typedef struct  audb_cursor  audb_cursor;
typedef struct  audb_queue  audb_queue;
typedef struct  audb_avl  audb_avl;
#define ec extern "C"

struct  audb_avl {
  sqlite3_int64 id;     /* Id of this entry in the table */
  int iGeneration;      /* Which generation is this entry part of */
   audb_avl *pList;   /* A linked list of nodes */
   audb_avl *pBefore; /* Other elements less than id */
   audb_avl *pAfter;  /* Other elements greater than id */
   audb_avl *pUp;     /* Parent element */
  short int height;     /* Height of this node.  Leaf==1 */
  short int imbalance;  /* Height difference between pBefore and pAfter */
};

/*
** A  audb virtual-table object MAY NOT NEED
*/
 struct  audb_vtab {
  sqlite3_vtab base;         /* Base class - must be first */
  char *zDb;                 /* Name of database.  (ex: "main") */
  char *zSelf;               /* Name of this virtual table */
  char *zTableName;          /* Name of table holding parent/child relation */
  char *zIdColumn;           /* Name of ID column of zTableName */
  char *zParentColumn;       /* Name of PARENT column in zTableName */
  sqlite3 *db;               /* The database connection */
  int nCursor;               /* Number of pending cursors */
};

/* A  audb cursor object  MAY NOT NEED*/
 struct  audb_cursor {
  sqlite3_vtab_cursor base;  /* Base class - must be first */
   audb_vtab *pVtab;       /* The virtual table this cursor belongs to */
   int *pCurrent;     /* Current element of output */
   
};
/*
** xDisconnect/xDestroy method for the  audb module.
*/
ec   int  audbDisconnect(sqlite3_vtab *pVtab){
	//printf("Entered audbDisconnect\n");
  return SQLITE_OK;
}
/*
** xConnect/xCreate method for the audb module. Arguments are:
**
**   argv[0]    -> module name  ("transitive_audb")
**   argv[1]    -> database name
**   argv[2]    -> table name
**   argv[3...] -> arguments
*/

ec   int audbConnect(
  sqlite3 *db,
  void *pAux,
  int argc, const char *const*argv,
  sqlite3_vtab **ppVtab,
  char **pzErr
){
	// when called as xCreate this must create
	// a new sqlite3_vtab object and return pointer to it
	// in *ppVTab and invoke sqlite3_declare_vtab(sqlite3 *db, const char *zCreateTable)
	//printf("Entered audbConnect\n");
 int rc = SQLITE_OK;              /* Return code */
  audb_vtab *pNew = 0;          /* New virtual table */
  const char *zDb = argv[1];
  const char *zVal;
  int i;

  (void)pAux;
  *ppVtab = 0;
  pNew = (audb_vtab*)sqlite3_malloc( sizeof(*pNew) );
  if( pNew==0 ) return SQLITE_NOMEM;
  rc = SQLITE_NOMEM;
  memset(pNew, 0, sizeof(*pNew));

  printf("Building meta-data table... schema is {id, path, song, artist}\n");
  rc = sqlite3_declare_vtab(db, "CREATE TABLE x(id, path, song, artist)");

  *ppVtab = &pNew->base;
  return rc;

}

/*
** Open a new audb cursor.
*/
 ec  int audbOpen(sqlite3_vtab *pVTab, sqlite3_vtab_cursor **ppCursor){
	//printf("Entered open\n");
	//successful invocation allocates memory for sqlite3_vtab_cursor
	// initializes and makes *ppCursor point to the new object
  audb_vtab *p = (audb_vtab*)pVTab;
  audb_cursor *pCur;
  pCur = (audb_cursor*)sqlite3_malloc( sizeof(*pCur) );
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
 ec  void audbClearCursor(audb_cursor *pCur){
	//printf("entered clearcursor\n");
}

/*
** Close a audb cursor.
*/
 ec int audbClose(sqlite3_vtab_cursor *cur){
	//printf("Entered close\n");
  return SQLITE_OK;
}

/*
** Advance a cursor to its next row of output
*/
 ec int audbNext(sqlite3_vtab_cursor *cur){
	//printf("Entered next\n");
	audb_cursor* cursor = (audb_cursor*)cur;
	*(cursor->pCurrent) = -1;
	sqlite3_free(cursor->pCurrent);
	cursor->pCurrent = NULL;
  return SQLITE_OK;
}



/*
** Called to "rewind" a cursor back to the beginning so that
** it starts its output over again.  Always called at least once
** prior to any audbColumn, audbRowid, or audbEof call.
**
** This routine actually computes the audb.
**
** See the comment at the beginning of audbBestIndex() for a
** description of the meaning of idxNum.  The idxStr parameter is
** not used.
*/
bool VP_INDEX_BUILT = false;
 ec  int audbFilter(
  sqlite3_vtab_cursor *pVtabCursor,
  int idxNum, const char *idxStr,
  int argc, sqlite3_value **argv
){
  //printf("Entered filter\n");
  if (!VP_INDEX_BUILT){
	VP_INDEX_BUILT = true;
	printf("Building index...\n");
	BuildIndex();
	printf("Index built...\n");
  }
  audb_cursor* pCur = (audb_cursor*)pVtabCursor;
  int rc = SQLITE_OK;
  string path((const char*)sqlite3_value_text(argv[0]));
 // printf("argv[0] value is %s\n", path.c_str());
  char* pathArg = (char*)malloc(path.length() + 1);
  strcpy(pathArg, path.c_str());
  //printf("Searching for song that matches sample found in %s\n", pathArg);
  pCur->pCurrent  = (int*)sqlite3_malloc(sizeof(int*));
  *(pCur->pCurrent) = SearchIndexForKey(pathArg);
 // printf("Seg is after index returns\n");
  free(pathArg);

  return rc;
}

/*
**   
** 
*/
 ec int audbColumn(sqlite3_vtab_cursor *cur, sqlite3_context *ctx, int i){
	//printf("Entered column\n");
	audb_cursor* cursor = (audb_cursor*)cur;
	int keyNum = *(cursor->pCurrent);
	char* res;
	switch (i){
		case 0: // id
		//printf("%i\n", keyNum);
		sqlite3_result_int(ctx, keyNum);
		break;
		case 1: // path
		//printf("%s\n",keys_paths[keyNum].c_str());
  		res = sqlite3_mprintf(keys_paths[keyNum].c_str());
		sqlite3_result_text(ctx, res,-1, sqlite3_free);
		break;
		case 2: // song name
		//printf("%s\n",keys_names[keyNum].c_str());
  		res = sqlite3_mprintf(keys_names[keyNum].c_str());
		sqlite3_result_text(ctx, res,-1, sqlite3_free);
		//sqlite3_result_text(ctx, keys_names[keyNum].c_str(),-1,NULL);
		break;
		case 3: // artist
		//printf("%s\n",keys_artists[keyNum].c_str());
		res = sqlite3_mprintf(keys_artists[keyNum].c_str());
		//sqlite3_result_text(ctx, keys_artists[keyNum].c_str(),-1, NULL);
		sqlite3_result_text(ctx, res,-1, sqlite3_free);
		break;

	}
  return SQLITE_OK;
}

/*
** The rowid.  For the audb table, this is the same as the "id" column.
*/
  ec int audbRowid(sqlite3_vtab_cursor *cur, sqlite_int64 *pRowid){
	//printf("Entered rowid\n");
	audb_cursor* cursor = (audb_cursor*)cur;
	*pRowid = *(cursor->pCurrent);
  return SQLITE_OK;
}

ec int audbUpdate(sqlite3_vtab *pVTab, int argc, sqlite3_value **argv, sqlite_int64 *pRowid) {
   // printf("update called with argc = %i\n", argc);
	// new row inserted with rowid argv[1] and column values in argv[2] and following
	// if argv[1] is SQL NULL, generate new unique rowid	
    if (argc > 1  && (sqlite3_value_type(argv[0]) ==SQLITE_NULL || argv[0] == NULL)){ 
	int id = sqlite3_value_int(argv[2]);
	//const char* path = (const char*)sqlite3_value_text(argv[3]);
	//const char* song = (const char*)sqlite3_value_text(argv[4]);
	//const char* artist = (const char*)sqlite3_value_text(argv[5]);
	VP_INDEX_BUILT = false;
	string path((const char*)sqlite3_value_text(argv[3]));
	string song((const char*)sqlite3_value_text(argv[4]));
	string artist((const char*)sqlite3_value_text(argv[5]));
	char* pathArg = (char*)malloc(path.length() + 1);
	char* songArg = (char*)malloc(song.length() + 1);
	char* artistArg = (char*)malloc(artist.length() + 1);
	strcpy(pathArg, path.c_str());
	strcpy(songArg, song.c_str());
	strcpy(artistArg, artist.c_str());
	Add_Key(pathArg, songArg, artistArg);
	printf("Adding new key {id:%i path: %s song: %s artist: %s\n", id, pathArg, songArg, artistArg);
	free(pathArg);
	free(songArg);
	free(artistArg);
    }
    return SQLITE_OK;
}

/*
** EOF indicator
*/
  ec int audbEof(sqlite3_vtab_cursor *cur){
	//printf("Entered EOF\n");
  audb_cursor *pCur = (audb_cursor*)cur;
	int rtn = (pCur->pCurrent == NULL);
  return rtn;
}

/*
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
  ec int audbBestIndex(
  sqlite3_vtab *pTab,             /* The virtual table */
  sqlite3_index_info *pIdxInfo    /* Information about the query */
){
 // printf("Entered bestindex\n");
  if (pIdxInfo->nConstraint != 1)
	printf("Invalid sql command, needs equality constraint\n");
  const sqlite3_index_info::sqlite3_index_constraint *constraint;
  constraint = pIdxInfo->aConstraint;
  if (constraint->op != SQLITE_INDEX_CONSTRAINT_EQ)
	printf("Invalid constraint op\n");
  pIdxInfo->aConstraintUsage->argvIndex = 1;
  pIdxInfo->aConstraintUsage->omit = 1;
  return SQLITE_OK;
}

/*
** A virtual table module that implements the "transitive_audb".
*/
 sqlite3_module audbModule = {
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
** Register the audb virtual table
*/
#ifdef _WIN32
__declspec(dllexport)
#endif
 ec int sqlite3_audb_init(
  sqlite3 *db,
  char **pzErrMsg,
  const sqlite3_api_routines *pApi
){
  //printf("Entered audb_init\n");
  int rc = SQLITE_OK;
  SQLITE_EXTENSION_INIT2(pApi);
  (void)pzErrMsg;
#ifndef SQLITE_OMIT_VIRTUALTABLE
  rc = sqlite3_create_module(db, "audb_tree", &audbModule, 0);
#endif /* SQLITE_OMIT_VIRTUALTABLE */
  return rc;
}
 ec int sqlite3_extension_init(sqlite3 *db, char **pzErrMsg, const sqlite3_api_routines *pApi){
	//printf("Entered extension_init\n");
	return sqlite3_audb_init(db, pzErrMsg, pApi);
}

