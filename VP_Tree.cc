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
#include <iostream>
using namespace std;

// Note: these typedef are used so that
// any changes to how we represent the feature can
// happen here
int FEATURE_LENGTH = 25;
string FOLDER = "Music/";
vector<string> keys_paths;
vector<string> sample_paths;
vector<string> keys_names;
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
	//static int hitCount = 0;
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
	static int index = 0;
	static bool cacheReset = false;
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
void AddKey(string path, string name = ""){
	AddtoPaths(path, &keys_paths, &keys_names, name);
}
void SetFolder(string folder){
	FOLDER = folder;
}
void BuildIndex(){
	vector<Feature> keys = GenerateFeatures(keys_paths);
	FeatureSet kSet = BuildFeatureSet(keys);
	index_root = Build_VP(kSet);
}
// return type needs to be determined
// based on what we want to return
// currently returns the name of the matched song
string SearchIndex(string sample_path){
	// the sample_path may contain multiple
	// feature windows, so we will find the one
	// that contains the best match
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
	return keys_names[min.Key_ID];
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
int main(){
	return AudioDBMain();
}
