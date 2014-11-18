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

// Note: this typedef is used so that
// any changes to how we represent the feature can
// happen here

typedef long Feature_Dist;
Feature_Dist MAX_FEATURE_DIST = (~(unsigned long) 0) >> 1;
typedef vector<int> Feature;
struct Feature_t {
	Feature* feature;
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
void Free_VP(VP_Node* root){
	// not currently being used
	if (root == NULL)
		return;
	Free_VP(root->Left);
	Free_VP(root->Right);
	delete root;
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

Feature_Dist GetEnergy(FeatureSet S, Feature P, Feature_Dist radius){
	// using variance as energy
	Feature_Dist energy, t;
	energy = 0; t = 0;
	for (int i = 0; i < S.size(); i++){
		t = Distance(*(S[i].feature), P) - radius;
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
Feature_Dist GetMedian(Feature vp, FeatureSet S){
	// Return median Distance of vp to features
	// in set S
	vector<Feature_Dist> dist;
	for (int i = 0; i < S.size(); i++){
		dist.push_back(Distance(vp, *(S[i].feature)));
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
	Feature_Dist cDist = Distance(sample, *(root->VantagePoint.feature));
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
		cRad = GetMedian(*(P[i].feature), D);
		spread = GetEnergy(D, *(P[i].feature), cRad);
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
	node.Radius = GetMedian(*(node.VantagePoint.feature), S); // Get radius for vantage point
	FeatureSet left, right;
	// Partition into two sets based on if there are within radius or not
	for (int i = 0; i < S.size(); i++){
		// Might need to consider not including
		// checking vantage point with itself
		// (it is still in the set)
		Feature_Dist d = Distance(*(node.VantagePoint.feature), *(S[i].feature));
		if (d< node.Radius || d == 0){
			left.push_back(S[i]);
		}
		else {
			right.push_back(S[i]);
		}
	}
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
int FEATURE_LENGTH = 15;
int dist(vector<int> d1, vector<int> d2, int min = -1){
	min= min == -1 ? (d1.size() < d2.size()) ? d1.size() : d2.size() : min;
	int dist = 0;
	for (int i = 0; i < min; i++){
		unsigned t = d1[i] ^ d2[i];
		while (t != 0){
			if (t & 1)
				dist++;
			t >>= 1;
		}
	}
	return dist;

}
vector<int> readfile(const char* file){
	ifstream fs(file);

	string t;

	getline(fs, t);
	getline(fs, t);
	// start at index 9
	int len = stoi(t.substr(9, t.size() - 9));
	vector<int> out;
	getline(fs, t, ',');
	t = t.substr(12, t.size() - 12);
	out.push_back(stoi(t));
	while (getline(fs, t, ',')){
		out.push_back(stoi(t));
	}
	return out;
}
int min_dist(vector<int> key, vector<int> sample){
	int mdist = (~((unsigned) 0)) >> 1;
	int upper = key.size() - sample.size();
	for (int i = 0; i <= upper; i++){
		vector<int> t(sample.size());
		for (int j = 0; j < sample.size(); j++){
			t[j] = key[i + j];
		}
		int td = dist(t, sample);
		if (td < mdist)
			mdist = td;
	}
	return mdist;
}
FeatureSet BuildFeatureSet(vector<vector<int> > keys){
	FeatureSet S;
	for (int i = 0; i < keys.size(); i++){
		for (int k = 0; k <= keys[i].size() - FEATURE_LENGTH; k++){
			vector<int> f(FEATURE_LENGTH);
			for (int j = 0; j < FEATURE_LENGTH; j++){
				f[j] = keys[i][k + j];
			}
			Feature_t t;
			t.feature = new Feature();
			*(t.feature) = f;
			t.Key_ID = i;
			S.push_back(t);
		}
	}
	return S;
}
int AudioDBMain(){
	int KEYS = 5;
	int SAMPLES = 10;
	string FOLDER = "Music/";
	string keys_paths [] = { "\"GorillazWav.wav\"", "\"Martin Garrix - Animals.wav\"", "\"Everything at Once.mp3\"", "\"RHCP - Look Around.mp3\"", "\"Aoki - Boneless.mp3\"" };
	string sample_paths [] = { "\"Martin2.wav\"", "\"Gorillaz2.wav\"", "\"Martin3.wav\"", "\"Gorillaz3.wav\"", "\"Everything2.wav\"", "\"Everything3.wav\"", "\"RHCP2.wav\"", "\"RHCP3.wav\"", "Aoki2.wav", "Aoki3.wav" };
	string keys_names [] = { "Gorillaz", "Animals", "Everything at Once", "Look Around", "Boneless" };
	string sample_names [] = { "Animals", "Gorillaz", "Animals", "Gorillaz", "Everything at Once", "Everything at Once", "Look Around", "Look Around", "Boneless", "Boneless" };
	vector<vector<int> > keys;
	for (int i = 0; i < KEYS; i++){
		string tm = "fpcalc -raw " + FOLDER + keys_paths[i] + " > out_temp.txt ";
		char* ar = new char[tm.size() + 1];
		tm.copy(ar, tm.size());
		ar[tm.size()] = 0;
		system(ar);
		keys.push_back(readfile("out_temp.txt"));
		delete [] ar;
	}
	vector<vector<int> > samples;
	for (int i = 0; i < SAMPLES; i++){
		string tm = "fpcalc -raw " + FOLDER + sample_paths[i] + " > out_temp.txt ";
		char* ar = new char[tm.size() + 1];
		tm.copy(ar, tm.size());
		ar[tm.size()] = 0;
		system(ar);
		samples.push_back(readfile("out_temp.txt"));

		delete [] ar;
	}
	for (; FEATURE_LENGTH <= 30; FEATURE_LENGTH += 5){
		cout << "\n\n*** Feature Length " << FEATURE_LENGTH << "***" << endl;
		FeatureSet kSet = BuildFeatureSet(keys);
		FeatureSet sSet = BuildFeatureSet(samples);
		clock_t t1;
		time_t t2;
		long total = sSet.size();
		long hit = 0;
		long miss = 0;
		t1 = clock();
		for (int i = 0; i < sSet.size(); i++){
			int minDist = (~((unsigned) 0)) >> 1;
			int ind = -1;
			for (int j = 0; j < kSet.size(); j++){
				int t = dist(*(kSet[j].feature), *(sSet[i].feature));
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
			found = Search(&root, *(sSet[i].feature), &d);
			//cout << "Sample: " << sample_names[sSet[i].Key_ID] << " Song: " << keys_names[found.Key_ID] << " Distance: " << d<< endl;
			if (sample_names[sSet[i].Key_ID] == keys_names[found.Key_ID])
				hit++;
			else
				miss++;
		}
		t2 = clock();
		cout << "{Indexed Search} Execution Time: " << t2 - t1 << "   Total Samples: " << total << "   Accuracy: " << ((double) hit) / total << endl;
		//Free_VP(&root);
	}
	return 0;
}

int main() {
    AudioDBMain();
}
