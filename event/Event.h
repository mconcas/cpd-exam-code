//////////////////////////////////////////////////////////////////////////////
// Event class to store event-related data:
// Each event consists of:
//
//  - An ID number.
//  - Montecarlo generated vertex.
//  - 7 layers (represented by vectors) containig "hits": array of 3D points.
//    with their errors, and an "alpha parameter".
//

#ifndef EVENT_H
#define EVENT_H

#include <vector>

using std::vector;

struct Layer {
  vector<float> x;
  vector<float> y;
  vector<float> z;
  vector<float> ex;
  vector<float> ey;
  vector<float> ez;
  vector<float> alpha;
};

class Event {

	public:
		Event(int Id);
		int GetId() const { return fId; }

    float GetMCVertex(int i) { return fMcvtx[i]; }

    void SetVertex(float x, float y, float z);
    void SetRadii(int iL, float radius) { fRadii[iL] = radius; }
 		void PrintVertex();
 		void PushHitToLayer(int id, float x, float y ,float z,
 			float ex, float ey, float ez, float alpha);

    Layer& GetLayer(const int idx) { return fLayers[idx]; }
    float AvgRadii(int i) { return fRadii[7]; }
    void  Dump(int=5);

	private:
		int fId;             // Id number
		float fMcvtx[3];     // Monte Carlo truth.
    float fRadii[7];     // Mean radii
		Layer fLayers[7];    // Array of layers
};

#endif
