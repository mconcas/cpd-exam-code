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
#include <array>
#include "definitions.h"

using std::array;
using std::vector;

struct Layer {
  vector<float> x;
  vector<float> y;
  vector<float> z;
  vector<float> phi;
  vector<int> mcl;
};

class Event {

  public:
    Event(int Id);
    int GetId() const { return fId; }
    void SetVertex(float x, float y, float z);
    void PrintVertex();
    array<float,3> GetVertex() const { return fMcvtx; }
    void PushHitToLayer(int id, float x, float y ,float z, float ex, float ey, float ez,
      float alpha, int mcl0, int mcl1);
    Layer& GetLayer(int layer) { return fLayers[layer]; }
    void Dump(int=5);

  private:
    int fId;
    array<float, 3> fMcvtx;
    array<Layer, 7> fLayers;

};

vector<Event> load_data(char* fname);
int numCluster(int* lut, int iPhi);
int numTracklets(int* lut0, int* lut1, int nphi);
int firstTracklet(int* lut0, int* lut1, int nphi0, int nz0);
int numClustersByZ(int* lut, int iPhi, int iZ);
#endif
