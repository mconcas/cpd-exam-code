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
};

class Event {

  public:
    Event(int Id);
    int GetId() const { return fId; }
    void SetVertex(float x, float y, float z);
    void PrintVertex();
    void PushHitToLayer(int id, float x, float y ,float z, float ex, float ey, float ez, float alpha);
    Layer& GetLayer(int layer) { return fLayers[layer]; }
    void Dump(int=5);

  private:
    int fId;
    array<float, 3> fMcvtx;
    array<Layer, 7> fLayers;

};

vector<Event> load_data(char* fname);

#endif
