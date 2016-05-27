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
#include <omp.h>

using std::array;
using std::vector;

#pragma offload_attribute(push,target(mic))
struct cluster {
  float fX;
  float fY;
  float fZ;
  float fP;
};
#pragma offload_attribute(pop)

class Event {

  public:
    Event(int Id);
    int GetId() const { return fId; }
    void SetVertex(float x, float y, float z);
    void PrintVertex();
    void PushHitToLayer(int id, float x, float y ,float z, float ex, float ey, float ez, float alpha);
    vector<cluster>& GetClustersFromLayer(int layer);
    void Dump(int=5);                    
    
  private:
    int fId;                            
    array<float, 3>          fMcvtx;     
    vector<cluster> fLayers[7];          

};

#endif