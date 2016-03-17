#include "Event.h"

#ifdef DEBUG
#include <iostream>
using std::cout;
using std::endl;
#endif

Event::Event(int Id):
  fId(Id) { };

Event::~Event() {
  /* yee */
};

void Event::SetVertex(float x, float y, float z) {
  fMcvtx = {x, y, z};
};

void Event::PushHitToLayer(int id, float x, float y, float z,
    float ex, float ey, float ez, float alpha) {
  fLayers[id].push_back( {x, y, z, ex, ey, ez, alpha} );
};

void Event::Dump(int lines) {
#ifdef DEBUG
  cout<<"Dumping event n° "<<fId<<":"<<endl;
  cout<<"\tVertex cordinates:"<<endl;
  cout<<"\t\tx = "<<fMcvtx[0]<<" y = "<<fMcvtx[1]<<" z = "<<fMcvtx[2]<<endl;
  for(int i=0; i<7; i++) {
    cout<<"\tFirst "<< lines <<" hits data on layer "<<i<<":"<<endl;
    for(int j=0; j<lines; j++) {
      cout<<"\t\t";
      for( float x : fLayers[i][j] ) cout<< x <<"\t";
      cout<<endl;
    }
  }
#endif
  return;
};

void Event::PrintVertex() {
#ifdef DEBUG
  cout<<"-1\t"<<fMcvtx[0]<<"\t"<<fMcvtx[1]<<"\t"<<fMcvtx[2]<<endl;
#endif
  return;
};

