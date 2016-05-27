#include "Event.h"
#include <omp.h>
#include <math.h>
#include <numeric>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;

Event::Event(int Id):
  fId{Id},
  fMcvtx{},
  fRadii{},
  fLayers{}
{}

void Event::SetVertex(float x, float y, float z)
{
  fMcvtx[0] = x;
  fMcvtx[1] = y;
  fMcvtx[2] = z;
};

void Event::PushHitToLayer(int id, float x, float y, float z,
    float ex, float ey, float ez, float alpha)
{
  fLayers[id].x.push_back(x);
  fLayers[id].y.push_back(y);
  fLayers[id].z.push_back(z);
  fLayers[id].ex.push_back(ex);
  fLayers[id].ey.push_back(ey);
  fLayers[id].ez.push_back(ez);
  fLayers[id].alpha.push_back(alpha);
};

void Event::Dump(int lines)
{
/*#ifdef DEBUG
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
  cout<<endl;
#endif*/
  return;
};

void Event::PrintVertex() {
/*#ifdef DEBUG
  cout<<"-1\t"<<fMcvtx[0]<<"\t"<<fMcvtx[1]<<"\t"<<fMcvtx[2]<<endl;
#endif*/
  return;
};

