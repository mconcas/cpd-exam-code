#include "Event.h"
#define DEBUG 0

Event::Event(int Id):
  fId(Id) { };

Event::~Event() 
{
  /* yee */
};

void Event::SetVertex(float x, float y, float z) 
{
  fMcvtx = {x, y, z};
};

void Event::PrintVertex() 
{
  std::cout<<"-1\t"<<fMcvtx[0]<<"\t"<<fMcvtx[1]<<"\t"<<fMcvtx[2]<<std::endl;
};

void Event::PushHitToLayer(int id, float x, float y, float z, 
  float ex, float ey, float ez, float alpha) 
{
  fLayers[id].push_back( {x, y, z, ex, ey, ez, alpha} );  
};

void Event::Dump(const int lines) 
{ 
  std::cout<<"Dumping event n° "<<fId<<":"<<std::endl;
  std::cout<<"\tVertex cordinates:"<<std::endl;
  std::cout<<"\t\tx = "<<fMcvtx[0]<<" y = "<<fMcvtx[1]<<" z = "<<fMcvtx[2]<<std::endl;
  for(int i=0; i<7; i++) {
    std::cout<<"\tFirst "<< lines <<" hits data on layer "<<i<<":"<<std::endl;
    for(int j=0; j<lines; j++) {
      std::cout<<"\t\t";
      for( float x : fLayers[i][j] ) std::cout<< x <<"\t";
      std::cout<<std::endl;
    }
  }
  std::cout<<std::endl;
  return; 
};

