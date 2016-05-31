#include "Event.h"
#include <omp.h>
#include <math.h>
#include <numeric>

#include <sstream>
#include <string>
#include <fstream>

#include <iostream>
using std::cout;
using std::endl;

Event::Event(int Id):
  fId(Id) { };

void Event::SetVertex(float x, float y, float z) {
  fMcvtx[0] = x;
  fMcvtx[1] = y;
  fMcvtx[2] = z;
};

void Event::PushHitToLayer(int id, float x, float y, float z,
    float ex, float ey, float ez, float alpha) {
  array<float, 7> Layer;
  float phi = (float)atan2(-y,-x) + kPi;
  fLayers[id].push_back( {x, y, z, phi} ); // At this point phi is note evaluated yet.
};

vector<cluster>& Event::GetClustersFromLayer(int layer) {
  return fLayers[layer];
};
void Event::Dump(int lines) {
  // cout<<"Dumping event n° "<<fId<<":"<<endl;
  // cout<<"\tVertex cordinates:"<<endl;
  // cout<<"\t\tx = "<<fMcvtx[0]<<" y = "<<fMcvtx[1]<<" z = "<<fMcvtx[2]<<endl;
  // for(int i=0; i<7; i++) {
  //   cout<<"\tFirst "<< lines <<" hits data on layer "<<i<<":"<<endl;
  //   for(int j=0; j<lines; j++) {
  //     cout<<"x: "<<fLayers[i][j].fX<<endl;
  //     cout<<"y: "<<fLayers[i][j].fY<<endl;
  //     cout<<"z: "<<fLayers[i][j].fZ<<endl;
  //     cout<<"ex: "<<fLayers[i][j].fEx<<endl;
  //     cout<<"ey: "<<fLayers[i][j].fEy<<endl;
  //     cout<<"ez: "<<fLayers[i][j].fEz<<endl;
  //     cout<<"alpha: "<<fLayers[i][j].fAlpha<<endl;

  //     cout<<endl;
  //   }
  // }
  // cout<<endl;

  return;
};

void Event::PrintVertex() {
/* #ifdef DEBUG
  cout<<"-1\t"<<fMcvtx[0]<<"\t"<<fMcvtx[1]<<"\t"<<fMcvtx[2]<<endl;
 #endif*/
  return;
};

vector<Event> load_data(char* fname)
{
  vector<Event> evector;
  std::ifstream instream;
  std::cout<<"Opening: "<<fname<<std::endl;
  instream.open(fname);
  int counter = -1;
  int id;
  float x, y, z, ex, ey, ez, alpha;
  std::string line;
  std::cout<<"Reading data and filling events vector."<<std::endl;
  while(std::getline(instream, line)) {

    // read line by line, due to the data heterogeneity.
    std::istringstream inss(line);
    if(!(inss >> id >> x >> y >> z >> ex >> ey >> ez >> alpha)) {
      if(id == -1) {
        counter++;
        Event event(counter);
        evector.push_back(event);
        evector[counter].SetVertex(x,y,z);
      }
    } else {
      evector[counter].PushHitToLayer(id, x, y, z, ex, ey, ez, alpha);
    }
  }
  std::cout<<"Events vector filled."<<std::endl;

  return evector;
}
