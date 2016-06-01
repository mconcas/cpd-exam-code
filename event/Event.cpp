#include "Event.h"
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
  float phi = (float)atan2(-y,-x) + kPi;
  fLayers[id].x.push_back(x);
  fLayers[id].y.push_back(y);
  fLayers[id].z.push_back(z);
  fLayers[id].phi.push_back(phi);
};

void Event::PrintVertex() {
  cout<<"-1\t"<<fMcvtx[0]<<"\t"<<fMcvtx[1]<<"\t"<<fMcvtx[2]<<endl;
  return;
};

vector<Event> load_data(char* fname)
{
  vector<Event> evector;
  std::ifstream instream;
  std::cout<<"Opening: "<<fname<<std::endl;
  instream.open(fname);
  int id;
  float x, y, z, ex, ey, ez, alpha;
  std::string line;
  std::cout<<"Reading data and filling events vector."<<std::endl;
  while(std::getline(instream, line)) {
    std::istringstream inss(line);
    if(!(inss >> id >> x >> y >> z >> ex >> ey >> ez >> alpha)) {
      if(id == -1) {
        Event event(evector.size());
        evector.push_back(event);
        evector.back().SetVertex(x,y,z);
      }
    } else {
      evector.back().PushHitToLayer(id, x, y, z, ex, ey, ez, alpha);
    }
  }
  std::cout<<"Events vector filled."<<std::endl;

  return evector;
}
