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
    float ex, float ey, float ez, float alpha, int mcl0, int mcl1) {
  float phi = (float)atan2(-y,-x) + kPi;
  fLayers[id].x.push_back(x);
  fLayers[id].y.push_back(y);
  fLayers[id].z.push_back(z);
  fLayers[id].phi.push_back(phi);
  fLayers[id].mcl.push_back(mcl0);
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
  int counter_MC = 0;
  int id, MCl0, MCl1;
  float x, y, z, ex, ey, ez, alpha;
  std::string line;
  std::cout<<"Reading data and filling events vector."<<std::endl;
  while(std::getline(instream, line)) {
    std::istringstream inss(line);
    if(!(inss >> id >> x >> y >> z >> ex >> ey >> ez >> alpha >> MCl0 >> MCl1 )) {
      if (MCl1 > 0) counter_MC ++;
      if(id == -1) {
        Event event(evector.size());
        evector.push_back(event);
        evector.back().SetVertex(x,y,z);
      }
    } else {
      evector.back().PushHitToLayer(id, x, y, z, ex, ey, ez, alpha, MCl0, MCl1);
    }
  }
  std::cout<<"Events vector filled. MCl1 > 0: "<<counter_MC<<endl;

  return evector;
}

int numCluster(int* lut, int iPhi) {
    iPhi &= (kNphi - 1);
    return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int numTracklets(int* lut0, int* lut1, int nphi) {
    int n = 0;
    for (int i = 0; i < nphi; ++i)
        n += numCluster(lut0,i) * (numCluster(lut1,i + 1) + numCluster(lut1,i) + numCluster(lut1,i - 1));
    return n;
};


int firstTracklet(int* lut0, int* lut1, int nphi0, int nz0) {
  int n = numTracklets(lut0,lut1,nphi0);
  int mult = (numCluster(lut1,nphi0 + 1) + numCluster(lut1,nphi0) + numCluster(lut1,nphi0 - 1));
  for (int iZ0 = 0; iZ0 < nz0; ++iZ0) 
    n += numClustersByZ(lut0,nphi0,iZ0) * mult; 
  return n;
};

int numClustersByZ(int* lut, int iPhi, int iZ) {
  iPhi &= (kNphi - 1);
  return lut[iPhi*kNz + iZ + 1] - lut[iPhi*kNz + iZ];
};

