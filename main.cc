#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Event.h"
#include "MeanRadius.h"
#include "Definitions.h"

using std::vector;
using std::cout;
using std::endl;

void parse_args(int argc, char* argv[]);
vector<Event> load_data(char* fname);


#define SIZE (1024)

int main(int argc, char** argv) {

  parse_args(argc, argv);
  vector<Event> events( load_data(argv[1]) );
  Event &e = events[0];
  for ( Event& e : events ) {
    cluster* layer_clusters[7] = {nullptr};
    vector<float> phiv[7];
    float* phi[7] = {nullptr};
    float radius[7] = {0.f};
    int   size[7] = {0};

    auto index = [&](const cluster& i, int l) {
      return (i.fP * kInvDphi) * kNz + ((i.fZ + kZ[l]) * kInvDz[l]);
    };

    for (int iL = 0; iL < 7; ++iL) {
      layer_clusters[iL] = e.GetClustersFromLayer(iL).data();
      size[iL] = e.GetClustersFromLayer(iL).size();
      /* Cluster sort */
      std::sort(e.GetClustersFromLayer(iL).begin(),e.GetClustersFromLayer(iL).end(),[&](const cluster& i, const cluster& j) {
            return index(i,iL) < index(j,iL);
          });
      /* Lookup table fill */
#pragma offload_transfer target(mic:0) in(layer_clusters[iL] : length(size[iL]) ALLOC RETAIN) signal(layer_clusters[iL])
#pragma offload target(mic:0) nocopy(layer_clusters[iL] : length(size[iL]) REUSE RETAIN) wait(layer_clusters[iL])
      {
        for (int iC = 0; iC < size[iL]; ++iC) {
          radius[iL] += sqrt(layer_clusters[iL][iC].fX * layer_clusters[iL][iC].fX + layer_clusters[iL][iC].fY * layer_clusters[iL][iC].fY);
        }
      }
    }
  }

  return 0;
}


void parse_args(int argc, char** argv)
{
  if( argv[1] == NULL ) {
    std::cerr<<"Please, provide a data file."<<std::endl;
    exit(EXIT_FAILURE);
  }
}

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
