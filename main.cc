#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "Event.h"
#include "VertexCandidate.h"
#include "VertexerFast.h"

using std::vector;

void parse_args(int argc, char* argv[]);
vector<Event> load_data(char* fname);

int main(int argc, char** argv) {
  parse_args(argc, argv);
	vector<Event> events( load_data(argv[1]) );
  VertexerFast vfast = VertexerFast(64, 256, 2*13.545);
	for( Event e : events ) {
    array<float, 7> avrad = e.AvgRadii();
    vfast.LoadClustersFromEvent( e );
    vfast.SetRadii( avrad[0], avrad[1], avrad[2] );
    float xyz[3];
    vfast.FindVertex(xyz);
    #ifdef DEBUG
    std::cout<<"Event: "<<e.GetId()<<" vertex -> x: "<<xyz[0]<<" y: "<<xyz[1]<<" z: "<<xyz[2]<<std::endl;
    #endif
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