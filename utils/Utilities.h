#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "Event.h"

using std::vector;
using std::cout;
using std::endl;

// void parse_args(int argc, char* argv[]);
// vector<Event> load_data(char* fname);

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
