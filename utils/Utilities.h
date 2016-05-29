#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "Event.h"

using std::vector;
using std::cout;
using std::endl;
using std::setw;

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

template<typename T> void print_elm(T t, const int& width)
{
    cout<< left << setw(width) << t;
} 

void DumpLUT(const vector<int>& LUT, const int size_x) {
  for (int i = 0; i < LUT.size(); ++i)
  { 
    if ( i !=0 && i % size_x == 0 ) cout<<endl;
    cout<< setw(10) << LUT[i];
  }
  cout<<endl<<endl<<endl;
}
