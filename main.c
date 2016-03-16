#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "Event.h"

#define DEBUG 1

int main(int argc, char** argv) {
  if( argv[1] == NULL ) {
    std::cerr<<"Please, provide a data file."<<std::endl;
    return 1;
  }

	std::ifstream instream;
	std::cout<<"Opening: "<<argv[1]<<std::endl;
	instream.open(argv[1]);

	std::vector<Event> events;
	int counter = -1;
	int id;
	float x, y, z, ex, ey, ez, alpha;
	std::string line;
  std::cout<<"Reading data and filling events vector."<<std::endl;
    while(std::getline(instream, line)) {
    	// read line by line, because of the data heterogeneity.
    	std::istringstream inss(line);
    	if(!(inss >> id >> x >> y >> z >> ex >> ey >> ez >> alpha)) {
   			if(id == -1) {
   				counter++;
  				Event event(counter);
  				events.push_back(event);
  				events.at(counter).SetVertex(x,y,z);
   			} 
      } else {
   				events.at(counter).PushHitToLayer(id, x, y, z, ex, ey, ez, alpha);
    		}
    	}
  std::cout<<"Events vector filled."<<std::endl;

  for( int i=0; i<(int)events.size(); i++) {
    events[i].Dump();
    std::cout<<std::endl;
  }  

	return 0;
  
}