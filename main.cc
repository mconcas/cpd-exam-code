#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "Event.h"

#ifdef __APPLE__
#include <OpenCL/cl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif
#define DEBUG 1

void parse_arguments(int argc, char *argv[]);

int main(int argc, char** argv) {
  parse_arguments(argc, argv);

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
        events[counter].SetVertex(x,y,z);
   		} 
    } else {
   		events[counter].PushHitToLayer(id, x, y, z, ex, ey, ez, alpha);
   	}
  }
  std::cout<<"Events vector filled."<<std::endl;

  // Dummy debug print.
  for(Event e : events) e.Dump();

	return 0;
}

void parse_arguments(int argc, char *argv[])
{
  if( argv[1] == NULL ) {
    std::cerr<<"Please provide a valid data file."<<std::endl;
    exit(EXIT_FAILURE);
  }
}