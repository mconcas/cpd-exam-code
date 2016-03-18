//////////////////////////////////////////////////////////////////////////////
// Event class to store event-related data:
// Each event consists of:
//
//  - An ID number.
//  - Montecarlo generated vertex.
//  - 7 layers (represented by vectors) containig "hits": array of 3D points.
//    with their errors, and an "alpha paramter".
//

#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include <array>

using std::array;
using std::vector;

class Event {

	public:
		Event(int Id);
		virtual ~Event();
 		void SetVertex(float x, float y, float z);
 		void PrintVertex();
 		void PushHitToLayer(int id, float x, float y ,float z,
 			float ex, float ey, float ez, float alpha);            // TODO: find what does alpha represent.
    void Dump(int=5);                                              // Dump event.

	private:
		int fId;                                                   // Id number
		array<float,3>          fMcvtx;                            // Monte Carlo truth.
		vector<array<float, 7>> fLayers[7];                        // Array of layers

};

#endif
