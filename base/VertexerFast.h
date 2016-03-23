#ifndef VERTEXERFAST_H
#define VERTEXERFAST_H

#include "VertexCandidate.h"
#include "Event.h"
#include <omp.h>
#include <vector>


using std::vector;
#endif

struct cluster
{
	float fX;
	float fY;
	float fZ;
	float fR;
	float fP;

	cluster() : fX(), fY(), fZ(), fR(), fP() {};
	cluster(float x, float y, float z, float r, float p): 
		fX(x), fY(y), fZ(z), fR(r), fP(p) {}
};

int PartPhi(vector<float>&, int lo, int hi);

class VertexerFast
{
	public:
		VertexerFast();
		VertexerFast(int phigran, int zgran, float layersize);
		virtual ~VertexerFast();
		void LoadClustersFromEvent(Event &event);
		static void QSortPhi(vector<cluster>&, int, int);  // QuickSort using phi as key.
		// void LoadClusters(TClonesArray *clusters[3]);
		void FindVertex(float *xyz);
		void SetRadii(float r0, float r1, float r2) { fRLayer[0]=r0; fRLayer[1]=r1; fRLayer[2]=r2; }
		void SetNumberOfThread(int n) { omp_set_num_threads(n); }	
		
	
	private:
		float IntersectCylinder(Line &l, float r);
		int fGranPhi;
		int fGranZ;
		float fSizePhi;
		float fSizeZ;             
		float fHalfLayerLength;   
		float fRLayer[3];
		vector<bool>  fMask;
		vector<int>   fLUT[3];
		vector<cluster> fClusters[3];
		
};