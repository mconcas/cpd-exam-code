#include "VertexerFast.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include "Event.h"

using std::cout;
using std::endl;
using std::vector;


int PartPhi(vector<cluster>& V, int lo, int hi)
{
  float x = V[lo].fP;
  int i = lo;
  int j;
  for(j=lo+1; j<hi; j++) {
    if(V[j].fP <= x) {
      i=i+1;
      std::swap(V[i],V[j]);
    }
  }
  std::swap(V[i],V[lo]);
  return i;
}


void VertexerFast::QSortPhi(vector<cluster>& V, int lo, int hi) 
{
  if(lo < hi) {
    int p = PartPhi(V, lo, hi);
    QSortPhi(V, lo, p);  
    QSortPhi(V, p+1, hi);
  }
}

VertexerFast::VertexerFast():
  fGranPhi(64), 
  fGranZ(256),
  fSizePhi(fGranPhi / 6.28318530718f),
  fSizeZ(),
  fHalfLayerLength(),
  fRLayer(),          
  fLUT(),
  fClusters() { }

VertexerFast::VertexerFast(int phigran, int zgran, float layersize):
  fGranPhi(phigran), 
  fGranZ(zgran),
  fSizePhi(phigran / 6.28318530718f), 
  fSizeZ(zgran / layersize),
  fHalfLayerLength(0.5f * layersize),
  fRLayer(),
  fLUT(),
  fClusters() { }

VertexerFast::~VertexerFast() {
  /* wow */
};

void VertexerFast::LoadClustersFromEvent( Event &event )
{
  fClusters[0].clear();
  fClusters[1].clear();
  fClusters[2].clear(); 
  fLUT[0].clear();
  fLUT[1].clear();
  fLUT[2].clear();

  // Selecting data.
#pragma omp parallel for

  for( int iL=0; iL<3; ++iL ) {
    vector<array<float, 7>> layerhits(event.GetLayerHits(iL));
    
    for( size_t iC=0; iC<layerhits.size(); ++iC ) {
      float xyz[3] = { layerhits[iC][0], layerhits[iC][1], layerhits[iC][2] };
      fClusters[iL].push_back({ xyz[0], xyz[1], xyz[2], (float)sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]), 
        (float)atan2(-xyz[1],-xyz[0])+3.14159265359f });
    }
    // Sort fCluster[j]
    VertexerFast::QSortPhi(fClusters[iL], 0, fClusters[iL].size());
    float size = 6.28318530718f / fGranPhi;
    vector<int> &tLUT = fLUT[iL];  
    tLUT.push_back(0);

    for (size_t iC = 0; iC < fClusters[iL].size(); ++iC) {
      while (fClusters[iL][iC].fP > size * tLUT.size()) {
        tLUT.push_back(iC);
      }
    }
    tLUT.push_back(fClusters[iL].size());
  }

}


void VertexerFast::FindVertex(float* xyz) {
	/// Trackleting
  cout << "Trackleting" << endl;
  int nThread = 1;
  VertexCandidate vertexCandidate;

  int nTot = 0;
  #pragma omp parallel 
  {
    nThread = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int n = 0;
    VertexCandidate candidate;

    #pragma omp for
    for (size_t iC0 = 0; iC0 < fClusters[0].size(); ++iC0) {
      int phiI = int(fClusters[0][iC0].fP * fSizePhi);
      for (int iC1 = fLUT[1][phiI]; iC1 < fLUT[1][phiI + 1]; ++iC1) {

        if (fabs(fClusters[0][iC0].fP - fClusters[1][iC1].fP) < 0.002) { 
          Line l;
          l.x[0] = fClusters[0][iC0].fX;
          l.x[1] = fClusters[0][iC0].fY;        
          l.x[2] = fClusters[0][iC0].fZ;
          l.c[0] = fClusters[1][iC1].fX - fClusters[0][iC0].fX;
          l.c[1] = fClusters[1][iC1].fY - fClusters[0][iC0].fY;
          l.c[2] = fClusters[1][iC1].fZ - fClusters[0][iC0].fZ;
          const float z = IntersectCylinder(l,fRLayer[2]);
          const int zI = int((z + fHalfLayerLength) * fSizeZ);

          for (int iC2 = fLUT[2][phiI]; iC2 < fLUT[2][phiI + 1]; ++iC2) {
            if (fabs(fClusters[2][iC2].fP - fClusters[1][iC1].fP) < 0.002 && fabs(fClusters[2][iC2].fZ - z) < 0.03) {
              candidate.Add(l);
            }
          }  
        }
      }
    }
    double localVert[3];
    #pragma omp critical
    {
      vertexCandidate.Add(candidate);      
    }
  }

  vertexCandidate.ComputeClusterCentroid();
  vertexCandidate.GetVertex(xyz);
}


float VertexerFast::IntersectCylinder(Line &l, float r) {
  const float a = l.c[0] * l.c[0] + l.c[1] * l.c[1];
  const float b = l.c[0] * l.x[0] + l.c[1] * l.x[1];
  const float c = l.x[0] * l.x[0] + l.x[1] * l.x[1] - r * r;
  const float t = (sqrt(b * b - a * c) - b) / a;
  return l.x[2] + l.c[2] * t;
}