#ifndef VERTEXCANDIDATE_H
#define VERTEXCANDIDATE_H
#include <iostream>

struct Line {
  double x[3];                     /// Start point of the line
  double c[3];                     /// Direction of the line
};

float IntersectCylinder(Line &l, float r);

class VertexCandidate {
  public:
    VertexCandidate();

    void Add(Line &line);
    void Add(VertexCandidate &cand);
    void ComputeClusterCentroid();
    inline int GetSize() const { return fSize; }
    inline void GetA(double a[6]) { for (int i = 0; i < 6; ++i) a[i] = fA[i]; }
    inline void GetB(double b[3]) { for (int i = 0; i < 3; ++i) b[i] = fB[i]; }
    void GetCovMatrix(double cov[6]);
    template<typename T> inline void GetVertex(T p[3]) { for(int i=0; i<3; ++i) p[i]=fV[i]; }
    inline void GetWeight(double w[9]) { for(int i=0; i<9; ++i) w[i]=fW[i]; }

  protected:
    double fA[6];         // AX=B weight matrix
    double fB[3];         // AX=B
    double fV[3];         // vertex candidate
    double fW[9];         // weight matrix
    int    fSize;
};

#endif
