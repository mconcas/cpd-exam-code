#ifndef UTILITIES_H
#define UTILITIES_H

#define kPi    (3.14159265359f)
#define kTwoPi (2.f * 3.14159265359f)

#define kNz   16
#define kNphi 256

struct cluster {
  float fX;
  float fY;
  float fZ;
  float fP;
};

struct tracklet {
  int i0,i1;
  float dzdr;
  int magic;
};

/*int get_nclusters(const int* lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int n_tracklets(const int* lut0, const int* lut1, int nphi = kNphi) {
  int n = 0;
  for (int i = 0; i < nphi; ++i)
    n += get_nclusters(lut0,i) * (get_nclusters(lut1,i + 1) + get_nclusters(lut1,i) + get_nclusters(lut1,i - 1));
  return n;
};*/

#endif
