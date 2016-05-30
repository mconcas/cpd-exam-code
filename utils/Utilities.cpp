#include "Utilities.h"
#include <iostream>

/*#include <iomanip>
using std::cout;
using std::endl;
using std::setw;*/

int get_nclusters(const int* lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int n_tracklets(const int* lut0, const int* lut1, int nphi) {
  int n = 0;
  for (int i = 0; i < nphi; ++i)
    n += get_nclusters(lut0,i) * (get_nclusters(lut1,i + 1) + get_nclusters(lut1,i) + get_nclusters(lut1,i - 1));
  return n;
};

void parse_args(int argc, char** argv)
{
  if( argv[1] == NULL ) {
    std::cerr<<"Please, provide a data file."<<std::endl;
    exit(EXIT_FAILURE);
  }
}

/*void DumpLUT(const vector<int>& LUT, const int size_x) {
  for (int i = 0; i < LUT.size(); ++i)
  {
    if ( i !=0 && i % size_x == 0 ) cout<<endl;
    cout<< setw(10) << LUT[i];
  }
  cout<<endl<<endl<<endl;
}*/
