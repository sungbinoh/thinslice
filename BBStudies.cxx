#include "BetheBloch.h"
#include <iostream>

using namespace std;

int main(){

  BetheBloch bb(13);

  cout<<bb.meandEdx(1000)<<" "<<bb.MPVdEdx(1000, 0.5)<<endl;

  const int np = 13;
  double spline_KE[np] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000};
  for (int i = 0; i<np; ++i){
    cout<<spline_KE[i]<<" "<<bb.RangeFromKE(spline_KE[i])<<endl;
  }

  return 0;

}
