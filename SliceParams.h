#ifndef SLICEPARAMS_H
#define SLICEPARAMS_H

//const int nwires_in_slice = 20;
//const int nslices = 480/nwires_in_slice;

namespace pi{
  const double thinslicewidth = 10; //cm
  const double Eslicewidth = 50; //MeV
  const double plim = 1000;
  const int nthinslices = 20;
  
  const int nbinse=100;
  const int nbinsthickness = 100;
}

namespace p{
  const double thinslicewidth = 4; //cm
  const int nthinslices = 24;
  
  const int nbinse=12; 
  const int nbinsthickness = 100;
}


#endif
