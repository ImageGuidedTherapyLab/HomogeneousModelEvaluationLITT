#include "brain3d.hpp"

class wrapbrain3d : public brain3d {
public:
  /// wrap base class constructor
  wrapbrain3d(setupAide setup):brain3d(setup){}
  /// Load temperature data from cl_mem to host
  #if OCCA_USE_OPENCL==1
  int GetNumberOfElements(){
    return elements*(Nq-1)*(Nq-1)*(Nq-1);
  }
  int GetNumberOfNodes(){
    return elements*Nq3;
  }
  void GetNodes(float *coords){
    // copied from vtu3d::addVertices
    for(int nodeid=1; nodeid<=globalX.size(); nodeid++){
      // TODO: note that fmatrix is 1-base indexing 
      // TODO: reshape the coordinates in driver routine
      coords[ 3*(nodeid-1)+0] =  globalX(nodeid);
      coords[ 3*(nodeid-1)+1] =  globalY(nodeid);
      coords[ 3*(nodeid-1)+2] =  globalZ(nodeid);
    }
    return;
  }
  void GetElements(int *connectivity){
    // copied from vtu3d::addElements
    // TODO: reshape the connectivities in driver routine
    const int numpoints = 8 ; 
    int offset[numpoints] = {0,1, Nq+1, Nq, Nq2, Nq2+1, Nq2+Nq+1, Nq2+Nq};
    int pos;
    int elementcounter = 0 ;
    for(int e=1; e<=elements; e++)
      for(int t=1; t<Nq; t++)
        for(int s=1; s<Nq; s++)
          for(int r=1; r<Nq; r++){
            pos = (r-1) + Nq*(s-1) + Nq2*(t-1) + Nq3*(e-1);
            connectivity[(numpoints+1)*elementcounter] =  numpoints;
            for(int v=0; v<8; v++)
              connectivity[(numpoints+1)*elementcounter+v+1] =  (pos + offset[v]);
            elementcounter = elementcounter + 1 ;
          }
    return;
  }
  void getHostTemperature(size_t sz, void *dest){
    assert (sz == brain_u.size);
    brain_u.toHost(sz, dest);
  }
  void setDeviceTemperature(size_t sz, void *dest){
    assert (sz == brain_u.size);
    brain_u.toDevice(sz, dest);
  }
  void getHostForcing(size_t sz, void *dest){
    assert (sz == brain_forcing.size);
    brain_forcing.toHost(sz, dest);
  }
  void setDeviceForcing(size_t sz, void *dest){
    assert (sz == brain_forcing.size);
    brain_forcing.toDevice(sz, dest);
  }
  intptr_t getTemperaturePointer() {
    return reinterpret_cast<intptr_t>(brain_u.clMem);
  }
  intptr_t setTemperaturePointer(intptr_t NewData) {
    // copy data to avoid memory leak
    cl_mem oldData =  brain_u.clMem;
    brain_u.clMem = reinterpret_cast< cl_mem >(NewData);
    return reinterpret_cast<intptr_t>(oldData);
  }
  #endif

};
