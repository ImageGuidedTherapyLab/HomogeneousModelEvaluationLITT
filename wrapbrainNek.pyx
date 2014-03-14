include "wrapbrainNek.pxi"
# --------------------------------------------------------------------

cdef inline int CHKERR(int ierr) except -1:
    if ierr != 0: raise RuntimeError
    return 0  # no error

# -------------------- wrapper for SetupAide --------------------- 
cdef class PySetupAide:
    """Interface to SetupAide class
        store a local pointer in cython
        to the the underlying C++ data structure
    """
    cdef setupAide *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self,PyStringSetupFile):
        cl_list_all_devices();
        cdef char* setupfile= PyStringSetupFile
        self.thisptr = new setupAide(string(setupfile))
        self.thisptr.append( (cython.operator.dereference(self.thisptr))[string("CASE FILE")] ) 
# -------------------- wrapper for brain3d --------------------- 
cdef class PyBrain3d:
    cdef wrapbrain3d *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self,PySetupAide setup not None):
        """
         wrapper to brain3d class
        """
        self.thisptr = new wrapbrain3d( cython.operator.dereference(setup.thisptr) )
    def timeStep( self,currentTime):
        """
         advance solution 1 timestep
        """
        return self.thisptr.timeStep(currentTime)
    def heatStep( self,currentTime):
        """
         advance solution 
        """
        self.thisptr.heatStep(currentTime)
    def GetNumberOfNodes( self):
        """
         get vtu number of nodes
        """
        return self.thisptr.GetNumberOfNodes()
    def GetNumberOfElements( self):
        """
         get vtu number of elements
        """
        return self.thisptr.GetNumberOfElements()
    def screenshot( self, double CurrentTime ):
        """
         write vtu file
        """
        self.thisptr.screenshot( <float> CurrentTime ) 
    def dt( self):
        """
         time step
        """
        return self.thisptr.dt
    def getHostTemperature(  self,np.ndarray[float, ndim=1, mode="c"] Temperature not None):
        """
        transfer data to host arrary
        """
        assert Temperature.dtype == np.float32 
        cdef size_t databyte = Temperature.shape[0] * 4
        self.thisptr.getHostTemperature(databyte,&Temperature[0])
    def setDeviceTemperature(self,np.ndarray[float, ndim=1, mode="c"] Temperature not None):
        """
        transfer data to device arrary
        """
        assert Temperature.dtype == np.float32 
        cdef size_t databyte = Temperature.shape[0] * 4
        self.thisptr.setDeviceTemperature(databyte,&Temperature[0])
    def getHostForcing(  self,np.ndarray[float, ndim=1, mode="c"] Forcing not None):
        """
        transfer data to host arrary
        """
        assert Forcing.dtype == np.float32 
        cdef size_t databyte = Forcing.shape[0] * 4
        self.thisptr.getHostForcing(databyte,&Forcing[0])
    def setDeviceForcing(self,np.ndarray[float, ndim=1, mode="c"] Forcing not None):
        """
        transfer data to device arrary
        """
        assert Forcing.dtype == np.float32 
        cdef size_t databyte = Forcing.shape[0] * 4
        self.thisptr.setDeviceForcing(databyte,&Forcing[0])
    def GetNodes(self,np.ndarray[float, ndim=1, mode="c"] nodesarray not None):
        """
        transfer nodes array as 1-d array
        """
        self.thisptr.GetNodes(&nodesarray[0])
    def GetElements(self,np.ndarray[int, ndim=1, mode="c"] elementarray not None):
        """
        transfer element connectivity as 1-d array
        """
        self.thisptr.GetElements(&elementarray[0])
    def PrintSelf(self):
        """
        print class infor
        """
        self.thisptr.PrintSelf()

