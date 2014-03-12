cimport cython 

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
import numpy as np
cimport numpy as np


#from libcpp.vector   cimport vector
#from libc.stdint  cimport intptr_t

#TODO need to use datafloat type
ctypedef float brainNekdatafloat
# -------------------- std::string interface ----------------------- 
cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()
# -------------------- wrapper for opencl utilities ----------------------- 
cdef extern from "occa.hpp": 
    cdef void cl_list_all_devices()
# -------------------- wrapper for setupAide ----------------------- 
cdef extern from "setupAide.hpp": 
    cdef cppclass setupAide:
        setupAide(string)
        void append(string)
        string operator[](string)
# ------------------- wrapper for brain3d ----------------------- 
cdef extern from "wrapbrainNek.h": 
    cdef cppclass wrapbrain3d:
        wrapbrain3d(setupAide) 
        int GetNumberOfElements()
        int GetNumberOfNodes()
        void GetElements(int*)
        void GetNodes(brainNekdatafloat*)
        int timeStep(brainNekdatafloat)
        brainNekdatafloat dt
        void screenshot(  brainNekdatafloat )
        void getHostTemperature(  size_t , void *)
        void setDeviceTemperature(size_t , void *)
        void getHostForcing(  size_t , void *)
        void setDeviceForcing(size_t , void *)
        void PrintSelf( )
#        # http://documen.tician.de/pyopencl/misc.html#interoperability-with-other-opencl-software
#        intptr_t getTemperaturePointer()
#        intptr_t setTemperaturePointer(intptr_t)
