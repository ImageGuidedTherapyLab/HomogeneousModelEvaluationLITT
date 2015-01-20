#
# file : Makefile (UNIX)
#
#  mexGPUArray   from 
# 
#     $(MATLABROOT)/toolbox/distcomp/gpu/extern/src/mex/mexGPUExample.cu
#
# You can invoke this Makefile using 
#  make -f Makefile MATLABROOT=[directory where MATLAB is installed]
#
# If you do not want to specify the MATLABROOT at the gmake command line 
# every single time, you can define it by uncommenting the line below
# and assigning the correct root of MATLAB (with no trailing '/') on
# the right hand side.
#
MATLABROOT	:= /opt/apps/matlab/2013a/
#

#
# Defaults
#

MEX=$(MATLABROOT)/bin/mex
MCC=$(MATLABROOT)/bin/mcc
CUDADIR = /opt/apps/cuda/5.0/cuda
CUDADIR = /opt/apps/cuda/5.5/

# The following are the definitions for each target individually.

steadyStatePennesLaser.ptx:   steadyStatePennesLaser.cu
	$(CUDADIR)/bin/nvcc -g -G -ptx -gencode=arch=compute_20,code=sm_20   $<

tags:
	ctags -R  --langmap=c++:+.cu $(MATLABROOT) . 

