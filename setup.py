from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import subprocess
import os

MacOSXBuild = False 
MacOSXBuild = True 

#FIXME - can we automate the paths ? 
if ( MacOSXBuild ):
  brainnek_dir = '/Users/fuentes/MyProjects/braincode/tym1'
  print "export PYTHONPATH=/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5-build/Wrapping/Python/:/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5-build/bin/:/Users/fuentes/github/DakotaApplications/PlanningValidation/build.tmp/"
else:
  brainnek_dir = '/workarea/fuentes/braincode/tym1'


build_dir   = "build.tmp" 

# local includes
brainnek_include = []
brainnek_include.append( "."  )
brainnek_include.append( "%s/src"      % brainnek_dir )
brainnek_include.append( "%s/include"  % brainnek_dir )
brainnek_include.append( "%s/libocca"  % brainnek_dir )


if ( MacOSXBuild ):
  #FIXME - can we automate the paths ? 
  brainnek_include.append( "/Applications/Xcode.app/Contents//Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk//System/Library/Frameworks/OpenCL.framework/" )
  brainnek_include.append( "/Applications/Xcode.app/Contents//Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include/")
  brainnek_include.append( "/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5/Common/" ) 
  brainnek_include.append( "/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5/Filtering/" ) 
  brainnek_include.append( "/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5/IO/" ) 
  brainnek_include.append( "/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5/Hybrid/" ) 
  brainnek_include.append( "/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/VTKv5-build/" ) 
else:
  brainnek_include.append( "/opt/apps/khronos/1.1")
  brainnek_include.append( "/opt/apps/EPD/epd-7.3-1-rh5-x86_64/lib/python2.7/site-packages/numpy/core/include")
  brainnek_include.append( "/opt/apps/EPD/epd-7.3-1-rh5-x86_64/include/vtk-5.6")


# cxx flags
brainnek_cxxflags = []
# FIXME why is this flag needed ??
brainnek_cxxflags.append( "-DCYTHON_CCOMPLEX=0"     )

brainnek_cxxflags.append( "-DMPI_ENABLED=0"      )
brainnek_cxxflags.append( "-Ddatafloat=float"    )
brainnek_cxxflags.append( "-DGL_ENABLED=0"       )
brainnek_cxxflags.append( "-DSCREENSHOT=1"       )
brainnek_cxxflags.append( "-DVERBOSE=0"          )
if ( MacOSXBuild ):
  brainnek_cxxflags.append( "-DOS_OSX=1"         )
else:
  brainnek_cxxflags.append( "-DOS_LINUX=1"         )
brainnek_cxxflags.append( "-DOCCA_USE_OPENCL=1"  )
brainnek_cxxflags.append( "-DOCCA_USE_CUDA=0"    )
brainnek_cxxflags.append( "-DOCCA_USE_ALL=0"     )
brainnek_cxxflags.append( "-DOCCA_USE_CPU=1"     )
brainnek_cxxflags.append( "-O0"     )

if ( MacOSXBuild ):
  brainnek_cxxflags.append( "-Wno-error=unused-command-line-argument-hard-error-in-future")

brainnek_libraries=[
"lapack",
"blas",
"m",
#"gfortran",
#"mpichcxx",
#"mpich",
#"opa",
"pthread",
#"rt",
"dl",
#"glut",
#"GL",
#"GLU",
#"GLEW",
"OpenCL"
               ]
if ( MacOSXBuild ):
  brainnek_library_dirs=[
  "/Applications/Xcode.app/Contents//Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk//System/Library/Frameworks/OpenCL.framework/"
      ] 
else:
  brainnek_library_dirs=[]

# source filenames
source_files = ["wrapbrainNek.pyx",
                "%s/src/brain3d.cpp"            % (brainnek_dir ),
                "%s/src/brainMaterials.cpp"     % (brainnek_dir ),
                "%s/src/buttonWidget.cpp"       % (brainnek_dir ),
                "%s/src/matrix.cpp"             % (brainnek_dir ),
                "%s/src/mesh3d.cpp"             % (brainnek_dir ),
                "%s/src/scp3d.cpp"              % (brainnek_dir ),
                "%s/src/sem3d.cpp"              % (brainnek_dir ),
                "%s/src/setupAide.cpp"          % (brainnek_dir ),
                "%s/src/sliderWidget.cpp"       % (brainnek_dir ),
                "%s/src/surfaceCutter.cpp"      % (brainnek_dir ),
                "%s/src/visTools.cpp"           % (brainnek_dir ),
                "%s/src/visualizer.cpp"         % (brainnek_dir ),
                "%s/src/vtk3d.cpp"              % (brainnek_dir ),
                "%s/src/vtu3d.cpp"              % (brainnek_dir ),
                "%s/src/widget.cpp"             % (brainnek_dir ),
                "%s/libocca/occa.cpp"           % (brainnek_dir )
               ]

dependency_files = list( source_files) # copy
dependency_files.extend(
               ["wrapbrainNek.pxi",  
               ] )  # header filenames
print source_files
ext_modules_list=[
    Extension("brainNekLibrary", # name of extension
    source_files , # source filenames
    language="c++",              # this causes Cython to create C++ source
    include_dirs= brainnek_include, # include directories
    extra_compile_args = brainnek_cxxflags, # cxx flags
    library_dirs =brainnek_library_dirs, 
    #runtime_library_dirs=brainnek_library_dirs ,
    libraries=brainnek_libraries
    ) # Unix-like specific
]


#Cython will generate and compile the signalmodel.cpp file (from the
#signalmodel.pyx), then it will compile signalmodel.cxx (implementation
#of the SignalModel class) and link both objects files together into
#signalmodel.so, which you can then import in Python using import
#signalmodel (if you forget to link the SignalModel.o, you will get
#missing symbols while importing the library in Python).

setup(
  name = "brainNekInterface",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules_list
)

