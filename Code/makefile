obj_path     = build.tmp
MATLABROOT	:= /opt/MATLAB/R2011a/
MCC=$(MATLABROOT)/bin/mcc

#======================================================================
## $ mpic++ -show
## c++ -g -O2 -g -Wall -O2 -Wl,-Bsymbolic-functions -I/usr/include/mpich2 -L/usr/lib -L/usr/lib -Wl,-rpath -Wl,/usr/lib -lmpichcxx -lmpich -lopa -lpthread -lrt
#
## mpic++ -o main  -DMPI_ENABLED=1 -Ddatafloat=float -DGL_ENABLED=0 -DSCREENSHOT=1 -DVERBOSE=0 -D OS_LINUX=1 -DOCCA_USE_OPENCL=1 -DOCCA_USE_CUDA=0 -DOCCA_USE_ALL=0 -DOCCA_USE_CPU=1 obj/brain3d.o obj/brainMaterials.o obj/buttonWidget.o obj/matrix.o obj/mesh3d.o obj/scp3d.o obj/sem3d.o obj/setupAide.o obj/sliderWidget.o obj/surfaceCutter.o obj/visTools.o obj/visualizer.o obj/vtk3d.o obj/vtu3d.o obj/widget.o libocca/occa.o main.cpp -I./src -I./include -I./libocca -I./include -I/opt/apps/khronos/1.1/ -I/usr/local/cuda/include -m64 -L. -L/usr/local/cuda/lib -llapack -lblas -lm -lgfortran -lrt -ldl -lglut -lGL -lGLU -lGLEW -lOpenCL

.PHONY: $(obj_path)/brainNekLibrary.so tags

$(obj_path)/brainNekLibrary.so:
	/Users/fuentes/MyProjects/Slicer4-SuperBuild-Debug/Cython-0.16/../python-install/bin/python setup.py build_ext -g  -v --build-lib=$(obj_path) --build-temp=$(obj_path)
	#python setup.py build_ext -g  -v --build-lib=$(obj_path) --build-temp=$(obj_path)
	chmod 755 $(obj_path)
	chmod 755 $(obj_path)/brainNekLibrary.so

clean:
	rm -rf $(obj_path)/brainNekLibrary.so ./wrapbrainNek.cpp

BrainNek_SOURCE = /workarea/fuentes/braincode/tym1/

tags:
	ctags -R --langmap=c++:+.tpp --langmap=c++:+.occa --langmap=c++:+.cl $(BrainNek_SOURCE) ./wrapbrainNek.* PyLITTPlan analytic

hello:  hello.m 
	$(MCC) -m $^ -o $@
