# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.
# DAKOTA will execute this script 

# necessary python modules
import sys
import re
import os
import ConfigParser

# numerical support
import numpy

# vis support
import vtk
import vtk.util.numpy_support as vtkNumPy 
print "using vtk version", vtk.vtkVersion.GetVTKVersion()

if( os.getenv("GPUWORKDIR") ) :
  workDirectory   = os.getenv("GPUWORKDIR") 
else:
  workDirectory   = 'optpp_pds/1'
os.system('mkdir -p %s' % workDirectory )

#FIXME global vars
globalconfig = ConfigParser.SafeConfigParser({})
globalconfig.read('./global.ini')
databaseDIR     = globalconfig.get('exec','databaseDIR')
c3dexe          = globalconfig.get('exec','c3dexe')
brainNekDIR     = globalconfig.get('exec','brainNekDIR')
outputDirectory = globalconfig.get('exec','outputDirectory')
MatlabDriver    = globalconfig.getboolean('exec','MatlabDriver')

# FIXME quick hack for 180deg flip 
FIXMEHackTransform = vtk.vtkTransform()
FIXMEHackTransform.RotateY( 180. )

# $ ls database workdir/
# database:
# Patient0002/  Patient0003/  Patient0004/  Patient0005/  Patient0006/  Patient0007/  Patient0008/
# 
# workdir/:
# Patient0002/  Patient0003/  Patient0004/  Patient0005/  Patient0006/  Patient0007/  Patient0008/
# $ ls database/Patient0002/ workdir/Patient0002
# database/Patient0002/:
# 000/  001/  010/  017/  020/  021/  022/
# 
# workdir/Patient0002:
# 000/  001/  010/  017/  020/  021/  022/

caseFunctionTemplate = \
"""
// Prototypes
occaDeviceFunction datafloat laserPower(datafloat time);
occaDeviceFunction datafloat initialTemperature(datafloat x, datafloat y, datafloat z);
occaDeviceFunction datafloat sourceFunction(datafloat x , datafloat y , datafloat z, 
			 datafloat x0, datafloat y0, datafloat z0, 
			 datafloat volumeFraction, 
			 datafloat muA_muTr      , datafloat muEff);
occaDeviceFunction datafloat DirichletTemp(unsigned int bcTag, datafloat time,
			datafloat x, datafloat y, datafloat z);
occaDeviceFunction datafloat RobinCoeff(unsigned int bcTag, datafloat x, datafloat y, datafloat z,
		     datafloat kappa, datafloat h);
occaDeviceFunction datafloat NeumannDeriv(unsigned int bcTag, datafloat time, 
		       datafloat x, datafloat y, datafloat z,
		       datafloat nx, datafloat ny, datafloat nz);
occaDeviceFunction datafloat exactSolution(datafloat x, datafloat y, datafloat z, datafloat time,
			datafloat kappa, datafloat lambda);

/*
 * Compile-time definitions
 *   - bodyTemperature    = ambient body temperature
 *   - coolantTemperature = probe coolant temperature
 *   - laserMaxPower      = reference laser power
 */

/// Initial temperature
/**
 * Boundary conditions will be enforced afterwards
 * @param x
 * @param y
 * @param z
 * @param bodyTemperature
 * @return initial temperature
 */
occaDeviceFunction datafloat initialTemperature(datafloat x, datafloat y, datafloat z) {
  return bodyTemperature;
}

/// Heating at a point due to a region of the laser tip
/**
 * @param x
 * @param y
 * @param z
 * @param x0 x-coordinate of centroid of laser tip region
 * @param y0 y-coordinate of centroid of laser tip region
 * @param z0 z-coordinate of centroid of laser tip region
 * @param volumeFraction volume fraction of laser tip region relative to the 
 *          entire laser tip
 * @param mu_a absorption coefficient of laser light in tissue
 * @param mu_eff effective absorption (\f$\mu_\text{eff}=\sqrt{3\mu_a\mu_{tr}}\f$)
 * @param mu_tr transport coefficient (\f$\mu_{tr}=\mu_a + \mu_s (1-g)\f$)
 * @return contribution of source point to heating function
 */
occaDeviceFunction datafloat sourceFunction(datafloat x , datafloat y , datafloat z, 
			 datafloat x0, datafloat y0, datafloat z0, 
			 datafloat volumeFraction, 
			 datafloat muA_muTr      , datafloat muEff) {
  // Distance between point and source point
  datafloat dist = (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0);
  dist = sqrt(dist);

  // Choose minimum distance to avoid dividing by zero
  if(dist < 1e-6)
    return 0;

  // Return contribution to forcing function
  return 0.75*M_1_PI*muA_muTr*volumeFraction*exp(-muEff*dist)/dist;
}


/// Returns the temperature corresponding to a Dirichlet boundary condition
/**
 * @param bcTag type of boundary condition
 *          - 1 = body temperature Dirichlet boundary condition
 *          - 2 = coolant temperature Dirichlet boundary condition
 * @param x
 * @param y
 * @param z
 * @param time
 * @return Dirichlet boundary condition temperature
 */
occaDeviceFunction datafloat DirichletTemp(unsigned int bcTag, datafloat time, 
			datafloat x, datafloat y, datafloat z) {
  switch(bcTag) {
  case 1:  return bodyTemperature;
  case 2:  return coolantTemperature;
  default: break;
  }
  
  return bodyTemperature;
}

/// Returns the coefficient corresponding to a Robin boundary condition
/**
 * We assume a Robin boundary condition of the form 
 *   \f[\kappa\frac{\partial u}{\partial n}=-\alpha\left(u-u_b\right)\f]
 * @param bcTag type of boundary condition
 *          - 3 = Neumann boundary condition (\f$\alpha=0\f$)
 *          - 4 = Robin condition at probe
 * @param x
 * @param y
 * @param z
 * @param kappa thermal conductivity
 * @param h heat transfer coefficient
 * @return \f$\alpha\f$
 */
occaDeviceFunction datafloat RobinCoeff(unsigned int bcTag, 
		     datafloat x, datafloat y, datafloat z,
		     datafloat kappa, datafloat h) {
  switch(bcTag) {
  case 3:
    return 0;
  case 4:
    return h; // Heat transfer coefficient
  default: return 0;
  }
}		  


/// Returns the derivative corresponding to a Neumann boundary condition
/**
 * Note: not currently used
 * @param bcTag type of boundary condition
 * @param time
 * @param x
 * @param y
 * @param z
 * @param nx x-coordinate of surface normal vector
 * @param ny y-coordinate of surface normal vector
 * @param nz z-coordinate of surface normal vector
 */
occaDeviceFunction datafloat NeumannDeriv(unsigned int bcTag, datafloat time, 
		       datafloat x, datafloat y, datafloat z,
		       datafloat nx, datafloat ny, datafloat nz){
  // Homogeneous Neumann
  return 0;
}

/// Analytic solution
/**
 * Note: an analytic solution is not known for this case. 
 * The numerical solution can be compared with the analytic solution, if it 
 *   is known
 * @param kappa tissue thermal conductivity
 * @param lambda (tissue density)*(tissue specific heat)/dt 
 *                 + (perfusion)*(blood specific heat)
 * @param x
 * @param y
 * @param z
 * @param time
 * @return temperature
 */
occaDeviceFunction datafloat exactSolution(datafloat x, datafloat y, datafloat z, datafloat time,
			datafloat kappa, datafloat lambda){
  return bodyTemperature;
}
"""

setuprcTemplate = \
"""
[THREAD MODEL]
OpenCL

[CASE FILE]
%s/case.%04d.setup

[MESH FILE]
meshes/cooledConformMesh.inp

[MRI FILE]
./mridata.setup

[POLYNOMIAL ORDER]
3

[DT]
0.25

[FINAL TIME]
%f

[PCG TOLERANCE]
1e-6

[GPU PLATFORM]
0

[GPU DEVICE]
%d

[SCREENSHOT OUTPUT]
%s

[SCREENSHOT INTERVAL]
%f
"""

caseFileTemplate = \
"""
# case setup file
# all physical quantities should be in MKS units and degrees Celsius

[FUNCTION FILE]
%s/casefunctions.%04d.occa

[HAS EXACT SOLUTION]
0

# first row is center of probe, second row specifies direction
[PROBE POSITION]
%s

[LASER MAXIMUM POWER]
15

[BODY TEMPERATURE]
%s

[BLOOD TEMPERATURE]
%s

[COOLANT TEMPERATURE]
%s

[BLOOD SPECIFIC HEAT]
%12.5e

[DAMAGE FREQUENCY FACTOR]
1e70

[DAMAGE ACTIVATION ENERGY]
4e5

[GAS CONSTANT]
8.314

[PROBE HEAT TRANSFER COEFFICIENT]
0

[MESH BLOCKS - COMPUTATIONAL DOMAIN]
brain

[MESH BLOCKS - LASER TIP]
laserTip

[MESH BLOCKS - DIRICHLET (BODY TEMPERATURE)]

[MESH BLOCKS - DIRICHLET (COOLANT TEMPERATURE)]

[MESH SIDE SETS - DIRICHLET (BODY TEMPERATURE)]
regionBoundary

[MESH SIDE SETS - DIRICHLET (COOLANT TEMPERATURE)]

[MESH SIDE SETS - NEUMANN]
probeSurface

[MESH SIDE SETS - ROBIN]

[BRAIN TYPES - DIRICHLET (BODY TEMPERATURE)]

[BRAIN MATERIAL DATA FILE]
%s

[BRAIN MATERIAL PROPERTIES FILE]
%s/material_types.%04d.setup

# Currently has material properties of water
[PROBE MATERIAL PROPERTIES]
# Name,   Density, Specific Heat, Conductivity, Absorption, Scattering, Anisotropy
catheter  1.0      4180           0.5985        500         14000       0.88
laserTip  1.0      4180           0.5985        500         14000       0.88
"""
# Convenience Routine
def GetMinJobID(FileNameTemplate):
    MinID      = 1  
    MinObjVal  = 1.e99
    MinDiceVal = 1.e99
    # get a list of all output files in the directory 
    DirectoryLocation = FileNameTemplate.split('/')
    FileTypeID = DirectoryLocation.pop() 
    DirectoryLocation = '/'.join(DirectoryLocation)
    DirectoryOutList = filter(lambda x: len( x.split("%s.out" % FileTypeID) ) == 2 , os.listdir(DirectoryLocation ))
    for dakotaoutfile in DirectoryOutList:
      obj_fn_data = numpy.loadtxt('%s/%s'  % (DirectoryLocation ,dakotaoutfile ) )
      #print '%s/%s'  % (DirectoryLocation, dakotaoutfile), obj_fn_data 
      # FIXME: find the best one, ignore errors
      if( obj_fn_data.size > 1 ):
       if(obj_fn_data[0] < MinObjVal ): 
         MinObjVal  = obj_fn_data[0]
         MinDiceVal = obj_fn_data[1]
         MinID     = int(dakotaoutfile.split(".").pop()) 
    return (MinID,MinObjVal,MinDiceVal )

# Convenience Routine
def DiceTxtFileParse(DiceInputFilename):
  # (1) split on ':' (2)  filter lists > 1 (3) convert to dictionary
  c3doutput = dict(filter( lambda x: len(x) > 1,[line.strip().split(':') for line in open(DiceInputFilename) ] ))
  return float(c3doutput['Dice similarity coefficient'])
  
# Convenience Routine
def WriteVTKOutputFile(vtkImageData,VTKOutputFilename):
    vtkImageDataWriter = vtk.vtkDataSetWriter()
    vtkImageDataWriter.SetFileTypeToBinary()
    print "writing ", VTKOutputFilename 
    vtkImageDataWriter.SetFileName( VTKOutputFilename )
    vtkImageDataWriter.SetInput(vtkImageData)
    vtkImageDataWriter.Update()

# Convenience Routine
def WriteJPGOutputFiles(**visargs):
    print 'opening' , visargs['magnitudefilename'] 
    vtkMagnImageReader = vtk.vtkDataSetReader() 
    vtkMagnImageReader.SetFileName( visargs['magnitudefilename'] )
    vtkMagnImageReader.Update() 
    # display VOI outline
    vtkVOIExtract = vtk.vtkExtractVOI() 
    vtkVOIExtract.SetInput( vtkMagnImageReader.GetOutput() ) 
    vtkVOIExtract.SetVOI( visargs['voi'] ) 
    vtkVOIExtract.Update()
    
    vtkOutline = vtk.vtkOutlineFilter() 
    vtkOutline.SetInput( vtkVOIExtract.GetOutput() ) 
    vtkOutline.Update()
    
    # slight rotation fixes vis bug/error
    AffineTransform = vtk.vtkTransform()
    AffineTransform.RotateX( 1.0 )
    OutlineRegister = vtk.vtkTransformFilter()
    OutlineRegister.SetInput( vtkOutline.GetOutput() )
    OutlineRegister.SetTransform(AffineTransform)
    OutlineRegister.Update()

    # create actor to render VOI
    outlineMapper = vtk.vtkPolyDataMapper(); 
    outlineMapper.SetInput(OutlineRegister.GetOutput());
    outlineActor = vtk.vtkActor(); 
    outlineActor.SetMapper(outlineMapper); 
    outlineActor.GetProperty().SetColor(1,1,1);
    outlineActor.GetProperty().SetLineWidth(2);
    outlineActor.GetProperty().SetRepresentationToWireframe();
    outlineActor.GetProperty().SetRepresentationToPoints();
    outlineActor.GetProperty().SetRepresentationToSurface();

    # Start by creating a black/white lookup table.
    bwLut = vtk.vtkLookupTable()
    bwLut.SetTableRange (0, 300);
    bwLut.SetSaturationRange (0, 0);
    bwLut.SetHueRange (0, 0);
    bwLut.SetValueRange (0, 1);
    bwLut.Build(); #effective built
    # color table
    # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/ImageProcessing/Python/ImageSlicing.py
    hueLut = vtk.vtkLookupTable()
    hueLut.SetNumberOfColors (256)
    #FIXME: adjust here to change color  range
    hueLut.SetRange ( 30.,80.)  
    #hueLut.SetSaturationRange (0.0, 1.0)
    #hueLut.SetValueRange (0.0, 1.0)
    hueLut.SetHueRange (0.667, 0.0)
    hueLut.SetRampToLinear ()
    hueLut.Build()

    # color table
    # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/ImageProcessing/Python/ImageSlicing.py
    arrLut = vtk.vtkLookupTable()
    arrLut.SetNumberOfColors (256)
    #FIXME: adjust here to change color  range
    arrLut.SetRange ( 0., 1.)  
    #arrLut.SetSaturationRange (0.0, 1.0)
    #arrLut.SetValueRange (0.0, 1.0)
    arrLut.SetHueRange (0.667, 0.0)
    arrLut.SetRampToLinear ()
    arrLut.Build()
    # plot mrti, fem, and magn
    for (lookuptable,legendname,sourceimage,outputname) in [  
                           (hueLut,"SEM" ,visargs['roisem']     ,'roisem' ),
                           (hueLut,"MRTI",visargs['roimrti']    ,'roimrti'),
                           (arrLut,"PDAM",visargs['roisemdose'] ,'roisemdose' ),
                           (arrLut,"MDAM",visargs['roimrtidose'],'roimrtidose'),
                           (bwLut ,"Magn",vtkMagnImageReader.GetOutput(),"magn")]:
      # colorbar
      # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
      scalarBar = vtk.vtkScalarBarActor()
      scalarBar.SetTitle(legendname)
      scalarBar.SetNumberOfLabels(4)
      scalarBar.SetLookupTable(lookuptable)

      # mapper
      #mapper = vtk.vtkDataSetMapper()
      mapper = vtk.vtkImageMapToColors()
      mapper.SetInput(  sourceimage )
      # set echo to display
      mapper.SetActiveComponent( 0 )
      mapper.SetLookupTable(lookuptable)
  
      # actor
      actor = vtk.vtkImageActor()
      actor.SetInput(mapper.GetOutput())
       
      # assign actor to the renderer
      ren = vtk.vtkRenderer()
      ren.AddActor(actor)
      ren.AddActor(outlineActor)
      ren.AddActor2D(scalarBar)
      renWin = vtk.vtkRenderWindow()
      renWin.AddRenderer(ren)
      renWin.SetSize(512,512)
      renWin.Render()

      windowToImage = vtk.vtkWindowToImageFilter() 
      windowToImage.SetInput(renWin)
      windowToImage.Update()
      jpgWriter     = vtk.vtkJPEGWriter() 
      jpgWriter.SetFileName( visargs['jpgoutnameformat'] % (outputname) )
      jpgWriter.SetInput(windowToImage.GetOutput())
      jpgWriter.Write()

##################################################################
##################################################################
##################################################################
class ImageDoseHelper:
  """ Class for output of arrhenius dose...  """
  def __init__(self,VOISizeInfo,DeltaT,ImageFileNameExample ):
    print " class constructor called \n\n" 
    # damage paramters
    self.ActivationEnergy     = 3.1e98
    self.FrequencyFactor      = 6.28e5
    self.GasConstant          = 8.314472
    self.BaseTemperature      = 273.0
    self.DeltaT               = DeltaT               

    # open image and extract VOI for a reference
    vtkImageReader = vtk.vtkDataSetReader() 
    vtkImageReader.SetFileName(ImageFileNameExample )
    vtkImageReader.Update() 
    
    # extract voi for QOI
    vtkVOIExtract = vtk.vtkExtractVOI() 
    vtkVOIExtract.SetInput( vtkImageReader.GetOutput() ) 
    vtkVOIExtract.SetVOI( VOISizeInfo ) 
    vtkVOIExtract.Update()
    # VOI Origin should be at the lower bound
    VOIBounds = vtkVOIExtract.GetOutput().GetBounds()
    self.origin               = (VOIBounds[0],VOIBounds[2],VOIBounds[4])
    self.spacing              = vtkVOIExtract.GetOutput().GetSpacing()

    # initialize dose map
    self.dimensions = [(VOISizeInfo[1] - VOISizeInfo[0]+1) , 
                       (VOISizeInfo[3] - VOISizeInfo[2]+1) ,
                       (VOISizeInfo[5] - VOISizeInfo[4]+1) ,1 ]
    numpyimagesize = self.dimensions[0]*self.dimensions[1]*self.dimensions[2]
    # store as double precision
    self.PredictedDamage = numpy.zeros( numpyimagesize, dtype=numpy.float ) 

  def UpdateDoseMap(self,NumpyTemperatureData ):
    """ update dose map with temperature """
    #  input should be temperature in degC
    #  convert to Kelvin (using double precision)
    TemperatureKelvin = NumpyTemperatureData.astype(numpy.float) + self.BaseTemperature

    #  A exp ( - E_a/ (RT)  ) \Delta t
    self.PredictedDamage = self.PredictedDamage + self.ActivationEnergy * numpy.exp(- self.FrequencyFactor/self.GasConstant * numpy.reciprocal( TemperatureKelvin )) * self.DeltaT;

    # return vtk format for write (as single precision)
    return self.ConvertNumpyVTKImage(self.PredictedDamage.astype(numpy.float32))

  # write a numpy data to disk in vtk format
  def ConvertNumpyVTKImage(self,NumpyImageData):
    # Create initial image
    dim = self.dimensions
    # imports raw data and stores it.
    dataImporter = vtk.vtkImageImport()
    #  indexing is painful.... reshape to dimensions and transpose 2d dimensions only
    ## try:
    ##   numpytmp = NumpyImageData.reshape(self.dimensions[0:3],order='F').transpose(1,0,2)
    ## except:
    ##   numpytmp = NumpyImageData.reshape(self.dimensions[0:2],order='F').transpose(1,0)
    data_string = NumpyImageData.tostring()
    # array is converted to a string of chars and imported.
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToFloat()
    # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
    # must be told this is the case.
    dataImporter.SetNumberOfScalarComponents(dim[3])
    # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
    # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
    # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
    # VTK complains if not both are used.
    dataImporter.SetDataExtent( 0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)
    dataImporter.SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)
    dataImporter.SetDataSpacing( self.spacing )
    dataImporter.SetDataOrigin(  self.origin )
    dataImporter.Update()
    dataImporter.SetScalarArrayName("arrayname")
    return dataImporter.GetOutput()


def ForwardSolve(**kwargs):
  ObjectiveFunction = 0.0
  # Debugging flags
  DebugObjective = True
  DebugObjective = False

  # initialize brainNek
  # CYTHON AND VTK NEED TO BUILD with the same PYTHON INCLUDE and LIB 
  import brainNekLibrary
  # setuprc file
  outputSetupRCFile = '%s/setuprc.%04d' % (workDirectory,kwargs['fileID'])
  setup = brainNekLibrary.PySetupAide(outputSetupRCFile )
  brainNek = brainNekLibrary.PyBrain3d(setup);

  # setup vtkUnstructuredGrid
  hexahedronGrid   = vtk.vtkUnstructuredGrid()
  numPoints = brainNek.GetNumberOfNodes( ) 
  numElems  = brainNek.GetNumberOfElements( ) 
  # initialize nodes and connectivity
  numHexPts = 8 
  bNekNodes         = numpy.zeros(numPoints * 3,dtype=numpy.float32)
  bNekConnectivity  = numpy.zeros(numElems  * (numHexPts +1),dtype=numpy.int32)
  print "setting up hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # get nodes and connectivity from brainnek
  brainNek.GetNodes(   bNekNodes)       ;
  brainNek.GetElements(bNekConnectivity);
  # reshape for convenience
  bNekNodes        = bNekNodes.reshape(      numPoints , 3)

  # TODO : check if deepcopy needed
  DeepCopy = 1

  #hexahedronGrid.DebugOn()
  # setup points
  hexahedronPoints = vtk.vtkPoints()
  vtkNodeArray = vtkNumPy.numpy_to_vtk( bNekNodes, DeepCopy)
  hexahedronPoints.SetData(vtkNodeArray)
  hexahedronGrid.SetPoints(hexahedronPoints);

  # setup elements
  aHexahedron = vtk.vtkHexahedron()
  HexCellType = aHexahedron.GetCellType()
  vtkTypeArray     = vtkNumPy.numpy_to_vtk( HexCellType * numpy.ones(  numElems) ,DeepCopy,vtk.VTK_UNSIGNED_CHAR) 
  #TODO: off by 1 indexing from npts, ie
  #TODO: note vtkIdType vtkCellArray::InsertNextCell(vtkIdList *pts) 
  #TODO:    this->InsertLocation += npts + 1;   (line 264)
  vtkLocationArray = vtkNumPy.numpy_to_vtk( numpy.arange(0,numElems*(numHexPts+1),(numHexPts+1)) ,DeepCopy,vtk.VTK_ID_TYPE) 
  vtkCells = vtk.vtkCellArray()
  vtkElemArray     = vtkNumPy.numpy_to_vtk( bNekConnectivity  , DeepCopy,vtk.VTK_ID_TYPE)
  vtkCells.SetCells(numElems,vtkElemArray)
  hexahedronGrid.SetCells(vtkTypeArray,vtkLocationArray,vtkCells) 
  print "done setting hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # setup solution
  bNekSoln = numpy.zeros(numPoints,dtype=numpy.float32)
  brainNek.getHostTemperature(bNekSoln )
  vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
  vtkScalarArray.SetName("bioheat") 
  hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);

  ## # dbg 
  ## brainNek.screenshot( 0.0 )

  # get registration parameters
  variableDictionary = kwargs['cv']

  # if ( DebugObjective ):
  #    vtkSEMReader = vtk.vtkXMLUnstructuredGridReader()
  #    SEMDataDirectory = outputDirectory % kwargs['UID']
  #    SEMtimeID = 0 
  #    vtufileName = "%s/%d.vtu" % (SEMDataDirectory,SEMtimeID)
  #    print "reading ", vtufileName 
  #    vtkSEMReader.SetFileName( vtufileName )
  #    vtkSEMReader.SetPointArrayStatus("Temperature",1)
  #    vtkSEMReader.Update()
  #    fem_point_data= vtkSEMReader.GetOutput().GetPointData() 
  #    tmparray = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Temperature')) 

  # loop over time steps
  tstep = 0
  currentTime = 0.0

  # setup MRTI data read
  MRTItimeID  = 0
  MRTIInterval = 5.0

  # setup screen shot interval 
  screenshotNum = 1;
  screenshotTol = 1e-10;

  ## loop over time
  PowerLambdaFunction = kwargs['lambdacode']
  while( brainNek.timeStep(tstep * brainNek.dt(), PowerLambdaFunction(currentTime )  ) ) :
    tstep = tstep + 1
    currentTime = tstep * brainNek.dt()

    if(currentTime+screenshotTol >= MRTItimeID * MRTIInterval):
      
      #print mrti_array
      #print type(mrti_array)

      # get brainNek solution 
      brainNek.getHostTemperature( bNekSoln )
      vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
      vtkScalarArray.SetName("bioheat") 
      hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);
      hexahedronGrid.Update()

      #print fem_array 
      #print type(fem_array )

      # FIXME  should this be different ?  
      SEMDataDirectory = outputDirectory % kwargs['UID']

      # write output
      if ( DebugObjective ):
        vtkSEMWriter = vtk.vtkXMLUnstructuredGridWriter()
        semfileName = "%s/semtransform.%04d.vtu" % (SEMDataDirectory,MRTItimeID)
        print "writing ", semfileName 
        vtkSEMWriter.SetFileName( semfileName )
        vtkSEMWriter.SetInput(hexahedronGrid)
        #vtkSEMWriter.SetDataModeToAscii()
        vtkSEMWriter.Update()

      # update counter
      MRTItimeID = MRTItimeID + 1;

  # read landmarks to transform
  SlicerLMReader   = vtk.vtkPolyDataReader()
  ParaviewLMReader = vtk.vtkPolyDataReader()
  ParaviewLMReader.SetFileName(  kwargs['target_landmarks']  )
  SlicerLMReader.SetFileName(    kwargs['transformlandmarks'])
  ParaviewLMReader.Update()
  SlicerLMReader.Update(  )
  # initialize transform
  LaserLineTransform = vtk.vtkLandmarkTransform()
  LaserLineTransform.SetSourceLandmarks(ParaviewLMReader.GetOutput().GetPoints())
  LaserLineTransform.SetTargetLandmarks(SlicerLMReader.GetOutput().GetPoints())
  LaserLineTransform.SetModeToRigidBody()
  LaserLineTransform.Update()
  print LaserLineTransform.GetMatrix()

  # apply slicer transform
  rastransformFilter = vtk.vtkTransformFilter()
  rastransformFilter.SetInput( hexahedronGrid ) 
  rastransformFilter.SetTransform(LaserLineTransform) 
  rastransformFilter.Update()

  # scale to millimeter for vis
  ScaleAffineTransform = vtk.vtkTransform()
  ScaleAffineTransform.Translate([ 0.0,0.0,0.0])
  ScaleAffineTransform.RotateZ( 0.0  )
  ScaleAffineTransform.RotateY( 0.0  )
  ScaleAffineTransform.RotateX( 0.0  )
  ScaleAffineTransform.Scale([1000.,1000.,1000.])

  scaletransformFilter = vtk.vtkTransformFilter()
  scaletransformFilter.SetInput( rastransformFilter.GetOutput() ) 
  scaletransformFilter.SetTransform(ScaleAffineTransform ) 
  scaletransformFilter.Update()

  # setup contour filter
  vtkContour = vtk.vtkContourFilter()
  vtkContour.SetInput( scaletransformFilter.GetOutput() )
  # TODO: not sure why this works...
  # set the array to process at the temperature == bioheat 
  vtkContour.SetInputArrayToProcess(0,0,0,0,'bioheat')
  ## contourValuesList  = eval(newIni.get('exec','contours'))
  contourValuesList  = [57. , 62.]
  vtkContour.SetNumberOfContours( len(contourValuesList ) )
  print "plotting array:", vtkContour.GetArrayComponent( )
  for idContour,contourValue in enumerate(contourValuesList):
     print "plotting contour:",idContour,contourValue
     vtkContour.SetValue( idContour,contourValue )
  vtkContour.Update( )

  # write stl file
  stlWriter = vtk.vtkSTLWriter()
  stlWriter.SetInput(vtkContour.GetOutput( ))
  stlWriter.SetFileName("fem.stl")
  stlWriter.SetFileTypeToBinary()
  stlWriter.Write()
  # write support file to signal full stl file is written
  #  Slicer module will wait for this file to be written
  #  before trying to open stl file to avoid incomplete reads
  with open('./fem.finish', 'w') as signalFile:
    signalFile.write("the stl file has been written\n")
  return ObjectiveFunction 
# end def ForwardSolve:
##################################################################
def ComputeObjective(**kwargs):
  ObjectiveFunction = 0.0
  dicevalue = 0.0
  # Debugging flags
  DebugObjective = False
  DebugObjective = True

  # initialize brainNek
  # CYTHON AND VTK NEED TO BUILD with the same PYTHON INCLUDE and LIB 
  import brainNekLibrary
  # setuprc file
  outputSetupRCFile = '%s/setuprc.%04d' % (workDirectory,kwargs['fileID'])
  setup = brainNekLibrary.PySetupAide(outputSetupRCFile )
  brainNek = brainNekLibrary.PyBrain3d(setup);


  # setup vtkUnstructuredGrid
  hexahedronGrid   = vtk.vtkUnstructuredGrid()
  numPoints = brainNek.GetNumberOfNodes( ) 
  numElems  = brainNek.GetNumberOfElements( ) 
  # initialize nodes and connectivity
  numHexPts = 8 
  bNekNodes         = numpy.zeros(numPoints * 3,dtype=numpy.float32)
  bNekConnectivity  = numpy.zeros(numElems  * (numHexPts +1),dtype=numpy.int32)
  print "setting up hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # get nodes and connectivity from brainnek
  brainNek.GetNodes(   bNekNodes)       ;
  brainNek.GetElements(bNekConnectivity);
  # reshape for convenience
  bNekNodes        = bNekNodes.reshape(      numPoints , 3)

  ## # setup elements
  ## bNekConnectivityreshape = bNekConnectivity.reshape(numElems , numHexPts +1)
  ## for ielem in range(numElems ):
  ##   aHexahedron = vtk.vtkHexahedron()
  ##   # print 'number of nodes %d ' % bNekConnectivityreshape[ielem][0]
  ##   aHexahedron.GetPointIds().SetId(0,bNekConnectivityreshape[ielem][1])
  ##   aHexahedron.GetPointIds().SetId(1,bNekConnectivityreshape[ielem][2])
  ##   aHexahedron.GetPointIds().SetId(2,bNekConnectivityreshape[ielem][3])
  ##   aHexahedron.GetPointIds().SetId(3,bNekConnectivityreshape[ielem][4])
  ##   aHexahedron.GetPointIds().SetId(4,bNekConnectivityreshape[ielem][5])
  ##   aHexahedron.GetPointIds().SetId(5,bNekConnectivityreshape[ielem][6])
  ##   aHexahedron.GetPointIds().SetId(6,bNekConnectivityreshape[ielem][7])
  ##   aHexahedron.GetPointIds().SetId(7,bNekConnectivityreshape[ielem][8])
  ##   hexahedronGrid.InsertNextCell(aHexahedron.GetCellType(),
  ##                                 aHexahedron.GetPointIds())

  # TODO : check if deepcopy needed
  DeepCopy = 1

  #hexahedronGrid.DebugOn()
  # setup points
  hexahedronPoints = vtk.vtkPoints()
  vtkNodeArray = vtkNumPy.numpy_to_vtk( bNekNodes, DeepCopy)
  hexahedronPoints.SetData(vtkNodeArray)
  hexahedronGrid.SetPoints(hexahedronPoints);

  # setup elements
  aHexahedron = vtk.vtkHexahedron()
  HexCellType = aHexahedron.GetCellType()
  vtkTypeArray     = vtkNumPy.numpy_to_vtk( HexCellType * numpy.ones(  numElems) ,DeepCopy,vtk.VTK_UNSIGNED_CHAR) 
  #TODO: off by 1 indexing from npts, ie
  #TODO: note vtkIdType vtkCellArray::InsertNextCell(vtkIdList *pts) 
  #TODO:    this->InsertLocation += npts + 1;   (line 264)
  vtkLocationArray = vtkNumPy.numpy_to_vtk( numpy.arange(0,numElems*(numHexPts+1),(numHexPts+1)) ,DeepCopy,vtk.VTK_ID_TYPE) 
  vtkCells = vtk.vtkCellArray()
  vtkElemArray     = vtkNumPy.numpy_to_vtk( bNekConnectivity  , DeepCopy,vtk.VTK_ID_TYPE)
  vtkCells.SetCells(numElems,vtkElemArray)
  hexahedronGrid.SetCells(vtkTypeArray,vtkLocationArray,vtkCells) 
  print "done setting hex mesh with %d nodes %d elem"  % (numPoints,numElems)

  # setup solution
  bNekSoln = numpy.zeros(numPoints,dtype=numpy.float32)
  brainNek.getHostTemperature(bNekSoln )
  vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
  vtkScalarArray.SetName("bioheat") 
  hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);

  MonteCarloSource = True
  MonteCarloSource = False
  if ( MonteCarloSource ):
    # Read In Fluence Source
    vtkForcingImageReader = vtk.vtkDataSetReader() 
    vtkForcingImageReader.SetFileName('./MC_PtSource.0000.vtk')
    vtkForcingImageReader.Update() 
    # Project Fluence Source to gll nodes
    print 'resampling fluence' 
    vtkForcingResample = vtk.vtkCompositeDataProbeFilter()
    vtkForcingResample.SetInput( hexahedronGrid )
    vtkForcingResample.SetSource( vtkForcingImageReader.GetOutput() ) 
    vtkForcingResample.Update()
    resampledForcingMesh = vtkForcingResample.GetOutput() 
    # test registration
    if ( DebugObjective ):
       # compare to old forcing
       bNekForcing = numpy.zeros(numPoints,dtype=numpy.float32)
       brainNek.getHostForcing( bNekForcing )
       vtkForcingArray = vtkNumPy.numpy_to_vtk( bNekForcing , DeepCopy) 
       vtkForcingArray.SetName("oldforcing") 
       # FIXME should be able to write all arrays to single mesh w/o copy ??
       oldForcingCopy = vtk.vtkUnstructuredGrid()
       oldForcingCopy.DeepCopy(resampledForcingMesh)
       oldForcingCopy.GetPointData().AddArray(vtkForcingArray);

       vtkDbgMeshWriter = vtk.vtkDataSetWriter() 
       vtkDbgMeshWriter.SetFileName('oldforcing.vtk')
       vtkDbgMeshWriter.SetInput( oldForcingCopy )
       vtkDbgMeshWriter.Update() 

    # Memory Copy projected solution
    forcing_point_data= resampledForcingMesh.GetPointData() 
    forcing_array = vtkNumPy.vtk_to_numpy(forcing_point_data.GetArray('scalars')) 
    brainNek.setDeviceForcing( forcing_array )
  
  ## # dbg 
  ## brainNek.screenshot( 0.0 )

  # get registration parameters
  variableDictionary = kwargs['cv']

  # register the SEM data to MRTI
  AffineTransform = vtk.vtkTransform()
  AffineTransform.Translate([ 
    float(variableDictionary['x_displace']),
    float(variableDictionary['y_displace']),
    float(variableDictionary['z_displace'])
                            ])
  # FIXME  notice that order of operations is IMPORTANT
  # FIXME   translation followed by rotation will give different results
  # FIXME   than rotation followed by translation
  # FIXME  Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  AffineTransform.RotateZ( float(variableDictionary['z_rotate'  ] ) ) 
  AffineTransform.RotateY( float(variableDictionary['y_rotate'  ] ) )
  AffineTransform.RotateX( float(variableDictionary['x_rotate'  ] ) )
  AffineTransform.Scale([1.e0,1.e0,1.e0])

  # FIXME  should this be different ?  
  SEMDataDirectory = outputDirectory % kwargs['UID']

  ## vtkSEMReader = vtk.vtkXMLUnstructuredGridReader()
  ## SEMDataDirectory = outputDirectory % kwargs['UID']
  ## SEMtimeID = 0 
  ## vtufileName = "%s/%d.vtu" % (SEMDataDirectory,SEMtimeID)
  ## print "reading ", vtufileName 
  ## vtkSEMReader.SetFileName( vtufileName )
  ## vtkSEMReader.SetPointArrayStatus("Temperature",1)
  ## vtkSEMReader.Update()
  ## fem_point_data= vtkSEMReader.GetOutput().GetPointData() 
  ## tmparray = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Temperature')) 

  print " brainNek deltat:", brainNek.dt()

  # setup MRTI data read
  MRTIInterval = fem_params['mrtideltat'] 
  MRTItimeID   = fem_params['timeinterval'][0]

  # initialize image dose
  semDose  = ImageDoseHelper(  kwargs['voi'], MRTIInterval ,'%s/temperature.0001.vtk' % (kwargs['mrti']))
  mrtiDose = ImageDoseHelper(  kwargs['voi'], MRTIInterval ,'%s/temperature.0001.vtk' % (kwargs['mrti']))

  # setup screen shot interval 
  screenshotNum = 1;
  screenshotTol = 1e-10;
  screenshotInterval = MRTIInterval ;

  # initialize temperature field with MRTI for cooling optimization
  if(kwargs['opttype'] == 'cooling'): 
    # load mrti for initial condition 
    mrtifilename = '%s/temperature.%04d.vtk' % ( kwargs['mrti'], MRTItimeID ) 
    print 'initial condition opening' , mrtifilename 
    vtkICImageReader = vtk.vtkDataSetReader() 
    vtkICImageReader.SetFileName(mrtifilename )
    vtkICImageReader.Update() 
    vtkICVOIExtract = vtk.vtkExtractVOI() 
    vtkICVOIExtract.SetInput( vtkICImageReader.GetOutput() ) 
    vtkICVOIExtract.SetVOI( kwargs['voi'] ) 
    vtkICVOIExtract.Update()
    # blur out of plane for physical temperature field
    imagevoiextents = list(vtkICVOIExtract.GetOutput().GetExtent())
    imagevoiextents [4] = imagevoiextents [4] - 1 
    imagevoiextents [5] = imagevoiextents [5] + 1 
    # NOTE to keep the same MRTI values in plane
    # NOTE   1. mirror pad out of plane in one pixel
    # NOTE   2. constant pad out of plane in one pixel
    # NOTE   3. gauss blur with 1 pixel std dev
    vtkMirrorPad = vtk.vtkImageMirrorPad() 
    vtkMirrorPad.SetOutputWholeExtent( imagevoiextents ) 
    vtkMirrorPad.SetInput( vtkICVOIExtract.GetOutput() ) 
    imagevoiextents [4] = imagevoiextents [4] - 1 
    imagevoiextents [5] = imagevoiextents [5] + 1 
    vtkImagePad = vtk.vtkImageConstantPad() 
    vtkImagePad.SetOutputWholeExtent( imagevoiextents ) 
    vtkImagePad.SetInput( vtkMirrorPad.GetOutput() ) 
    body_temp = float(kwargs['cv']['body_temp']) 
    vtkImagePad.SetConstant( body_temp ) 
    GaussSmooth = vtk.vtkImageGaussianSmooth()
    GaussSmooth.SetInput( vtkImagePad.GetOutput() ) 
    GaussSmooth.SetStandardDeviations( 1.e-6,1.e-6,1 ) 
    #GaussSmooth.SetRadiusFactors( 1,1,1.5 )
    # register and resample the MRTI onto the SEM mesh
    ICSEMRegister = vtk.vtkTransformFilter()
    ICSEMRegister.SetInput( hexahedronGrid )
    ICSEMRegister.SetTransform(AffineTransform)
    ICSEMRegister.Update()
    vtkICResample = vtk.vtkProbeFilter()
    vtkICResample.SetSource( GaussSmooth.GetOutput() )
    vtkICResample.SetInput( ICSEMRegister.GetOutput() ) 
    vtkICResample.Update()
    mrti_ic_point_data= vtkICResample.GetUnstructuredGridOutput().GetPointData() 
    mrti_ic_array = vtkNumPy.vtk_to_numpy(mrti_ic_point_data.GetArray('image_data')) 
    # threshold by body temp
    mrti_ic_array[ mrti_ic_array < body_temp ] = body_temp  
    brainNek.setDeviceTemperature( mrti_ic_array )
    # write output
    if ( DebugObjective ):
      # check gauss image
      WriteVTKOutputFile ( GaussSmooth.GetOutput()   ,"%s/gauss.%s.%04d.vtk"   % (SEMDataDirectory,kwargs['opttype'],MRTItimeID))

      vtkSEMWriter = vtk.vtkXMLUnstructuredGridWriter()
      semfileName = "%s/semtransform.%04d.vtu" % (SEMDataDirectory,MRTItimeID)
      print "writing ", semfileName 
      vtkSEMWriter.SetFileName( semfileName )
      #vtkSEMWriter.SetInput(hexahedronGrid )
      vtkSEMWriter.SetInput(vtkICResample.GetUnstructuredGridOutput())
      #vtkSEMWriter.SetDataModeToAscii()
      vtkSEMWriter.Update()

      # verify temperature on  brainNek data structures 
      brainNek.getHostTemperature( bNekSoln )
      vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
      vtkScalarArray.SetName("bioheat") 
      hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);
      verifSEMRegister = vtk.vtkTransformFilter()
      verifSEMRegister.SetInput( hexahedronGrid )
      verifSEMRegister.SetTransform(AffineTransform)
      verifSEMRegister.Update()

      verifSEMWriter = vtk.vtkXMLUnstructuredGridWriter()
      semfileName = "%s/verifysemtransform.%04d.vtu" % (SEMDataDirectory,MRTItimeID)
      print "writing ", semfileName 
      verifSEMWriter.SetFileName( semfileName )
      verifSEMWriter.SetInput(verifSEMRegister.GetUnstructuredGridOutput())
      verifSEMWriter.Update()

  # debugging info
  brainNek.PrintSelf()

  ## loop over time
  currentTime = fem_params['initialtime'] 
  ## FIXME timing errors
  PowerLambdaFunction = kwargs['lambdacode']
  for MRTItimeID in range(fem_params['timeinterval'][0]+1,fem_params['timeinterval'][1]):

    while( currentTime  < (MRTItimeID +1)*MRTIInterval ) :
      currentTime  = currentTime + brainNek.dt()
      currentPower = PowerLambdaFunction(currentTime   )
      brainNek.heatStep( currentTime,currentPower  )
      ## if(currentTime+screenshotTol >= screenshotNum*screenshotInterval):
      ##    brainNek.getHostTemperature( bNekSoln )
      ##    screenshotNum = screenshotNum + 1;
      ##    print "get host data",bNekSoln 

    # load image 
    mrtifilename = '%s/temperature.%04d.vtk' % (kwargs['mrti'], MRTItimeID) 
    if (os.path.isfile(mrtifilename ) ):
      print 'opening' , mrtifilename , currentTime
    else:
      print '#####NOT FOUND' , mrtifilename 
      print '#####USING DEFAULT at time 0' 
      mrtifilename = '%s/temperature.%04d.vtk' % (kwargs['mrti'], 0) 
    vtkImageReader = vtk.vtkDataSetReader() 
    vtkImageReader.SetFileName(mrtifilename )
    vtkImageReader.Update() 
    ## image_cells = vtkImageReader.GetOutput().GetPointData() 
    ## data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars')) 
    
    # extract voi for QOI
    vtkVOIExtract = vtk.vtkExtractVOI() 
    vtkVOIExtract.SetInput( vtkImageReader.GetOutput() ) 
    vtkVOIExtract.SetVOI( kwargs['voi'] ) 
    vtkVOIExtract.Update()
    mrti_point_data= vtkVOIExtract.GetOutput().GetPointData() 
    mrti_array = vtkNumPy.vtk_to_numpy(mrti_point_data.GetArray('image_data')) 
    # update dose
    vtkmrtiDose = mrtiDose.UpdateDoseMap(mrti_array)
    # x = vtkNumPy.vtk_to_numpy(vtkmrtiDose.GetPointData().GetArray('scalars')) 
    
    #print mrti_array
    #print type(mrti_array)

    # get brainNek solution 
    brainNek.getHostTemperature( bNekSoln )
    vtkScalarArray = vtkNumPy.numpy_to_vtk( bNekSoln, DeepCopy) 
    vtkScalarArray.SetName("bioheat") 
    hexahedronGrid.GetPointData().SetScalars(vtkScalarArray);
    hexahedronGrid.Update()

    # project SEM onto MRTI for comparison
    print 'resampling' 
    FixmeHackSEMRegister = vtk.vtkTransformFilter()
    FixmeHackSEMRegister.SetInput( hexahedronGrid )
    FixmeHackSEMRegister.SetTransform(FIXMEHackTransform)
    FixmeHackSEMRegister.Update()
    SEMRegister = vtk.vtkTransformFilter()
    SEMRegister.SetInput( FixmeHackSEMRegister.GetOutput() )
    SEMRegister.SetTransform(AffineTransform)
    SEMRegister.Update()
    vtkResample = vtk.vtkCompositeDataProbeFilter()
    vtkResample.SetSource( SEMRegister.GetOutput() )
    vtkResample.SetInput( vtkVOIExtract.GetOutput() ) 
    vtkResample.Update()

    fem_point_data= vtkResample.GetOutput().GetPointData() 
    fem_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('bioheat')) 
    # update dose
    vtksemDose  = semDose.UpdateDoseMap(  fem_array)
    print 'resampled' 
    #print fem_array 
    #print type(fem_array )


    # write output
    if ( False ):
      vtkSEMWriter = vtk.vtkXMLUnstructuredGridWriter()
      semfileName = "%s/semtransform.%04d.vtu" % (SEMDataDirectory,MRTItimeID)
      print "writing ", semfileName
      vtkSEMWriter.SetFileName( semfileName )
      vtkSEMWriter.SetInput(SEMRegister.GetOutput())
      #vtkSEMWriter.SetDataModeToAscii()
      vtkSEMWriter.Update()

    # write output
    # FIXME auto read ??
    if ( DebugObjective ):
       # write temperature
       WriteVTKOutputFile ( vtkResample.GetOutput()   ,"%s/roisem.%s.%04d.vtk"   % (SEMDataDirectory,kwargs['opttype'],MRTItimeID))
       WriteVTKOutputFile ( vtkVOIExtract.GetOutput() ,"%s/roimrti.%s.%04d.vtk"  % (SEMDataDirectory,kwargs['opttype'],MRTItimeID))
       # write dose
       WriteVTKOutputFile ( vtksemDose  ,"%s/roisemdose.%s.%04d.vtk"   % (SEMDataDirectory,kwargs['opttype'],MRTItimeID))
       WriteVTKOutputFile ( vtkmrtiDose ,"%s/roimrtidose.%s.%04d.vtk"  % (SEMDataDirectory,kwargs['opttype'],MRTItimeID))

       # write dice coefficient
       dicefilename = "%s/dice.%s.%04d.txt" % ( SEMDataDirectory,kwargs['opttype'],MRTItimeID)
       dicecmd = "%s -verbose %s/roisemdose.%s.%04d.vtk -thresh 1 inf 1 0 -type uchar -as SEM %s/roimrtidose.%s.%04d.vtk -thresh 1 inf 1 0 -type uchar -push SEM -overlap 1 > %s  2>&1" % (c3dexe,SEMDataDirectory,kwargs['opttype'],MRTItimeID,SEMDataDirectory,kwargs['opttype'],MRTItimeID,dicefilename)
       print dicecmd, dicefilename 
       os.system(dicecmd)
       if (  MRTItimeID == fem_params['maxheatid'] ):
         dicevalue = DiceTxtFileParse(dicefilename)
       ##if (MRTItimeID > 20):
       ##  raise 

    # Write JPG's for tex
    if ( kwargs['VisualizeOutput'] and MRTItimeID == fem_params['maxheatid'] ):
    #if ( kwargs['VisualizeOutput'] ):
       VisDictionary = {'voi'    :  kwargs['voi'] ,
                        'roisem' : vtkResample.GetOutput()   ,
                        'roimrti': vtkVOIExtract.GetOutput() ,
                    'roisemdose' : vtksemDose  ,
                    'roimrtidose': vtkmrtiDose ,
              'magnitudefilename':'%s/magnitude.%04d.vtk' % (kwargs['mrti'], MRTItimeID) ,
               'jpgoutnameformat':"%s/%%s%s%04d.jpg"  % (SEMDataDirectory,kwargs['opttype'],MRTItimeID)
                       }
       WriteJPGOutputFiles(**VisDictionary)

    # accumulate objective function
    diff =  numpy.abs(mrti_array-fem_array)
    diffsq =  diff**2
    ObjectiveFunction = ObjectiveFunction + diff.sum()

  return (ObjectiveFunction ,dicevalue)
# end def ComputeObjective:
##################################################################
def brainNekWrapper(**kwargs):
  """
  call brainNek code 
  """
  # occa case file
  outputOccaCaseFile = '%s/casefunctions.%04d.occa' % (workDirectory,kwargs['fileID'])
  print 'writing', outputOccaCaseFile 
  with file(outputOccaCaseFile, 'w') as occaCaseFileName: occaCaseFileName.write(caseFunctionTemplate )

  # setuprc file
  outputSetupRCFile = '%s/setuprc.%04d' % (workDirectory,kwargs['fileID'])
  print 'writing', outputSetupRCFile 
  fileHandle = file(outputSetupRCFile ,'w')
  semfinaltime = kwargs['finaltime']
  # make sure write directory exists
  os.system('mkdir -p %s' % outputDirectory % kwargs['UID'] )
  GPUDeviceID = int(workDirectory.split('/').pop())
  fileHandle.write(setuprcTemplate % (workDirectory,kwargs['fileID'] ,semfinaltime ,GPUDeviceID  ,  outputDirectory % kwargs['UID'] ,semfinaltime ) )
  fileHandle.flush(); fileHandle.close()

  # get variables
  variableDictionary = kwargs['cv']
  rho    = variableDictionary['rho'    ]   
  c_p    = variableDictionary['c_p'    ]   
  k_0    = variableDictionary['k_0'    ]   
  w_0    = variableDictionary['w_0'    ]   
  mu_a   = variableDictionary['mu_a'   ]   
  mu_s   = variableDictionary['mu_s'   ]   
  anfact = variableDictionary['anfact' ]   

  # materials
  outputMaterialFile = '%s/material_types.%04d.setup' % (workDirectory,kwargs['fileID'])
  print 'writing', outputMaterialFile 
  fileHandle = file(outputMaterialFile   ,'w')
  fileHandle.write('[MATERIAL PROPERTIES]\n'  )
  fileHandle.write('# Name,      Type index, Density, Specific Heat, Conductivity, Perfusion, Absorption, Scattering, Anisotropy\n'  )
  fileHandle.write('Brain     0           %12.5f     %12.5f           %12.5f        %12.5f     %12.5f      %12.5f      %12.5f \n' % ( rho, c_p, k_0, w_0, mu_a, mu_s, anfact ))
  fileHandle.write('Tumor     1           %12.5f     %12.5f           %12.5f        %12.5f     %12.5f      %12.5f      %12.5f \n' % ( rho, c_p, k_0, w_0, mu_a, mu_s, anfact ))
  fileHandle.write('CSF      25           %12.5f     %12.5f           %12.5f        %12.5f     %12.5f      %12.5f      %12.5f \n' % ( rho, c_p, k_0, 100.*w_0, mu_a, mu_s, anfact ))
  fileHandle.flush(); fileHandle.close()

  # case file
  outputCaseFile = '%s/case.%04d.setup' % (workDirectory,kwargs['fileID'])
  print 'writing', outputCaseFile 
  with file(outputCaseFile , 'w') as fileHandle: fileHandle.write(caseFileTemplate % (workDirectory,kwargs['fileID'],kwargs['target_landmarks'],variableDictionary['body_temp'],variableDictionary['body_temp'],variableDictionary['probe_init'],variableDictionary['c_blood'],kwargs['segment_file'],workDirectory,kwargs['fileID'])  )

  ## # build command to run brainNek
  ## brainNekCommand = "%s/main %s -heattransfercoefficient %s -coolanttemperature  %s > %s/run.%04d.log 2>&1 " % (brainNekDIR , outputSetupRCFile ,variableDictionary['robin_coeff'  ], variableDictionary['probe_init'   ], workDirectory ,kwargs['fileID'])

  ## # system call to run brain code
  ## print brainNekCommand 
  ## os.system(brainNekCommand )
# end def brainNekWrapper:
##################################################################
def ParseInput(paramfilename,VisualizeOutput):
  # ----------------------------
  # Parse DAKOTA parameters file
  # ----------------------------
  
  # setup regular expressions for parameter/label matching
  e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
  f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
  i = '-?\\d+'                                         # integer
  value = e+'|'+f+'|'+i                                # numeric field
  tag = '\\w+(?::\\w+)*'                               # text tag field
  
  # regular expression for aprepro parameters format
  aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
  # regular expression for standard parameters format
  standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')
  
  # open DAKOTA parameters file for reading
  paramsfile = open(paramfilename, 'r')

  fileID = int(paramfilename.split(".").pop())
  #fileID = int(os.getcwd().split(".").pop())
  
  # extract the parameters from the file and store in a dictionary
  paramsdict = {}
  for line in paramsfile:
      m = aprepro_regex.match(line)
      if m:
          paramsdict[m.group(1)] = m.group(2)
      else:
          m = standard_regex.match(line)
          if m:
              paramsdict[m.group(2)] = m.group(1)
  
  paramsfile.close()
  
  # crude error checking; handle both standard and aprepro cases
  num_vars = 0
  if ('variables' in paramsdict):
      num_vars = int(paramsdict['variables'])
  elif ('DAKOTA_VARS' in paramsdict):
      num_vars = int(paramsdict['DAKOTA_VARS'])
  
  num_fns = 0
  if ('functions' in paramsdict):
      num_fns = int(paramsdict['functions'])
  elif ('DAKOTA_FNS' in paramsdict):
      num_fns = int(paramsdict['DAKOTA_FNS'])
  
  # initialize dictionary
  fem_params =  {} 

  # -------------------------------
  # Convert and send to application
  # -------------------------------
  
  # set up the data structures the rosenbrock analysis code expects
  # for this simple example, put all the variables into a single hardwired array
  continuous_vars = {} 

  DescriptorList = ['robin_coeff','probe_init','mu_eff_healthy','body_temp','anfact_healthy', 'mu_a_healthy','mu_s_healthy','c_blood_healthy','c_p_healthy','rho_healthy','alpha_healthy','k_0_healthy','w_0_healthy','x_displace','y_displace','z_displace','x_rotate','y_rotate','z_rotate']
  for paramname in DescriptorList:
    try:
      continuous_vars[paramname  ] = paramsdict[paramname ]
    except KeyError:
      pass
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    try:
      active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
    except KeyError:
      active_set_vector = [ int(paramsdict['ASV_%d:obj_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  
  ################################
  # convert to uniform interface
  ################################
  #      mu_a_min               <      mu_a + (1-g) mu_s < mu_a_max + (1-g_min) mu_s_max
  #         5.e-1               <          mu_tr         < 600. + .3 * 50000. 
  #
  #  sqrt( 3 * 5.e-1 * 5.e-1 )  <  sqrt( 3 mu_a  mu_tr ) < sqrt( 3 * 600. * (600. + .3 * 50000.) ) 
  #  sqrt( 3 * 5.e-1 * 5.e-1 )  <        mu_eff          < sqrt( 3 * 600. * (600. + .3 * 50000.) ) 
  #            8.e-1            <        mu_eff          <    5.3e3
  import math
  mu_s   = float(continuous_vars['mu_s_healthy'])
  anfact = float(continuous_vars['anfact_healthy'])
  mu_s_p = mu_s * (1.-anfact) 
  # mu_tr  = mu_a + (1-g) mu_s 
  # mu_eff = sqrt( 3 mu_a  mu_tr )
  mu_eff = float(continuous_vars['mu_eff_healthy'])
  mu_a   =  0.5*( -mu_s_p + math.sqrt( mu_s_p * mu_s_p  + 4. * mu_eff * mu_eff  /3. ) )
  # alpha  == k / rho / c_p   W/m/K * m^3/kg * kg*K/W/s = m^2/s
  alpha   = float(continuous_vars['alpha_healthy'])
  rho     = float(continuous_vars['rho_healthy'])
  c_p     = float(continuous_vars['c_p_healthy'])
  c_blood = float(continuous_vars['c_blood_healthy'])
  k_0    = alpha * c_p * rho 
  w_0    = float(continuous_vars['w_0_healthy'])

  # store dakota vars
  continuous_vars['rho'    ]   =  rho   
  continuous_vars['c_p'    ]   =  c_p   
  continuous_vars['k_0'    ]   =  k_0   
  continuous_vars['w_0'    ]   =  w_0   
  continuous_vars['mu_a'   ]   =  mu_a  
  continuous_vars['mu_s'   ]   =  mu_s  
  continuous_vars['anfact' ]   =  anfact
  continuous_vars['c_blood']   =  c_blood 
  fem_params['cv']         = continuous_vars

  fem_params['asv']        = active_set_vector
  fem_params['functions']  = num_fns
  fem_params['fileID']     = fileID 
  fem_params['UID']        = int(paramfilename.split('/').pop(3))
  fem_params['opttype']    = paramfilename.split('.').pop(-3)
  fem_params['VisualizeOutput'] = VisualizeOutput 

  # parse file path
  locatemrti = paramfilename.split('/')
  locatemrti.pop()

  # database and run directory have the same structure
  fem_params['mrti']       = '%s/%s/%s/vtk/referenceBased/' % (databaseDIR,locatemrti[2],locatemrti[3])

  # get header info
  mrtifilename = '%s/temperature.%04d.vtk' % (fem_params['mrti'], 1) 
  print 'opening' , mrtifilename 
  vtkSetupImageReader = vtk.vtkDataSetReader() 
  vtkSetupImageReader.SetFileName(mrtifilename )
  vtkSetupImageReader.Update() 
  SetupImageData = vtkSetupImageReader.GetOutput() 
  fem_params['spacing']        = SetupImageData.GetSpacing()
  fem_params['dimensions']     = SetupImageData.GetDimensions()

  # get power file name
  inisetupfile  = "/".join(locatemrti)+"/setup.ini"
  config = ConfigParser.SafeConfigParser({})
  config.read(inisetupfile)
  if (not MatlabDriver):
       fem_params['lambdacode']       = eval(config.get('power','lambdacode'))
  fem_params['segment_file']     = config.get('exec','segment_file')
  fem_params['target_landmarks'] = config.get('exec','target_landmarks')
  #fem_params['powerhistory']     = config.get('power','history')
  fulltimeinterval               = eval(config.get('mrti','fulltime') )
  cooltimeinterval               = eval(config.get('mrti','cooling')  )
  heattimeinterval               = eval(config.get('mrti','heating')  )
  # select time interval selection from optimization type
  if(fem_params['opttype'] == 'heating'):
    timeinterval = heattimeinterval
  elif(fem_params['opttype'] == 'cooling'):
    timeinterval = cooltimeinterval
  else:
    timeinterval = fulltimeinterval
  fem_params['timeinterval'] = timeinterval
  fem_params['mrtideltat']   = config.getfloat('mrti','deltat') 
  fem_params['initialtime']  = timeinterval[0] * config.getfloat('mrti','deltat') 
  fem_params['finaltime']    = timeinterval[1] * config.getfloat('mrti','deltat') 
  fem_params['maxheatid']    = heattimeinterval[1]
  fem_params['voi']          = eval(config.get('mrti','voi'))

  print 'opttype',fem_params['opttype'],timeinterval ,'mrti data from' , fem_params['mrti'] , 'setupfile', inisetupfile  

  return fem_params
  ## ----------------------------
  ## Return the results to DAKOTA
  ## ----------------------------
  #
  #if (fem_results['rank'] == 0 ):
  #  # write the results.out file for return to DAKOTA
  #  # this example only has a single function, so make some assumptions;
  #  # not processing DVV
  #  outfile = open('results.out.tmp.%d' % fileID, 'w')
  #  
  #  # write functions
  #  for func_ind in range(0, num_fns):
  #      if (active_set_vector[func_ind] & 1):
  #          functions = fem_results['fns']    
  #          outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
  #  
  #  ## write gradients
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 2):
  #  #        grad = rosen_results['fnGrads'][func_ind]
  #  #        outfile.write('[ ')
  #  #        for deriv in grad: 
  #  #            outfile.write(str(deriv) + ' ')
  #  #        outfile.write(']\n')
  #  #
  #  ## write Hessians
  #  #for func_ind in range(0, num_fns):
  #  #    if (active_set_vector[func_ind] & 4):
  #  #        hessian = rosen_results['fnHessians'][func_ind]
  #  #        outfile.write('[[ ')
  #  #        for hessrow in hessian:
  #  #            for hesscol in hessrow:
  #  #                outfile.write(str(hesscol) + ' ')
  #  #            outfile.write('\n')
  #  #        outfile.write(']]')
  #  #
  #  outfile.close();outfile.flush
  #  #
  #  ## move the temporary results file to the one DAKOTA expects
  #  #import shutil
  #  #shutil.move('results.out.tmp.%d' % fileID, sys.argv[2])
# end def ParseInput:
##################################################################

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--run_fem","--param_file", 
                  action="store", dest="param_file", default=None,
                  help="run code with parameter FILE", metavar="FILE")
parser.add_option( "--vis_out", 
                  action="store_true", dest="vis_out", default=False,
                  help="visualise output", metavar="bool")
parser.add_option( "--accum_history", 
                  action="store", dest="accum_history", default=None,
                  help="accumulate output from opttype", metavar="bool")
parser.add_option( "--run_min", 
                  action="store", dest="run_min", default=None,
                  help="re-run the optimum", metavar="FILE")
parser.add_option( "--ini", 
                  action="store", dest="config_ini", default=None,
                  help="ini FILE containing setup info", metavar="FILE")
(options, args) = parser.parse_args()

if (options.param_file != None):
  # parse the dakota input file
  fem_params = ParseInput(options.param_file,options.vis_out)

  if(MatlabDriver):
    print fem_params
    import scipy.io as scipyio
    # write out for debug
    fem_params['patientID'] = options.param_file.split('/')[2]
    fem_params['UID']       = options.param_file.split('/')[3]
    #scipyio.savemat( '%s.mat' % options.param_file, MatlabDataDictionary )
    scipyio.savemat( './TmpDataInput.mat' , fem_params )
    # FIXME setup any needed paths
    # FIXME this nees to have a clean matlab env for dakmatlab
    # FIXME then setup ONCE
    #os.system( './analytic/dakmatlab setup workspace ' )
    matlabcommand  = './analytic/dakmatlab %s %s' %  (options.param_file,sys.argv[3])
    print matlabcommand  
    os.system( matlabcommand )
  else:
    # FIXME link needed directories
    linkDirectoryList = ['occa','libocca','meshes']
    for targetDirectory in linkDirectoryList:
      linkcommand = 'ln -sf %s/%s .' % (brainNekDIR,targetDirectory )
      print linkcommand 
      os.system(linkcommand )

    # execute the rosenbrock analysis as a separate Python module
    print "Running BrainNek..."
    brainNekWrapper(**fem_params)
    
    # write objective function back to Dakota
    objfunctionlist = ComputeObjective(**fem_params)

    print "current objective function: ",objfunctionlist 
    fileHandle = file(sys.argv[3],'w')
    for objfncvalue in objfunctionlist:
      fileHandle.write('%f\n' % objfncvalue )
    fileHandle.flush(); fileHandle.close();

# find the best point for each run
elif (options.accum_history ):
  resultfileList = [
  './workdir/Study0035/0530/',
  #'./workdir/Study0023/0433/',
  #'./workdir/Study0023/0428/',
  ##'./workdir/Study0023/0425/',
  './workdir/Study0030/0495/',
  './workdir/Study0030/0497/',
  './workdir/Study0030/0488/',
  './workdir/Study0030/0491/',
  './workdir/Study0030/0496/',
  './workdir/Study0030/0490/',
  './workdir/Study0017/0378/',
  ##'./workdir/Study0018/0388/',
  #'./workdir/Study0018/0402/',
  #'./workdir/Study0018/0389/',
  #'./workdir/Study0018/0385/',
  #'./workdir/Study0029/0476/',
  #'./workdir/Study0029/0477/',
  './workdir/Study0025/0438/',
  './workdir/Study0025/0435/',
  './workdir/Study0025/0440/',
  './workdir/Study0025/0436/',
  './workdir/Study0028/0466/',
  './workdir/Study0028/0468/',
  './workdir/Study0028/0471/',
  #'./workdir/Study0052/0725/',
  #'./workdir/Study0052/0720/',
  './workdir/Study0026/0447/',
  './workdir/Study0026/0457/',
  './workdir/Study0026/0455/',
  './workdir/Study0026/0453/',
  './workdir/Study0026/0450/',
  './workdir/Study0026/0451/',
  ##'./workdir/Study0057/0772/',
  ##'./workdir/Study0057/0769/',
  './workdir/Study0022/0418/',
  './workdir/Study0022/0417/',
  './workdir/Study0021/0409/',
  './workdir/Study0021/0414/',
  './workdir/Study0021/0415/',
  ##'./workdir/Study0054/0753/',
  ##'./workdir/Study0054/0756/',
  ##'./workdir/Study0053/0755/',
  ##'./workdir/Study0006/0183/',
  ]
  
  ## resultfileList = [
  ## './workdir/Study0035/0530/',
  ## './workdir/Study0030/0491/',
  ## ]
  
  texHandle  = open('datasummary.tex' , 'w') 
  fileHandle = open('datasummary.txt' , 'w') 
  # write header
  fileHandle.write("idstudy,iddata,idmin,mu_eff,alpha,robin,dice,obj\n")
  # loop over files and extract optimal value
  opttype = options.accum_history 
  for filenamebase in resultfileList:
    # get latex command
    config = ConfigParser.SafeConfigParser({})
    inisetupfile = '%s/opt/setup.ini' % (filenamebase)
    config.read(inisetupfile)
  
    # get min value
    (idmin,minobjval,dicevalue ) = GetMinJobID( '%s/opt/optpp_pds.%s' % (filenamebase,opttype))
    print (idmin,minobjval,dicevalue ) 
    
    studyid= int(filenamebase.split('/')[2].replace('Study',''))
    dataid = int(filenamebase.split('/')[3])
    # count the file lines
    dakotafilename = '%s/opt/optpp_pds.%s.in.%d' % (filenamebase,opttype,idmin)
    opt_fem_params = ParseInput(dakotafilename,False)
    simvariable = opt_fem_params['cv']     
    # get arrhenius dice value
    heattimeinterval               = eval(config.get('mrti','heating')  )
    SEMDataDirectory               = outputDirectory % int(filenamebase.split('/')[-2]) 
    #dataarray = numpy.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4,6)
    fileHandle.write("%05d,%05d,%05d,%s,%s,%s,%12.5e,%12.5e\n" %( studyid, dataid, idmin     ,
                                                                    simvariable['mu_eff_healthy'],
                                                                    simvariable['alpha_healthy'],
                                                                    simvariable['robin_coeff'],
                                                                    dicevalue,
                                                                    minobjval))
    # format latex ouput
    outputformat                   = config.get('latex','opttype')
    texFormat = outputformat % (opttype,heattimeinterval[1],opttype,heattimeinterval[1],opttype,heattimeinterval[1],opttype,heattimeinterval[1],opttype,heattimeinterval[1],minobjval,dicevalue)
    #print texFormat 
    texHandle.write("%s\n" %(texFormat))

  texHandle.close() 
  fileHandle.close()

# rerun the optimizer at the minimum
elif (options.run_min != None):

  templatefilename = options.run_min
  # get min value
  (idmin,minobjval,dicevalue) = GetMinJobID( templatefilename )
  print (idmin,minobjval,dicevalue) 

  # build execution command
  runcmd = "vglrun python ./brainsearch.py --param_file  %s.in.%d %s.out.%d --vis_out" % (templatefilename,idmin,templatefilename,idmin)
  print runcmd
  #FIXME not running ???
  os.system( runcmd )
  

# run planning solver w/ default options from ini file
elif (options.config_ini != None):

  # read config file
  config = ConfigParser.SafeConfigParser({})
  config.read(options.config_ini)
  
  fem_params = {}
  fem_params['UID']           =  0000
  fem_params['fileID']        = 0
  fem_params['segment_file']  = config.get('exec','segment_file')
  # store the entire configuration file for convienence
  fem_params['config_parser'] = config
  
  # write initial slicer.ini
  SlicerIniFilename = "./slicer.ini"
  initialconfig = ConfigParser.SafeConfigParser({})
  initialconfig.add_section("timestep")
  initialconfig.add_section("exec")
  initialconfig.set("timestep","power",config.get('timestep','power'))
  initialconfig.set("timestep","finaltime","10.0")
  initialconfig.set("exec","target_landmarks"  ,"./TargetLandmarksvtk.vtk")
  initialconfig.set("exec","transformlandmarks","./TargetLandmarksras.vtk")

  with open(SlicerIniFilename , 'w') as configfile:
    initialconfig.write(configfile)
  
  # time stamp
  import time
  timeStamp =0 
  while(True):
    if(os.path.getmtime( SlicerIniFilename ) > timeStamp):
        timeStamp = os.path.getmtime(SlicerIniFilename ) 
        slicerconfig = ConfigParser.SafeConfigParser({})
        slicerconfig.read( SlicerIniFilename )
        fem_params['deltat']        =  5.0
        fem_params['finaltime']     =  slicerconfig.getfloat('timestep','finaltime')
        # build lambda funtion for power history
        fem_params['lambdacode']    =  lambda time:  0.0 if time < fem_params['deltat'] else slicerconfig.getfloat('timestep','power') if time < fem_params['finaltime'] else 0.0
        fem_params['target_landmarks']   = slicerconfig.get('exec','target_landmarks'  )
        fem_params['transformlandmarks'] = slicerconfig.get('exec','transformlandmarks')
        # set tissue lookup tables
        k_0Table  = {"default":config.getfloat("thermal_conductivity","k_0_healthy")  ,
                     "vessel" :config.getfloat("thermal_conductivity","k_0_healthy")  ,
                     "grey"   :config.getfloat("thermal_conductivity","k_0_grey"   )  ,
                     "white"  :config.getfloat("thermal_conductivity","k_0_white"  )  ,
                     "csf"    :config.getfloat("thermal_conductivity","k_0_csf"    )  ,
                     "tumor"  :config.getfloat("thermal_conductivity","k_0_tumor"  )  }
        w_0Table  = {"default":config.getfloat("perfusion","w_0_healthy")  ,
                     "vessel" :config.getfloat("perfusion","w_0_healthy")  ,
                     "grey"   :config.getfloat("perfusion","w_0_grey"   )  ,
                     "white"  :config.getfloat("perfusion","w_0_white"  )  ,
                     "csf"    :config.getfloat("perfusion","w_0_csf"    )  ,
                     "tumor"  :config.getfloat("perfusion","w_0_tumor"  )  }
        mu_aTable = {"default":config.getfloat("optical","mu_a_healthy")  ,
                     "vessel" :config.getfloat("optical","mu_a_healthy")  ,
                     "grey"   :config.getfloat("optical","mu_a_grey"   )  ,
                     "white"  :config.getfloat("optical","mu_a_white"  )  ,
                     "csf"    :config.getfloat("optical","mu_a_csf"    )  ,
                     "tumor"  :config.getfloat("optical","mu_a_tumor"  )  }
        mu_sTable = {"default":config.getfloat("optical","mu_s_healthy")  ,
                     "vessel" :config.getfloat("optical","mu_s_healthy")  ,
                     "grey"   :config.getfloat("optical","mu_s_grey"   )  ,
                     "white"  :config.getfloat("optical","mu_s_white"  )  ,
                     "csf"    :config.getfloat("optical","mu_s_csf"    )  ,
                     "tumor"  :config.getfloat("optical","mu_s_tumor"  )  }
        anfactTable={"default":config.getfloat("optical","anfact_healthy")  ,
                     "vessel" :config.getfloat("optical","anfact_healthy")  ,
                     "grey"   :config.getfloat("optical","anfact_grey"   )  ,
                     "white"  :config.getfloat("optical","anfact_white"  )  ,
                     "csf"    :config.getfloat("optical","anfact_csf"    )  ,
                     "tumor"  :config.getfloat("optical","anfact_tumor"  )  }
        labelTable= {config.get("labels","greymatter" ):"grey" , 
                     config.get("labels","whitematter"):"white", 
                     config.get("labels","csf"        ):"csf"  , 
                     config.get("labels","tumor"      ):"tumor", 
                     config.get("labels","vessel"     ):"vessel"}
        labelCount= {"default":0,
                     "grey"   :0, 
                     "white"  :0, 
                     "csf"    :0, 
                     "tumor"  :0, 
                     "vessel" :0}
        # store constitutive data
        continuous_vars  = {}
        continuous_vars['rho'    ]   =  1045.
        continuous_vars['c_p'    ]   =  3640.
        continuous_vars['c_blood']   =  3840.
        continuous_vars['k_0'    ]   =  k_0Table["default"]
        continuous_vars['w_0'    ]   =  w_0Table["default"]
        continuous_vars['mu_a'   ]   =  mu_aTable["default"] 
        continuous_vars['mu_s'   ]   =  mu_sTable["default"]
        continuous_vars['anfact' ]   =  anfactTable["default"]
        continuous_vars['body_temp'] = config.getfloat("initial_condition","u_init"  ) 
        continuous_vars['probe_init'] = 21.0
        continuous_vars['x_displace'] = 0.0
        continuous_vars['y_displace'] = 0.0
        continuous_vars['z_displace'] = 0.0
        continuous_vars['x_rotate']   = 0.0
        continuous_vars['y_rotate']   = 0.0
        continuous_vars['z_rotate']   = 0.0
        fem_params['cv']         = continuous_vars
        # execute 
        print "Running BrainNek..."
        brainNekWrapper(**fem_params)
        
        # write objective function 
        objfunction = ForwardSolve(**fem_params)
    else:
      print "waiting on user input..",SlicerIniFilename 
      # echo lookup table
      print "lookup tables"
      #print "labeled %d voxels" % len(imageLabel)
      print "labels"      , labelTable
      print "counts"      , labelCount
      print "conductivity", k_0Table  
      print "perfusion"   , w_0Table  
      print "absorption"  , mu_aTable  
      print "scattering"  , mu_sTable  
      print "anfact"      , anfactTable  
      time.sleep(2)
else:
  parser.print_help()
  print options
