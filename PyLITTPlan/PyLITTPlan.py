import os
import unittest
import numpy
import ConfigParser
import time
from __main__ import vtk, qt, ctk, slicer

#
# PyLITTPlan
#

class PyLITTPlan:
  def __init__(self, parent):
    parent.title = "PyLITTPlan" # TODO make this more human readable by adding spaces
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["David Fuentes (MD Anderson)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc. and Steve Pieper, Isomics, Inc.  and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
    self.parent = parent

    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.
    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['PyLITTPlan'] = self.runTest

  def runTest(self):
    tester = PyLITTPlanTest()
    tester.runTest()

#
# qPyLITTPlanWidget
#

class PyLITTPlanWidget:
  def __init__(self, parent = None):
    self.developerMode = True # change this to true to get reload and test
    self.developerMode = False # change this to true to get reload and test
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.inputFiducialsNodeSelector.setMRMLScene(slicer.mrmlScene)
      self.parent.show()

  def setup(self):
    # Instantiate and connect widgets ...
    self.ApplicatorModel           = None
    self.ApplicatorTrajectoryModel = None
    self.DamageTemplateModel       = None
    self.SolverDamageModel         = None
    self.SourceLandmarkFileName    = "./SourceLandmarks.vtk"
    self.TargetLandmarkFileName    = "./TargetLandmarks.vtk"

    # define constant parameters
    #  diffusing tip is 10. mm axial
    #  diffusing tip is 1.5 mm radial
    self.DiffusingTipLength = 10. #mm
    self.DiffusingTipRadius = .75 #mm
    # In canine brain and transmissible venereal tumours, up to 18.11.4mm lesions were achieved. It is concluded
    self.AblationMinorAxis = 18.0/2. #mm
    self.AblationMajorAxis = 22.0/2. #mm
    #
    # Reload and Test area
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Main"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    if self.developerMode:
      # reload button
      # (use this during development, but remove it when delivering
      #  your module to users)
      self.reloadButton = qt.QPushButton("Reload")
      self.reloadButton.toolTip = "Reload this module."
      self.reloadButton.name = "PyLITTPlan Reload"
      reloadFormLayout.addWidget(self.reloadButton)
      self.reloadButton.connect('clicked()', self.onReload)

      # reload and test button
      # (use this during development, but remove it when delivering
      #  your module to users)
      self.reloadAndTestButton = qt.QPushButton("Reload and Test")
      self.reloadAndTestButton.toolTip = "Reload this module and then run the self tests."
      reloadFormLayout.addWidget(self.reloadAndTestButton)
      self.reloadAndTestButton.connect('clicked()', self.onReloadAndTest)

    # Input fiducials node selector
    inputFiducialsNodeSelector = slicer.qMRMLNodeComboBox()
    inputFiducialsNodeSelector.nodeTypes = ['vtkMRMLMarkupsFiducialNode', 'vtkMRMLAnnotationHierarchyNode', 'vtkMRMLFiducialListNode']
    inputFiducialsNodeSelector.objectName = 'inputFiducialsNodeSelector'
    inputFiducialsNodeSelector.selectNodeUponCreation = True
    inputFiducialsNodeSelector.noneEnabled = False
    inputFiducialsNodeSelector.addEnabled  = False
    inputFiducialsNodeSelector.removeEnabled = False
    inputFiducialsNodeSelector.showHidden = False
    inputFiducialsNodeSelector.showChildNodeTypes = False
    inputFiducialsNodeSelector.setMRMLScene( slicer.mrmlScene )
    self.inputFiducialsNodeSelector = inputFiducialsNodeSelector
    inputFiducialsNodeSelector.toolTip = "Select a fiducial list to define control points for the path."
    reloadFormLayout.addRow("Input Fiducials:", self.inputFiducialsNodeSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 1 )
    self.outputSelector.selectNodeUponCreation = False
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = False
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    reloadFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # Power Value
    #
    self.PowerValueSliderWidget = ctk.ctkSliderWidget()
    self.PowerValueSliderWidget.singleStep = 0.5
    self.PowerValueSliderWidget.minimum = 1.0
    self.PowerValueSliderWidget.maximum = 30.0
    self.PowerValueSliderWidget.value = 12.0
    self.PowerValueSliderWidget.setToolTip("Set Power Value")
    reloadFormLayout.addRow("Power [W]", self.PowerValueSliderWidget)

    #
    # Time Length Value
    #
    self.TimeValueSliderWidget = ctk.ctkSliderWidget()
    self.TimeValueSliderWidget.singleStep = 5.0
    self.TimeValueSliderWidget.minimum = 5.0
    self.TimeValueSliderWidget.maximum = 300.0
    self.TimeValueSliderWidget.value = 50.0
    self.TimeValueSliderWidget.setToolTip("Laser On Duration")
    reloadFormLayout.addRow("Time [s]", self.TimeValueSliderWidget)

    #
    # Damage Iso Value
    #
    self.DamageValueSliderWidget = ctk.ctkSliderWidget()
    self.DamageValueSliderWidget.singleStep = 1.0
    self.DamageValueSliderWidget.minimum = 0.01
    self.DamageValueSliderWidget.maximum = 100.0
    self.DamageValueSliderWidget.value = 57.0
    self.DamageValueSliderWidget.setToolTip("Set Damage Value")
    reloadFormLayout.addRow("Damage", self.DamageValueSliderWidget)

    #
    # Locate Applicator Button
    #
    self.referenceButton = qt.QPushButton("Locate Applicator")
    self.referenceButton.toolTip = "Locate the Applicator for Reference"
    self.referenceButton.enabled = False
    reloadFormLayout.addRow(self.referenceButton)

    #
    # Treat Button
    #
    self.applyButton = qt.QPushButton("Treat")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    reloadFormLayout.addRow(self.applyButton)

    #
    # Damage Iso Value
    #
    self.PullBackValueSliderWidget = ctk.ctkSliderWidget()
    self.PullBackValueSliderWidget.singleStep = 0.5
    self.PullBackValueSliderWidget.minimum = -80.0
    self.PullBackValueSliderWidget.maximum =  80.0
    self.PullBackValueSliderWidget.value = 0.0
    self.PullBackValueSliderWidget.setToolTip("Set Pullback Length")
    reloadFormLayout.addRow("Pullback [mm]", self.PullBackValueSliderWidget)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    parametersCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    #
    # scale factor for screen shots
    #
    self.PerfusionValueSliderWidget = ctk.ctkSliderWidget()
    self.PerfusionValueSliderWidget.singleStep = 1.0
    self.PerfusionValueSliderWidget.minimum = 0.0
    self.PerfusionValueSliderWidget.maximum = 50.0
    self.PerfusionValueSliderWidget.value = 6.0
    self.PerfusionValueSliderWidget.setToolTip("Set Perfusion Value")
    parametersFormLayout.addRow("Perfusion [kg/s/m^3]", self.PerfusionValueSliderWidget)

    #
    # scale factor for screen shots
    #
    self.AbsorptionValueSliderWidget = ctk.ctkSliderWidget()
    self.AbsorptionValueSliderWidget.singleStep = 0.5
    self.AbsorptionValueSliderWidget.minimum = 0.05
    self.AbsorptionValueSliderWidget.maximum = 5.e2
    self.AbsorptionValueSliderWidget.value = 5.0
    self.AbsorptionValueSliderWidget.setToolTip("Set mu_a Value")
    parametersFormLayout.addRow("mu_a [1/m]", self.AbsorptionValueSliderWidget)

    # connections
    self.inputFiducialsNodeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.referenceButton.connect('clicked(bool)', self.onReferenceButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Add vertical spacer
    self.layout.addStretch(1)

  #def enableOrDisableCreateButton(self):
  #  """Connected to both the fiducial and camera node selector. It allows to
  #  enable or disable the 'create path' button."""
  #  self.createPathButton.enabled = self.cameraNodeSelector.currentNode() != None and self.inputFiducialsNodeSelector.currentNode() != None

  def cleanup(self):
    pass

  def onSelect(self):
    #self.applyButton.enabled = self.inputFiducialsNodeSelector.currentNode() and self.outputSelector.currentNode()
    self.applyButton.enabled = self.inputFiducialsNodeSelector.currentNode() 
    self.referenceButton.enabled = self.applyButton.enabled 

  # write vtk points file
  def WriteVTKPoints(self,vtkpoints,OutputFileName):
     # get data structures to convert to pixel ijk space
     outnode = self.outputSelector.currentNode()
     ras2ijk = vtk.vtkMatrix4x4()
     outnode.GetRASToIJKMatrix(ras2ijk)

     # loop over points an store in vtk data structure
     # write in pixel coordinates
     scalevtkPoints = vtk.vtkPoints()
     vertices= vtk.vtkCellArray()
     for idpoint in range(vtkpoints.GetNumberOfPoints()):
         currentpoint = vtkpoints.GetPoint(idpoint)
         ijkpoint = ras2ijk.MultiplyPoint( (currentpoint[0],currentpoint[1],currentpoint[2],1.) )
         print ijkpoint 
         vertices.InsertNextCell( 1 ); vertices.InsertCellPoint( scalevtkPoints.InsertNextPoint(ijkpoint[0:3]) )
         #vertices.InsertNextCell( 1 ); vertices.InsertCellPoint( idpoint )
  
     # set polydata
     polydata = vtk.vtkPolyData()
     polydata.SetPoints(scalevtkPoints )
     polydata.SetVerts( vertices )
  
     # write to file
     print "WriteVTKPoints: writing",OutputFileName
     polydatawriter = vtk.vtkDataSetWriter()
     polydatawriter.SetFileName(OutputFileName)
     polydatawriter.SetInput(polydata)
     polydatawriter.Update()

  def GetDiffusingLaserTransform(self):
    fiducialsNode = self.inputFiducialsNodeSelector.currentNode();
    print fiducialsNode.GetClassName() 
    if fiducialsNode.GetClassName() == "vtkMRMLMarkupsFiducialNode":
      # slicer4 Markups node
      NumFid = fiducialsNode.GetNumberOfFiducials()
      if NumFid < 2:
        print("Two Fiducials Needed ")
        return None
      # get fiducial positions
      else:
        # slicer Points
        pointtip   = [0.0, 0.0, 0.0]
        fiducialsNode.GetNthFiducialPosition(0, pointtip   )
        pointentry = [0.0, 0.0, 0.0]
        fiducialsNode.GetNthFiducialPosition(1, pointentry)
    else:
      print("Unknown Class ")
      return None

    # diffusing applicator center at coordinate  (0,                        0., 0. ) mm
    # template laser distal ends  at coordinates (0, +/- DiffusingTipLength/2., 0. ) mm
    originalOrientation = vtk.vtkPoints()
    originalOrientation.SetNumberOfPoints(2)
    originalOrientation.SetPoint(0,0.,                        0.,0.)
    originalOrientation.SetPoint(1,0.,self.DiffusingTipLength/2.,0.)
    slicerLength   = numpy.linalg.norm( numpy.array(pointentry) - numpy.array(pointtip) )
    unitdirection  = 1./slicerLength * (numpy.array(pointentry) - numpy.array(pointtip) ) 
    pointtip    = pointtip - self.PullBackValueSliderWidget.value * unitdirection    
    pointscaled = pointtip - self.PullBackValueSliderWidget.value * unitdirection + self.DiffusingTipLength/2. * unitdirection
    print "points", pointentry, pointtip, pointscaled, slicerLength, numpy.linalg.norm( unitdirection  ), numpy.linalg.norm( pointscaled - pointtip ) 
    slicerOrientation   = vtk.vtkPoints()
    slicerOrientation.SetNumberOfPoints(2)
    slicerOrientation.SetPoint(0,pointtip[   0],pointtip[   1],pointtip[   2] )
    slicerOrientation.SetPoint(1,pointscaled[0],pointscaled[1],pointscaled[2] )

    LaserLineTransform = vtk.vtkLandmarkTransform()
    LaserLineTransform.SetSourceLandmarks(originalOrientation)
    LaserLineTransform.SetTargetLandmarks(slicerOrientation  )
    LaserLineTransform.SetModeToRigidBody()
    LaserLineTransform.Update()
    print LaserLineTransform.GetMatrix()

    ## # Get RAS transformation
    ## rasToXY = vtk.vtkMatrix4x4()


    ## layoutManager = slicer.app.layoutManager() 
    ## sliceWidget = layoutManager.sliceWidget('Red')
    ## sliceLogic = sliceWidget.sliceLogic()
    ## layerLogic = sliceLogic.GetBackgroundLayer()
    ## xyToIJK = layerLogic.GetXYToIJKTransform().GetMatrix() 
    ## rasToXY.DeepCopy( xyToIJK )
    ## #rasToXY.Invert()
    ## outnode = self.outputSelector.currentNode()
    ## ras2ijk = vtk.vtkMatrix4x4()
    ## ijk2ras = vtk.vtkMatrix4x4()
    ## outnode.GetRASToIJKMatrix(ras2ijk)
    ## outnode.GetIJKToRASMatrix(ijk2ras)
    ## print ras2ijk 
    ## print ijk2ras 
    ## print rasToXY
    ## print pointtip
    ## print rasToXY.MultiplyPoint( (pointtip[0],pointtip[1],pointtip[2],1.) )
    ## print ras2ijk.MultiplyPoint( (pointtip[0],pointtip[1],pointtip[2],1.) )

    return LaserLineTransform

  def AddApplicatorModel(self,scene,LineLandmarkTransform ):
      # check if already implemented
      if ( self.ApplicatorModel != None ):
        print "removing", self.ApplicatorModel.GetName()
        scene.RemoveNode(self.ApplicatorModel)

      # Define applicator tip
      vtkCylinder = vtk.vtkCylinderSource()
      vtkCylinder.SetHeight(self.DiffusingTipLength ); 
      vtkCylinder.SetRadius(self.DiffusingTipRadius );
      vtkCylinder.SetCenter(0.0, 0.0, 0.0);
      vtkCylinder.SetResolution(16);

      # transform
      slicertransformFilter = vtk.vtkTransformFilter()
      slicertransformFilter.SetInput(vtkCylinder.GetOutput() ) 
      slicertransformFilter.SetTransform( LineLandmarkTransform ) 
      slicertransformFilter.Update()
      apppolyData=slicertransformFilter.GetOutput();

      # add to scene
      self.ApplicatorModel = slicer.vtkMRMLModelNode()
      self.ApplicatorModel.SetScene(scene)
      #self.ApplicatorModel.SetName(scene.GenerateUniqueName("Applicator-%s" % fiducialsNode.GetName()))
      self.ApplicatorModel.SetName("Applicator" )
      self.ApplicatorModel.SetAndObservePolyData(apppolyData)

      # Create display node
      modelDisplay = slicer.vtkMRMLModelDisplayNode()
      modelDisplay.SetColor(1,0,0) # red
      modelDisplay.SetScene(scene)
      scene.AddNode(modelDisplay)
      self.ApplicatorModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

      # Add to scene
      modelDisplay.SetInputPolyData(self.ApplicatorModel.GetPolyData())
      scene.AddNode(self.ApplicatorModel)

  def AddApplicatorTrajectoryModel(self,scene,LineLandmarkTransform ):
      # check if already implemented
      if ( self.ApplicatorTrajectoryModel != None ):
        print "removing", self.ApplicatorTrajectoryModel.GetName()
        scene.RemoveNode(self.ApplicatorTrajectoryModel)

      # Define applicator tip
      vtkCylinder = vtk.vtkCylinderSource()
      vtkCylinder.SetHeight(10.*self.DiffusingTipLength); 
      vtkCylinder.SetRadius(0.1*self.DiffusingTipRadius);
      vtkCylinder.SetCenter(0.0,  10.*self.DiffusingTipLength/2.0, 0.0);
      vtkCylinder.SetResolution(16);

      # transform
      slicertransformFilter = vtk.vtkTransformFilter()
      slicertransformFilter.SetInput(vtkCylinder.GetOutput() ) 
      slicertransformFilter.SetTransform( LineLandmarkTransform ) 
      slicertransformFilter.Update()
      apppolyData=slicertransformFilter.GetOutput();

      # add to scene
      self.ApplicatorTrajectoryModel = slicer.vtkMRMLModelNode()
      self.ApplicatorTrajectoryModel.SetScene(scene)
      #self.ApplicatorTrajectoryModel.SetName(scene.GenerateUniqueName("Applicator-%s" % fiducialsNode.GetName()))
      self.ApplicatorTrajectoryModel.SetName("ApplicatorTrajectory" )
      self.ApplicatorTrajectoryModel.SetAndObservePolyData(apppolyData)

      # Create display node
      modelDisplay = slicer.vtkMRMLModelDisplayNode()
      modelDisplay.SetColor(1,0.5,0) # red
      modelDisplay.SetScene(scene)
      scene.AddNode(modelDisplay)
      self.ApplicatorTrajectoryModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

      # Add to scene
      modelDisplay.SetInputPolyData(self.ApplicatorTrajectoryModel.GetPolyData())
      scene.AddNode(self.ApplicatorTrajectoryModel)


  def AddDamageTemplateModel(self,scene,LineLandmarkTransform ):
      # check if already implemented
      if ( self.DamageTemplateModel != None ):
        print "removing", self.DamageTemplateModel.GetName()
        scene.RemoveNode(self.DamageTemplateModel)
      # create sphere and scale
      vtkSphere = vtk.vtkSphereSource()
      vtkSphere.SetRadius(1.) 
      vtkSphere.SetThetaResolution(16)
      vtkSphere.SetPhiResolution(16)
      ScaleAffineTransform = vtk.vtkTransform()
      ScaleAffineTransform.Scale([self.AblationMinorAxis,self.AblationMajorAxis,self.AblationMinorAxis])
      vtkEllipsoid= vtk.vtkTransformFilter()
      vtkEllipsoid.SetInput(vtkSphere.GetOutput() ) 
      vtkEllipsoid.SetTransform( ScaleAffineTransform ) 
      vtkEllipsoid.Update()

      slicertransformFilter = vtk.vtkTransformFilter()
      slicertransformFilter.SetInput(vtkEllipsoid.GetOutput() ) 
      slicertransformFilter.SetTransform( LineLandmarkTransform ) 
      slicertransformFilter.Update()
      dampolyData=slicertransformFilter.GetOutput();

      # Create model node
      self.DamageTemplateModel = slicer.vtkMRMLModelNode()
      self.DamageTemplateModel.SetScene(scene)
      #self.DamageTemplateModel.SetName(scene.GenerateUniqueName("Treatment-%s" % fiducialsNode.GetName()))
      self.DamageTemplateModel.SetName( scene.GenerateUniqueName("Treatment") )
      self.DamageTemplateModel.SetAndObservePolyData(dampolyData)

      # Create display node
      modelDisplay = slicer.vtkMRMLModelDisplayNode()
      modelDisplay.SetColor(1,1,0) # yellow
      modelDisplay.SetOpacity(0.3) 
      modelDisplay.SetScene(scene)
      scene.AddNode(modelDisplay)
      self.DamageTemplateModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

      # Add to scene
      modelDisplay.SetInputPolyData(self.DamageTemplateModel.GetPolyData())
      scene.AddNode(self.DamageTemplateModel)

  def onReferenceButton(self):
    logic = PyLITTPlanLogic()
    screenshotScaleFactor = int(self.PerfusionValueSliderWidget.value)
    print("Run the algorithm")

    # get fiducial based transform
    LineLandmarkTransform = self.GetDiffusingLaserTransform()
    if ( LineLandmarkTransform == None): 
      return
    print  "Model Params:", self.PowerValueSliderWidget.value, self.AbsorptionValueSliderWidget.value 

    # TODO logic.run for brainNek
    #logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), enableScreenshotsFlag,screenshotScaleFactor)

    # get slicer scene
    scene = slicer.mrmlScene

    # display applicator tip
    DisplayApplicatorTipOnLocate = False
    DisplayApplicatorTipOnLocate = True
    if (DisplayApplicatorTipOnLocate):
      self.AddApplicatorModel(scene,LineLandmarkTransform )
      self.AddApplicatorTrajectoryModel(scene,LineLandmarkTransform )

    # display applicator tip
    DisplayDamageTemplateOnLocate = False
    DisplayDamageTemplateOnLocate = True
    if (DisplayDamageTemplateOnLocate ):
      self.AddDamageTemplateModel(scene,LineLandmarkTransform )

  def onApplyButton(self):
    logic = PyLITTPlanLogic()
    screenshotScaleFactor = int(self.PerfusionValueSliderWidget.value)
    print("Run the algorithm")

    # get fiducial based transform
    LineLandmarkTransform = self.GetDiffusingLaserTransform()
    if ( LineLandmarkTransform == None): 
      return
    print  "Model Params:", self.PowerValueSliderWidget.value, self.AbsorptionValueSliderWidget.value 

    # get slicer scene
    scene = slicer.mrmlScene

    # display ellipse if solver not available
    GPUSolverAvailable = False
    GPUSolverAvailable = True
    if ( GPUSolverAvailable ) :
      # TODO logic.run for brainNek
      #logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), enableScreenshotsFlag,screenshotScaleFactor)
      self.WriteVTKPoints(LineLandmarkTransform.GetSourceLandmarks(),self.SourceLandmarkFileName )
      self.WriteVTKPoints(LineLandmarkTransform.GetTargetLandmarks(),self.TargetLandmarkFileName )
      SlicerIniFilename = "./slicer.ini"
      initialconfig = ConfigParser.SafeConfigParser({})
      initialconfig.add_section("timestep")
      initialconfig.add_section("exec")
      initialconfig.set("timestep","power","%s" % self.PowerValueSliderWidget.value)
      initialconfig.set("timestep","finaltime","%s" % self.TimeValueSliderWidget.value)
      initialconfig.set("exec","target_landmarks","./TargetLandmarks.vtk")
      with open(SlicerIniFilename , 'w') as configfile:
        initialconfig.write(configfile)
      timeStamp = time.time()
      WaitForSolver = True
      while(WaitForSolver ):
        if( os.path.getmtime( "./fem.finish" ) > timeStamp):
          if ( self.SolverDamageModel !=  None):
             print "removing", self.SolverDamageModel.GetName()
             scene.RemoveNode(self.SolverDamageModel)
          self.SolverDamageModel = slicer.util.loadModel("./fem.stl",True)[1]
          WaitForSolver  = False
        else:
          print "waiting on solver time ..",time.time(),os.path.getmtime( "./fem.finish" ) , timeStamp
          time.sleep(1)

  def onReload(self,moduleName="PyLITTPlan"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)

  def onReloadAndTest(self,moduleName="PyLITTPlan"):
    try:
      self.onReload()
      evalString = 'globals()["%s"].%sTest()' % (moduleName, moduleName)
      tester = eval(evalString)
      tester.runTest()
    except Exception, e:
      import traceback
      traceback.print_exc()
      qt.QMessageBox.warning(slicer.util.mainWindow(),
          "Reload and Test", 'Exception!\n\n' + str(e) + "\n\nSee Python Console for Stack Trace")


#
# PyLITTPlanLogic
#   TODO add brainNek Here
#

class PyLITTPlanLogic:
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """
  def __init__(self):
    pass

  def hasImageData(self,volumeNode):
    """This is a dummy logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      print('no volume node')
      return False
    if volumeNode.GetImageData() == None:
      print('no image data')
      return False
    return True

  def delayDisplay(self,message,msec=1000):
    #
    # logic version of delay display
    #
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    self.delayDisplay(description)

    if self.enableScreenshots == 0:
      return

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == -1:
      # full window
      widget = slicer.util.mainWindow()
    elif type == slicer.qMRMLScreenShotDialog().FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog().ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog().Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog().Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog().Green:
      # green slice window
      widget = lm.sliceWidget("Green")

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, self.screenshotScaleFactor, imageData)

  def run(self,inputVolume,outputVolume,enableScreenshots=0,screenshotScaleFactor=1):
    """
    Run the actual algorithm
    """

    self.delayDisplay('Running the aglorithm')

    self.enableScreenshots = enableScreenshots
    self.screenshotScaleFactor = screenshotScaleFactor

    self.takeScreenshot('PyLITTPlan-Start','Start',-1)

    return True


class PyLITTPlanTest(unittest.TestCase):
  """
  This is the test case for your scripted module.
  """

  def delayDisplay(self,message,msec=1000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PyLITTPlan1()

  def test_PyLITTPlan1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests sould exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        print('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        print('Loading %s...\n' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading\n')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = PyLITTPlanLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
