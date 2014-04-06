import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# PyLITTPlan
#

class PyLITTPlan:
  def __init__(self, parent):
    parent.title = "PyLITTPlan" # TODO make this more human readable by adding spaces
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["Jean-Christophe Fillion-Robin (Kitware), Steve Pieper (Isomics)"] # replace with "Firstname Lastname (Org)"
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
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Treat")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    reloadFormLayout.addRow(self.applyButton)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
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
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputFiducialsNodeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

  #def enableOrDisableCreateButton(self):
  #  """Connected to both the fiducial and camera node selector. It allows to
  #  enable or disable the 'create path' button."""
  #  self.createPathButton.enabled = self.cameraNodeSelector.currentNode() != None and self.inputFiducialsNodeSelector.currentNode() != None

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputFiducialsNodeSelector.currentNode() and self.outputSelector.currentNode()

  def onApplyButton(self):
    logic = PyLITTPlanLogic()
    screenshotScaleFactor = int(self.PerfusionValueSliderWidget.value)
    print("Run the algorithm")
    fiducialsNode = self.inputFiducialsNodeSelector.currentNode();
    if fiducialsNode.GetClassName() == "vtkMRMLMarkupsFiducialNode":
      # slicer4 Markups node
      NumFid = fiducialsNode.GetNumberOfFiducials()
      if NumFid < 2:
        print("Two Fiducials Needed ")
        return
      # get fiducial positions
      for iFid in xrange(NumFid):
        coord = [0.0, 0.0, 0.0]
        fiducialsNode.GetNthFiducialPosition(iFid, coord)
        print coord
    print  "Model Params:", self.PowerValueSliderWidget.value, self.AbsorptionValueSliderWidget.value 
    # TODO logic.run for brainNek
    #logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), enableScreenshotsFlag,screenshotScaleFactor)
    x = slicer.util.loadModel("/Users/fuentes/fem.stl",True)
    print x[1].GetName() 

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
