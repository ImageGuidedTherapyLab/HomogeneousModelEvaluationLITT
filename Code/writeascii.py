import vtk
import vtk.util.numpy_support as vtkNumPy
import numpy
import scipy.io as scipyio
input_filename = 'temperature.vtk'
extension = input_filename.split('.').pop()
vtkReader = None
if extension == 'vtk':
  vtkReader = vtk.vtkDataSetReader()
elif extension == 'vti':
  vtkReader = vtk.vtkXMLImageDataReader()
else:
  raise RuntimeError('unknown file type %s ' % input_filename)
vtkReader.SetFileName( "%s" % (input_filename) )
vtkReader.Update()
# get image info
imageDataVTK = vtkReader.GetOutput()
dimensions = imageDataVTK.GetDimensions()
spacing = imageDataVTK.GetSpacing()
origin  = imageDataVTK.GetOrigin()
# get data pointer
image_point_data = imageDataVTK.GetPointData()
image_data       = vtkNumPy.vtk_to_numpy( image_point_data.GetArray(0) )
 
print spacing, origin, dimensions
numpy.savetxt('temperature.txt',image_data)
