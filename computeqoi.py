# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.
# DAKOTA will execute this script 

# necessary python modules
import sys
import re
import os
import scipy.io as scipyio
import ConfigParser
import time
import pyopencl as cl
import numpy
import numpy.linalg as la

brainNekDIR     = '/workarea/fuentes/braincode/tym1' 
workDirectory   = 'optpp_pds'
outputDirectory = '/tmp/outputs/dakota/%04d'

# database and run directory have the same structure
databaseDIR     = 'database/'

##################################################################
class BrainNekWrapper:
  def __init__(self, SEMDataDirectory,variableDictionary  ):
     
     self.DataDictionary = {}
     self.DebugObjective = True
     self.DebugObjective = False
     self.ctx   = cl.create_some_context()
     self.queue = cl.CommandQueue(self.ctx)
     self.prg   = cl.Program(self.ctx, """
              __kernel void diff_sq(__global const float *a,
              __global const float *b, __global float *c)
              {
                int gid = get_global_id(0);
                c[gid] = (a[gid] - b[gid]) * (a[gid] - b[gid]);
              }
              """).build()

     # FIXME  should this be different ?  
     self.SEMDataDirectory = SEMDataDirectory 
  
     # FIXME vtk needs to be loaded AFTER kernel is built
     import vtk
     import vtk.util.numpy_support as vtkNumPy 
     print "using vtk version", vtk.vtkVersion.GetVTKVersion()
     print "read SEM data"

     start = time.clock()
     vtkSEMReader = vtk.vtkXMLUnstructuredGridReader()
     vtufileName = "%s/%d.vtu" % (self.SEMDataDirectory,0)
     vtkSEMReader.SetFileName( vtufileName )
     vtkSEMReader.SetPointArrayStatus("Temperature",1)
     vtkSEMReader.Update()
     elapsed = (time.clock() - start)
     print "read SEM data", elapsed
  
     # get registration parameters
  
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
     self.SEMRegister = vtk.vtkTransformFilter()
     self.SEMRegister.SetInput(vtkSEMReader.GetOutput())
     self.SEMRegister.SetTransform(AffineTransform)
     self.SEMRegister.Update()
  
     print "write transform output"
     if ( self.DebugObjective ):
        vtkSEMWriter = vtk.vtkDataSetWriter()
        vtkSEMWriter.SetFileTypeToBinary()
        semfileName = "%s/semtransform%04d.vtk" % (self.SEMDataDirectory,0)
        print "writing ", semfileName 
        vtkSEMWriter.SetFileName( semfileName )
        vtkSEMWriter.SetInput(self.SEMRegister.GetOutput())
        vtkSEMWriter.Update()



  ##################################################################
  def ComputeObjective(self,MRTIDataDirectory,VolumeOfInterest ):
    print self.SEMDataDirectory 
    ObjectiveFunction = 0.0
  
    # loop over time points of interest
    for idwrite,(SEMtimeID,MRTItimeID) in enumerate([(0,135),(0,136),(0,137),(0,138)]):
  
      mrtifilename = '%s/temperature.%04d.vtk' % (MRTIDataDirectory,MRTItimeID) 
      if MRTItimeID in self.DataDictionary:
        print "already loaded",  mrtifilename 
      else:
        print 'opening' , mrtifilename 
        # FIXME vtk needs to be loaded AFTER kernel is built
        import vtk
        import vtk.util.numpy_support as vtkNumPy 
        print "using vtk version", vtk.vtkVersion.GetVTKVersion()
        print "read SEM data"
        vtkImageReader = vtk.vtkDataSetReader() 
        vtkImageReader.SetFileName(mrtifilename )
        vtkImageReader.Update() 
        ## image_cells = vtkImageReader.GetOutput().GetPointData() 
        ## data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars')) 
        
        # extract voi for QOI
        vtkVOIExtract = vtk.vtkExtractVOI() 
        vtkVOIExtract.SetInput( vtkImageReader.GetOutput() ) 
        vtkVOIExtract.SetVOI( VolumeOfInterest  ) 
        vtkVOIExtract.Update()
        mrti_point_data= vtkVOIExtract.GetOutput().GetPointData() 
        mrti_array = vtkNumPy.vtk_to_numpy(mrti_point_data.GetArray('image_data')) 
        #print mrti_array
        #print type(mrti_array)
  
        print "project SEM onto MRTI for comparison"
        vtkResample = vtk.vtkCompositeDataProbeFilter()
        vtkResample.SetSource( vtkVOIExtract.GetOutput() )
        vtkResample.SetInput( self.SEMRegister.GetOutput() ) 
        vtkResample.Update()
  
        if ( self.DebugObjective ):
          roifileName = "%s/mrti.%04d.vtu" % (self.SEMDataDirectory,MRTItimeID)
          print "writing ", roifileName 
          start = time.clock()
          # setup mrti projection file
          vtkTemperatureWriter = vtk.vtkXMLUnstructuredGridWriter()
          vtkTemperatureWriter.SetFileName( roifileName )
          vtkTemperatureWriter.SetInput(vtkResample.GetOutput())
          vtkTemperatureWriter.Update()
          elapsed = (time.clock() - start)
          print "write output", elapsed
  
        print " accumulate objective function"
        fem_point_data= vtkResample.GetOutput().GetPointData() 
        self.DataDictionary[MRTItimeID] = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('image_data')) 
        #print fem_array 
        #print type(fem_array )

      h_mrti = self.DataDictionary[MRTItimeID] 
      mf = cl.mem_flags
      d_mrti = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=h_mrti )
      # TODO need to get cl buffer from brainNek
      d_fem  = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=h_mrti )

      # storage for output
      dest_buf = cl.Buffer(self.ctx, mf.WRITE_ONLY, h_mrti.nbytes)

      self.prg.diff_sq(self.queue, h_mrti.shape, None, d_mrti , d_fem , dest_buf)

      mrti_minus_fem_squared = numpy.empty_like(h_mrti)
      cl.enqueue_copy(self.queue, mrti_minus_fem_squared , dest_buf)
      ObjectiveFunction = ObjectiveFunction + mrti_minus_fem_squared.sum()

    return ObjectiveFunction 
  # end def ComputeObjective:
##################################################################
def ParseInput(paramfilename):
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

  DescriptorList = ['robin_coeff','probe_init','anfact_healthy', 'mu_a_healthy','mu_s_healthy','k_0_healthy','w_0_healthy','x_displace','y_displace','z_displace','x_rotate','y_rotate','z_rotate']
  for paramname in DescriptorList:
    try:
      continuous_vars[paramname  ] = paramsdict[paramname ]
    except KeyError:
      pass
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # store dakota vars
  fem_params['cv']         = continuous_vars
  fem_params['asv']        = active_set_vector
  fem_params['functions']  = num_fns
  fem_params['fileID']     = fileID 
  fem_params['UID']        = int(paramfilename.split('/').pop(3))

  # parse file path
  locatemrti = paramfilename.split('/')
  locatemrti.pop()

  # database and run directory have the same structure
  fem_params['mrti']       = '%s/%s/%s/vtk/referenceBased/' % (databaseDIR,locatemrti[2],locatemrti[3])

  # get power file name
  inisetupfile  = "/".join(locatemrti)+"/setup.ini"
  config = ConfigParser.SafeConfigParser({})
  config.read(inisetupfile)
  fem_params['ccode']        = config.get('power','ccode')
  fem_params['semwritetime'] = config.getfloat('mrti','deltat') * config.getfloat('mrti','maxheat')
  fem_params['maxheat']      = config.getfloat('mrti','maxheat')
  fem_params['voi']          = eval(config.get('mrti','voi'))

  print 'mrti data from' , fem_params['mrti'] , 'setupfile', inisetupfile  

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
(options, args) = parser.parse_args()



if (options.param_file != None):
  import brainNekLibrary

  # parse the dakota input file
  fem_params = ParseInput(options.param_file)

  brain = BrainNekWrapper(outputDirectory % fem_params['UID'],fem_params['cv'])

  setup = brainNekLibrary.PySetupAide("optpp_pds/setuprc.%04d" %  fem_params['fileID'] )
  print setup 
  brainNek = brainNekLibrary.PyBrain3d(setup);
  print 'intpointer', brainNek.getTemperaturePointer() 
  tstep = 0
  ## while( brainNek.timeStep(tstep * .25 ) ) :
  ##   print brainNek.dt
  ##   tstep = tstep + 1
  ## # write objective function back to Dakota
  ## objfunction = brain.ComputeObjective(fem_params['mrti'],fem_params['voi'])
  ## print "objective function 1", objfunction
  ## objfunction = brain.ComputeObjective(fem_params['mrti'],fem_params['voi'])
  ## print "objective function 2", objfunction
  ## objfunction = brain.ComputeObjective(fem_params['mrti'],fem_params['voi'])
  ## print "objective function 3", objfunction

else:
  parser.print_help()
  print options
