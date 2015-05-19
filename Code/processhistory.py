
import re

rundirectory = 'workdir/Patient0006/007/opt/'
maxfile = 126

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

  DescriptorList = ['robin_coeff','probe_init','mu_eff_healthy','anfact_healthy', 'mu_a_healthy','mu_s_healthy','k_0_healthy','w_0_healthy','x_displace','y_displace','z_displace','x_rotate','y_rotate','z_rotate']
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

  return fem_params

with file('%s/iterationhistory.txt' % rundirectory, 'w') as fileHandle: 
  fileHandle.write("iter,k,mu_eff,dz,u0,h,obj\n")
  for idfile in  range(1,maxfile):
   infile =  '%s/optpp_pds.in.%d'  % (rundirectory,idfile) 
   outfile=  '%s/optpp_pds.out.%d' % (rundirectory,idfile)
   fem_params = ParseInput(infile )
   cv =  fem_params['cv']
   with file(outfile, 'r') as objectivefnc: 
     print 
     fileHandle.write("%05d,%s,%s,%s,%s,%s,%s" %(idfile,cv['k_0_healthy'],cv['mu_eff_healthy'],cv['z_displace'],cv['probe_init'],cv['robin_coeff'],objectivefnc.read()))

