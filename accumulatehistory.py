import numpy
import os
import ConfigParser

resultfileList = [
'./workdir/Study0035/0530/',
'./workdir/Study0030/0495/',
##'./workdir/Study0023/0433/',
'./workdir/Study0030/0497/',
'./workdir/Study0030/0491/',
'./workdir/Study0030/0496/',
'./workdir/Study0030/0490/',
'./workdir/Study0017/0378/',
'./workdir/Study0025/0438/',
'./workdir/Study0025/0435/',
'./workdir/Study0025/0440/',
'./workdir/Study0025/0436/',
'./workdir/Study0028/0466/',
'./workdir/Study0028/0468/',
'./workdir/Study0028/0471/',
'./workdir/Study0026/0447/',
'./workdir/Study0026/0457/',
'./workdir/Study0026/0455/',
'./workdir/Study0026/0453/',
'./workdir/Study0026/0450/',
'./workdir/Study0026/0451/',
'./workdir/Study0022/0418/',
'./workdir/Study0022/0417/',
'./workdir/Study0021/0409/',
'./workdir/Study0021/0414/',
'./workdir/Study0021/0415/',
]

#resultfileList = [
#'./workdir/Study0035/0530/',
#]

outputDirectory = '/tmp/outputs/dakota/%04d'

def DiceTxtFileParse(DiceInputFilename):
  # (1) split on ':' (2)  filter lists > 1 (3) convert to dictionary
  c3doutput = dict(filter( lambda x: len(x) > 1,[line.strip().split(':') for line in open(DiceInputFilename) ] ))
  return float(c3doutput['Dice similarity coefficient'])

with file('datasummary.tex' , 'w') as texHandle: 
  with file('datasummary.txt' , 'w') as fileHandle: 
    # write header
    fileHandle.write("iddata,mu_eff,obj\n")
    # loop over files and extract optimal value
    opttype = 'heating'
    for filenamebase in resultfileList:
      # get latex command
      config = ConfigParser.SafeConfigParser({})
      inisetupfile = '%s/opt/setup.ini' % (filenamebase)
      config.read(inisetupfile)
  
      filename = '%s/opt/dakota_q_newton_%s.in.tabular.dat' % (filenamebase,opttype)
      dataid = int(filename.split('/')[3])
      mu_eff = numpy.loadtxt(filename,skiprows=1,usecols=(1,))
      obj_fn = numpy.loadtxt(filename,skiprows=1,usecols=(2,))
      if( len(obj_fn.shape) > 0 ) :
        if( obj_fn.shape[0]  > 1 ) :
          idmin = numpy.argmin(obj_fn)
          print idmin, mu_eff[idmin], obj_fn[idmin]
          mu_effopt =  mu_eff[idmin]
          minobjval =  obj_fn[idmin]
      else:
        print "single entry ??", mu_eff, obj_fn
        mu_effopt =  mu_eff
        minobjval =  obj_fn 
      #dataarray = numpy.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4,6)
      fileHandle.write("%05d,%12.5e,%12.5e\n" %(dataid,mu_effopt,minobjval))
      # FIXME
      runcmd = "vglrun python ./brainsearch.py --param_file  %s/opt/optpp_pds.%s.in.%d %s/opt/optpp_pds.%s.out.%d --vis_out" % (filenamebase,opttype,idmin,filenamebase,opttype,idmin)
      print runcmd
      os.system( runcmd )
  
      # get arrhenius dice value
      heattimeinterval               = eval(config.get('mrti','heating')  )
      SEMDataDirectory               = outputDirectory % int(filenamebase.split('/')[-2]) 
      dicefilename = "%s/dice.%s.%04d.txt" % (SEMDataDirectory,opttype,heattimeinterval[1])
      print dicefilename 
      dicevalue = DiceTxtFileParse(dicefilename)
  
      # format latex ouput
      #outputformat                   = config.get('latex',opttype)
      #texFormat = outputformat % (minobjval,dicevalue)
      #print texFormat 
      #texHandle.write("%s\n" %(texFormat))
