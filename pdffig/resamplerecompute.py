#! /usr/bin/env python
# $Id: demo.py 299 2007-03-30 12:52:17Z mhagger $

# Copyright (C) 1999-2003 Michael Haggerty <mhagger@alum.mit.edu>
#
# This file is licensed under the GNU Lesser General Public License
# (LGPL).  See LICENSE.txt for details.

"""demo.py -- Demonstrate the Gnuplot python module.

Run this demo by typing 'python demo.py'.  For a more complete test of
the Gnuplot package, see test.py.

"""

from numpy import *

# If the package has been installed correctly, this should work:
import Gnuplot, Gnuplot.funcutils,os

def TransferTimePlot(outputName):
    """Demonstrate the Gnuplot package."""

    # A straightforward use of gnuplot.  The `debug=1' switch is used
    # in these examples so that the commands that are sent to gnuplot
    # are also output on stderr.
    g = Gnuplot.Gnuplot(debug=1)
    #g.title('A simple example') # (optional)
    g('set term pslatex color size 7.5,5.25')
    g('set   output "%s.tex"' % outputName )
    #g.ylabel('blood perfusion $\\\\frac{kg}{s\\\\;m^3}$') # (optional)

    # initialize data dictionary
    data = {}

    # number of nodes in typical problem
    NumberOfNodes = [10000, 100000,1000000]
    # assume it takes 27 column entries per node/row  in the stiffness matrix 
    DataPerNode   = 27. * 8.e0 # assume doubles
    # data size in Bytes
    DataSizeList  = map(lambda nodenum: nodenum*DataPerNode, NumberOfNodes)
    DataSizeList  = array( DataSizeList , dtype='float_') 

    # bandwidth in Gb/s
    BandwidthTransfer = {'global':(5.e9,1), 'local':(200.e9,2), 'private':(500.e9,3)}
    for key,(bandwidth,color) in BandwidthTransfer.iteritems():
      print key,bandwidth,color
      transfertime = map(lambda datasize: datasize/bandwidth, DataSizeList)
      data[key]  = Gnuplot.Data(NumberOfNodes, transfertime , title=key , with_='l lw 7 lc %d lt 1 ' % color )

    # assume it takes 1000 floating point operations to assemble each row in the stiffness matrix 
    FlopPerNode   = 1000
    # total floating point operations to generate matrix
    FlopList  = map(lambda nodenum: nodenum*FlopPerNode , NumberOfNodes)
    FlopList  = array( FlopList , dtype='float_') 

    # Number GPU Cores on Kepler = 15 SMX units * 192 cors per SMX
    # http://www.nvidia.com/content/PDF/kepler/NV_DS_TeslaK_Family_May_2012_LR.pdf
    # format 'arch':(cores,ClockSpeed,FlopPerCycle,plottingcolor)
    ArchitectureSpecs   = {'fermi':(8*128,1.4e9,2,1), 'kepler':(15*192,1.4e9,2,2), 'MIC':(64,1.4e9,4,3),'CPU':(20,1.4e9,4,4)}

    for key,(cores,clockspeed,flopspercycle,color) in ArchitectureSpecs.iteritems():
      print key,cores,color
      # Max flops per second
      MaxFlops = flopspercycle * clockspeed * cores
      computetime = map(lambda floatops: floatops/MaxFlops, FlopList)
      data[key]  = Gnuplot.Data(NumberOfNodes, computetime, title=key , with_='l lw 7 lc %d lt 2 ' % color )

    g('set format x "%9.2e"' ) 
    g('set format y "%9.2e"' ) 
    g('set logscale xy' ) 
    g.xlabel('num nodes ') # (optional)
    g.ylabel('time $(s)$') # (optional)
    g('set key spacing 1.4 at 4.0e4,5.e-2')
    g('set label "{\\\\normalsize  compute time = $\\\\frac{FLOP}{FLOPS_{max}}$}" at 4.e4,5.e-2 ')
    g('set label "{\\\\normalsize transfer time = $\\\\frac{Data}{Bandwidth}$}" at 4.e5,5.e-2 ')
    g.plot(data['global'],data['local'],data['private'],data['fermi'],data['kepler'],data['MIC'],data['CPU'])

    # give gnuplot an arbitrary commands
    #g('set xr [20.:500.]' ) 
    #g('set grid') 
    #g('set yr [-1.:10.0]') 
    #g('set y2r [450:650]' ) 
    ## g('set y2tics  (40,50,60,70,80) ' ) 
    #g('set y2label  "optical absorption $\\\\frac{1}{m}$" ' ) # (optional)
    ## xtics = ','.join( [ '"%d" %f' % (time,time) for time in acquisitionTime] )
    ## g('set xtics (%s)' % xtics ) 
    # output to pdf
    #g('set term pslatex color size 7.5,5.25')
    #g('set   output "powerProfile.tex"         ')
    ## g('set key spacing 1.7 at 130,15')
    ## g('set label "{\\\\normalsize (2) nanoparticle mediated heating}"   at 370,12.1 ')
    ## g('set label "{\\\\small$\\\\delta t_{heating} = [204,390] $}" at 300,7.7 ')
    ## g('set label "{\\\\large $\\\\omega=\\\\omega^\\\\text{native} + \\\\frac{\\\\Omega}{\\\\ln{2}+\\\\Omega}(\\\\omega^\\\\text{coag} - \\\\omega^\\\\text{native})  $}" at 120.0,5.2')
    ## g('set label "{\\\\large $\\\\mu_a=\\\\mu_a^\\\\text{native} + \\\\frac{\\\\Omega}{\\\\ln{2}+\\\\Omega}(\\\\mu_a^\\\\text{coag} - \\\\mu_a^\\\\text{native})  $}" at 352.0,6.5')

    ## deltat=6.0
    ## temperature = arange(20,100,deltat, dtype='float_')
    ## temperature = [ 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 37.7, 43.3, 47.4, 50.3, 52.5, 53.3, 49.1, 46.2, 44.2, 42.8, 41.8, 41.0, 40.4, 39.9, 39.5, 39.1, 38.9, 38.6, 39.4, 45.8, 50.4, 53.6, 56.1, 58.0, 59.5, 60.7, 61.8, 62.6, 63.4, 64.0, 64.6, 65.0, 65.5, 65.8, 66.2, 66.5, 66.7, 66.9, 67.1, 67.3, 67.5, 67.6, 67.8, 67.9, 68.0, 68.1, 68.2, 68.3, 68.3, 67.4, 60.9, 56.2, 52.9, 50.4, 48.4, 46.9, 45.6, 44.5, 43.6, 42.8, 42.2, 41.6, 41.1, 40.7, 40.3, 39.9, 39.6, 39.4, 39.1]
    ## temperature = [ 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 36.9, 45.8, 50.4, 53.6, 56.1, 58.0, 59.5, 60.7, 61.8, 62.6, 63.4, 64.0, 64.6, 65.0, 65.5, 65.8, 66.2, 66.5, 66.7, 66.9, 67.1, 67.3, 67.5, 67.6, 67.8, 67.9, 68.0, 68.1, 68.2, 68.3, 68.3, 67.4, 60.9, 56.2, 52.9, 50.4, 48.4, 46.9, 45.6, 44.5, 43.6, 42.8, 42.2, 41.6, 41.1, 40.7, 40.3, 39.9, 39.6, 39.4, 39.1]
    ## npoints= len(temperature)
    ## time        = deltat * arange(npoints, dtype='float_')
    ## damageFrac  = ArrheniusDamageFraction(temperature,deltat)
    ## w_0  = 6.0
    ## w_1  = 1.0
    ## w_lb = 0.0
    ## w_ub = 8.0
    ## g('set ytics  ( "$\\\\omega^\\\\text{linear}_{lb}$" %f , "$\\\\omega^\\\\text{coag}$" %f , "$\\\\omega^\\\\text{native}$" %f, "$\\\\omega^\\\\text{linear}_{ub}$" %f ) nomirror ' % (w_lb,w_1,w_0,w_ub) ) 
    ## perfusion   =  map( lambda dmge : (w_0 + dmge*(w_1-w_0)), damageFrac )
    ## mu_a_0  = 500.0
    ## mu_a_1  = 600.0
    ## mu_a_lb = 475.
    ## mu_a_ub = 625.
    ## g('set y2tics  ( "$\\\\mu^\\\\text{linear}_{lb}$" %f ,"$\\\\mu_a^\\\\text{native}$" %f ,"$\\\\mu_a^\\\\text{coag}$" %f , "$\\\\mu^\\\\text{linear}_{ub}$" %f ) nomirror ' % (mu_a_lb ,mu_a_0,mu_a_1,mu_a_ub) ) 
    ## absorption  =  map( lambda dmge : (mu_a_0 + dmge*(mu_a_1-mu_a_0)), damageFrac )
    #data['local'] = Gnuplot.Data(time,            damageFrac, title='damage'     , with_='l lw 7 lc 2 lt 2 ',axes="x1y2")

    # Plot a list of (x, y) pairs (tuples or a numpy array would
    # also be OK):
    #g.plot([[0,1.1], [1,5.8], [2,3.3], [3,4.2]])
    ## g('set multiplot')
    ## g('set key spacing 1.2 at 160,9.7')
    ## # set margins
    ## g('set rmargin 9' ) # (optional)
    ## g('set lmargin 9' ) # (optional)
    ## g('set bmargin 0 ') # (optional)
    ## # size of this plot
    ## g('set size   1.0,0.7')
    ## g('set origin 0.0,0.3')
    ## g('set xtics   ( "" 50, "" 100, "" 150 , "" 200 , "" 250, "" 300 , "" 350, "" 400, "" 450 , "" 500 ) '  ) 
    ## # size of this plot
    ## g('set origin 0.0,0.0')
    ## g('set size   1.0,0.3')
    ## g('set key spacing 1.2 at 160,75')
    ## g('set bmargin   ' ) # (optional)
    ## g('set tmargin 0 ' ) # (optional)
    ## g('unset label ')
    ## g('set label "{damage=$\\\\frac{\\\\Omega}{\\\\ln{2}+\\\\Omega} $}" at 100.0,50.4')
    ## g('set label "{$\\\\Omega(t) = \\\\int_0^t \\\\! Ae^{\\\\frac{-E_A}{Ru(\\\\tau)}} \\\\, \\\\mathrm{d} \\\\tau$}" at 320.0,50.4')
    ## g('set yr  [35:80]' ) 
    ## g('set y2r [-0.2:1.2]') 
    ## g('set xtics   ( "" 50, "" 100, "" 150 , "0" 200 , "50" 250, "100" 300 , "150" 350, "200" 400, "250" 450 , "300" 500 ) '  ) 
    ## g('set ytics   ( "37" 37, "43" 43 , "51" 51, "61" 61, "75" 75 ) ' ) 
    ## g('set y2tics  (  "0"  0, "1"   1 ) ' ) 
    ## g('set ylabel   "temperature $^oC$"' ) # (optional)
    ## g('set y2label  "damage"' ) # (optional)
    ## g.xlabel('time  $(s)$ ') # (optional)
    ## g.plot(data['tempHistory'],data['damageFrac'])
    ## g('unset multiplot')
    #g.plot(data['damage'],data['(1) MRTI'],data['(2) MRTI'],data['(1) FEM'],data['(2) FEM'])
    g('unset output ') # flush the output buffer
    os.system("""latex -jobname %s                                     \
\\\\documentclass{article}                                             \
\\\\usepackage{amssymb,amsfonts,amsmath}                               \
\\\\usepackage{nopageno}                                               \
\\\\usepackage[left=.5in,right=.5in,top=1.0in,bottom=.5in]{geometry}   \
\\\\begin{document}                                                    \
\\\\LARGE                                                              \
\\\\begin{center}                                                      \
\\\\input{%s}                                                          \
\\\\end{center}                                                        \
\\\\end{document}                                                      \
    """ % (outputName,outputName) )
    os.system("dvipdf  %s.dvi" % outputName)
    os.system("pdfcrop %s.pdf" % outputName)
    os.system("mv      %s-crop.pdf %s.pdf" % (outputName,outputName) )

# when executed, just run powerProfile():
if __name__ == '__main__':
    TransferTimePlot("resamplerecompute")
