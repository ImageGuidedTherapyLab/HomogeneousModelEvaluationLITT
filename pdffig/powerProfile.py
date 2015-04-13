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

def PowerProfile(y):
    x = y/6.
    if(x<29):
     return 0
    elif(x<34):
     return 3.0
    elif(x<55):
     return 0
    elif(x<60):
     return 12.5
    else:
     return 0

def PlotPower(outputName):
    """Demonstrate the Gnuplot package."""

    # A straightforward use of gnuplot.  The `debug=1' switch is used
    # in these examples so that the commands that are sent to gnuplot
    # are also output on stderr.
    g = Gnuplot.Gnuplot(debug=1)
    #g.title('A simple example') # (optional)
    g('set term pslatex color size 7.5,5.25')
    g('set   output "%s.tex"' % outputName )
    g.xlabel('time  [s]') # (optional)
    g.ylabel('power [W]') # (optional)

    dataFiles = [["profiles/temperature0.0.csv","(2) FEM" ,(2,0),("x1y2",2,1)],
                 ["profiles/temperature1.0.csv","(1) FEM" ,(2,0),("x1y2",3,1)]
                ]
    deltat=6.0
    data = {}
    for filename,varName,(xcomp,ycomp),(axisplot,lc,lt)  in dataFiles:
       csvData = loadtxt(filename,delimiter=',',skiprows=1)
       data[varName] = Gnuplot.Data(deltat * csvData [:,xcomp] , csvData [:,ycomp],
                        title=varName,
                        with_='l lw 7 lc %d lt %d ' % (lc,lt) ,axes=axisplot )
    acquisitionImageID = [0,16,21,34,65,84]
    acquisitionImageID = [0,21,34,65,84]
    acquisitionTime =  map( lambda t : t*deltat, acquisitionImageID )
    maxTime = acquisitionTime[-1] 
    # give gnuplot an arbitrary commands
    g('set xr [*:%f]' % (maxTime)) 
    g('set grid') 
    g('set yr [*:25.0]') 
    g('set ytics  (3,5,10,15) nomirror ' ) 
    g('set y2r [30:88]' ) 
    g('set y2tics  (40,50,60,70,80) ' ) 
    g('set y2label  "temperature [$^o$C]"' ) # (optional)
    xtics = ','.join( [ '"%d" %f' % (time,time) for time in acquisitionTime] )
    g('set xtics (%s)' % xtics ) 
    # output to pdf
    #g('set term pslatex color size 7.5,5.25')
    #g('set   output "powerProfile.tex"         ')
    g('set key spacing 1.7 at 130,15')
    #g('set label "{\\\\normalsize (1) normal heating}" at 370,3.9 ')
    #g('set label "{\\\\normalsize (2) nanoparticle mediated heating}"   at 370,12.1 ')
    #g('set label "{\\\\small$\\\\delta t_{heating} = [204,390] $}" at 300,7.7 ')
    #g('set label "{\\\\small$\\\\delta t_{cooling} = [390,504] $}" at 450,6.7 ')


    time = arange(0,maxTime,5, dtype='float_')
    data['power'] = Gnuplot.Data(time, map(PowerProfile,time), title='power', with_='boxes fs solid 0.2 lc 1 ')

    # Plot a list of (x, y) pairs (tuples or a numpy array would
    # also be OK):
    #g.plot([[0,1.1], [1,5.8], [2,3.3], [3,4.2]])
    g.plot(data['power'],data['(1) FEM'],data['(2) FEM'])
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
    ## raw_input('Please press return to continue...\n')

    ## g.reset()
    ## # Plot one dataset from an array and one via a gnuplot function;
    ## # also demonstrate the use of item-specific options:
    ## x = arange(10, dtype='float_')
    ## y1 = x**2
    ## # Notice how this plotitem is created here but used later?  This
    ## # is convenient if the same dataset has to be plotted multiple
    ## # times.  It is also more efficient because the data need only be
    ## # written to a temporary file once.
    ## d = Gnuplot.Data(x, y1,
    ##                  title='calculated by python',
    ##                  with_='points 3 3')
    ## g.title('Data can be computed by python or gnuplot')
    ## g.xlabel('x')
    ## g.ylabel('x squared')
    ## # Plot a function alongside the Data PlotItem defined above:
    ## g.plot(Gnuplot.Func('x**2', title='calculated by gnuplot'), d)
    ## raw_input('Please press return to continue...\n')

    ## # Save what we just plotted as a color postscript file.

    ## # With the enhanced postscript option, it is possible to show `x
    ## # squared' with a superscript (plus much, much more; see `help set
    ## # term postscript' in the gnuplot docs).  If your gnuplot doesn't
    ## # support enhanced mode, set `enhanced=0' below.
    ## g.ylabel('x^2') # take advantage of enhanced postscript mode
    ## g.hardcopy('gp_test.ps', enhanced=1, color=1)
    ## print ('\n******** Saved plot to postscript file "gp_test.ps" ********\n')
    ## raw_input('Please press return to continue...\n')

    ## g.reset()
    ## # Demonstrate a 3-d plot:
    ## # set up x and y values at which the function will be tabulated:
    ## x = arange(35)/2.0
    ## y = arange(30)/10.0 - 1.5
    ## # Make a 2-d array containing a function of x and y.  First create
    ## # xm and ym which contain the x and y values in a matrix form that
    ## # can be `broadcast' into a matrix of the appropriate shape:
    ## xm = x[:,newaxis]
    ## ym = y[newaxis,:]
    ## m = (sin(xm) + 0.1*xm) - ym**2
    ## g('set parametric')
    ## g('set data style lines')
    ## g('set hidden')
    ## g('set contour base')
    ## g.title('An example of a surface plot')
    ## g.xlabel('x')
    ## g.ylabel('y')
    ## # The `binary=1' option would cause communication with gnuplot to
    ## # be in binary format, which is considerably faster and uses less
    ## # disk space.  (This only works with the splot command due to
    ## # limitations of gnuplot.)  `binary=1' is the default, but here we
    ## # disable binary because older versions of gnuplot don't allow
    ## # binary data.  Change this to `binary=1' (or omit the binary
    ## # option) to get the advantage of binary format.
    ## g.splot(Gnuplot.GridData(m,x,y, binary=0))
    ## raw_input('Please press return to continue...\n')

    ## # plot another function, but letting GridFunc tabulate its values
    ## # automatically.  f could also be a lambda or a global function:
    ## def f(x,y):
    ##     return 1.0 / (1 + 0.01 * x**2 + 0.5 * y**2)

    ## g.splot(Gnuplot.funcutils.compute_GridData(x,y, f, binary=0))
    ## raw_input('Please press return to continue...\n')

    ## # Explicit delete shouldn't be necessary, but if you are having
    ## # trouble with temporary files being left behind, try uncommenting
    ## # the following:
    ## #del g, d


# when executed, just run powerProfile():
if __name__ == '__main__':
    PlotPower("powerProfile")

