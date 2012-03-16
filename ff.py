#!/usr/bin/env python
"""
This is a simple finite field helper script. By default it prepares a set of
input files for finite field calculations of each xyz structure provided as an
input, in which case an input template is required. If none is found a standard
one is prepared. Invoking the script with -c option allows to retrieve the
properties, provided that the calculations are completed, in which case paths
to data dir(s) are expected on input.

Usage: ff.py [options] xyz file(s) or data dir(s)

Options:
  -h, --help       show this help
  -c, --calculate  calculate properties for each data dir
  -f, --base-field set base field ( >= 0.0001; defaults to 0.001 )
  -s, --gaussian   prepare gaussian inputs (expects gaussian.tmpl)
  -g, --gamess     prepare gamess inputs (expects gamess.tmpl)
  -m, --molcas     prepare molcas inputs (expects molcas.tmpl)
  -u, --units      set units; chose from: esu, si, asi 
                   (si multiplied by electric permittivity of free space),
                   or au which is the default
  --print-data     if '-c' option is set, the script prints also the extracted
                   energies and dipole moments (usefull for testing purposes)

Fields selection:
  --fields 25         25 fields sufficient for <g> (5pt FF formulas, enabled by default)
           13         13 fields sufficient for diagonal components (5pt)
           5x          5 fields along x (y or z) axis for selected longitudinal components (5pt)
           55         55 fields sufficient for <g> (7pt FF formulas)
           19         19 fields sufficient for diagonal components (7pt)
           7x          7 fields along x (y or z) axis for selected longitudinal components (7pt)
           GR[n,a,x]  Use for Generalized Romberg: n - the number of iterations
                                                   a - the quotient (e.g. 2 for RR)
                                                   x - the axis (0,1,2 for x,y or z)

  --print-fields      print the chosen selection of fields

"""

#     Copyright (C) 2011, Robert W. Gora (robert.gora@pwr.wroc.pl)
#  
#     This program is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#  
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#  
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

__author__  = "Robert Gora (robert.gora@pwr.wroc.pl)"
__version__ = filter(str.isdigit, "$Revision$")

# Import necessary modules
import os, sys, getopt, re

from string import Template
from numpy import *

# Regular expressions
reflags = re.DOTALL

#----------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------
def Usage():
    """Print usage information and exit."""
    print __doc__
    print "Machine epsilon is: ",finfo(float64).eps,"for float64 type\n"

    sys.exit()

#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------
def Main(argv):
    '''Parse commandline and loop throught the logs'''

    data = {}

    # Set up defaults
    fstep = 0.001
    calculate = 0
    density = 0
    show = 0
    gamess = 0
    gaussian = 0
    molcas = 0
    native = 0
    frange = range(25)
    grr = (0,0,0)
    units = SetUnits('au')

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hu:f:cdgsmn",
                                        ["help",
                                         "units=",
                                         "base-field=",
                                         "fields=",
                                         "print-fields",
                                         "calculate",
                                         "density",
                                         "print-data",
                                         "gamess",
                                         "gaussian",
                                         "molcas",
                                         "molcas-native",
                                         ])
    except getopt.GetoptError, error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
        elif opt in ("-u", "--units="):
            units=SetUnits(arg)
        elif opt in ("-f", "--base-field="):
            fstep = float64(arg)
        elif opt in ("-g", "--gamess"):
            gamess=1
        elif opt in ("-s", "--gaussian"):
            gaussian=1
        elif opt in ("-m", "--molcas"):
            molcas=1
        elif opt in ("-n", "--molcas-native"):
            native=1
        elif opt in ("-c", "--calculate"):
            calculate=1
        elif opt in ("-d", "--density"):
            density=1
        elif opt in ("--fields="):
            if arg == '25':
                frange = range(25)
            elif arg == '13':
                frange = range(13)
            elif arg == '5x':
                frange = [ 0, 1, 2, 7, 8]
            elif arg == '5y':
                frange = [ 0, 3, 4, 9,10]
            elif arg == '5z':
                frange = [ 0, 5, 6,11,12]
            elif arg == '55':
                frange = range(55)
            elif arg == '19':
                frange = range(19)
            elif arg == '7x':
                frange = [ 0, 1, 2, 7, 8,13,14]
            elif arg == '7y':
                frange = [ 0, 3, 4, 9,10,15,16]
            elif arg == '7z':
                frange = [ 0, 5, 6,11,12,17,18]
            elif arg[:2] == 'GR':
                grr = arg[3:-1].split(',')
                grr = ( int(grr[0]), float(grr[1]), int(grr[2]) )
            else:
                print "Unrecognized fields selection ... assuming the default."
        elif opt in ("--print-fields"):
            print INPUTS(None, fstep, frange, grr).PrintFields()
        elif opt in ("--print-data"):
            show=1

    # Parse each data file (with xyz coords) or dir (with the results)
    data_files = args

    for f in data_files:

        # Analyze outputs
        if calculate:

            if gamess:
                data[f]=GAMESS(f, 1, units, 'gamess', density)
            if molcas:
                data[f]=MOLCAS(f, -1, units, 'molcas', density)
            if gaussian:
                data[f]=GAUSSIAN(f, -1, units, 'gaussian', density)
            if show:
                data[f].PrintData()

        # Prepare inputs
        else:
            if gaussian:
                GAUSSIAN_INPUTS(f,fstep,frange,grr)
            if gamess:
                GAMESS_INPUTS(f,fstep,frange,grr)
            if molcas:
                MOLCAS_INPUTS(f,fstep,frange,grr)
            if native:
                MOLCAS_NATIVE_INPUTS(f,fstep,frange,grr)


#----------------------------------------------------------------------------
# Common parser routines
#----------------------------------------------------------------------------
class PARSER:
    """A common parser routines"""

    def __init__(self, logpath, sign, units, pkg, density):
        # initialize
        self.logpath = logpath
        self.sign = sign
        self.units = units
        self.pkg = pkg
        self.density = density
        self.ext = "\.log$"
        # input data
        self.energies = {}
        self.dipoles = {}
        self.densities = {}
        # results
        self.properties = {'E':{}, 'D':{}}
        # useful regexp's
        self.re_number=re.compile(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?")
        # parse files and calculate properties
        self.ParseFiles()
        self.Calculate()

    def ParseFiles(self):
        '''read logfiles in logpath'''
    
        self.fstep = 0

        logfiles=[]
        for f in os.listdir(self.logpath):
            if re.search(self.ext, f):
                logfiles.append(self.logpath+'/'+f)
        logfiles.sort()
    
        for log in logfiles:
            self.ParseFile(log)
    
    def ParseFile(self,filename):
        pass

    def SaveEnergy(self,label,conf_id,conf_no,field,value):
        if label not in self.energies:
            self.energies[label]={}
        if conf_id not in self.energies[label]:
            self.energies[label][conf_id]={}
        if field not in self.energies[label][conf_id]:
            self.energies[label][conf_id][field]=value

    def SaveDipole(self,label,conf_id,conf_no,field,value):
        if label not in self.dipoles:
            self.dipoles[label]={}
        if conf_id not in self.dipoles[label]:
            self.dipoles[label][conf_id]={}
        if field not in self.dipoles[label][conf_id]:
            self.dipoles[label][conf_id][field]=value

    def Calculate(self):
        """Calculate properties"""
        for e in self.energies.keys():

            self.properties['E'][e]={}
            npts = int(2*max(max(self.energies[e].keys()))/self.fstep+1)

            if self.pkg == 'gamess' and self.runtyp == 'eds':
                for i in self.energies[e].keys():
                    self.properties['E'][e][i]={}
                    label = e + " C(" + i + ") energy"
                    if npts >= 5:
                        self.properties['E'][e][i]['FF5'] = FDIFF_E_V(self.energies[e][i], self.fstep, self.units, label+" (5pt)")
                    if npts >= 7:
                        self.properties['E'][e][i]['FF7'] = FDIFF_E_VII(self.energies[e][i], self.fstep, self.units, label+" (7pt)")
                    self.properties['E'][e][i]['RR'] = ROMBERG(self.energies[e][i], -1, label)
            else:
                label = e + " energy"
                if npts >= 5:
                    self.properties['E'][e]['FF5'] = FDIFF_E_V(self.energies[e], self.fstep, self.units, label+" (5pt)")
                if npts >= 7:
                    self.properties['E'][e]['FF7'] = FDIFF_E_VII(self.energies[e], self.fstep, self.units, label+" (7pt)")
                self.properties['E'][e]['RR'] = ROMBERG(self.energies[e], -1, label)

        for d in self.dipoles.keys():

            self.properties['D'][d]={}

            if self.pkg == 'gamess' and self.runtyp == 'eds':
                for i in self.dipoles[d].keys():
                    self.properties['D'][d][i]={}
                    label = d + " C(" + i + ") dipole"
                    self.properties['D'][d][i]['FF'] = FDIFF_D_V(self.dipoles[d][i], self.fstep, self.units, label)
                    self.properties['D'][d][i]['RR'] = {
                        'X': ROMBERG(self.SelectDipoles(self.dipoles[d][i],0), 1, label+" (x component)"),
                        'Y': ROMBERG(self.SelectDipoles(self.dipoles[d][i],1), 1, label+" (y component)"),
                        'Z': ROMBERG(self.SelectDipoles(self.dipoles[d][i],2), 1, label+" (z component)") }
            else:
                label = d + " dipole"
                self.properties['D'][d]['FF'] = FDIFF_D_V(self.dipoles[d], self.fstep, self.units, label)
                self.properties['D'][d]['RR'] = {
                    'X': ROMBERG(self.SelectDipoles(self.dipoles[d],0), 1, label+" (x component)"),
                    'Y': ROMBERG(self.SelectDipoles(self.dipoles[d],1), 1, label+" (y component)"),
                    'Z': ROMBERG(self.SelectDipoles(self.dipoles[d],2), 1, label+" (z component)") }

        if self.pkg == 'gaussian' and self.density == 1:
            FDIFF_DD(self.densities, self.fstep, self.units, self.cubedata)


    def SelectDipoles(self,D,axis):
        """Return selected elements of dipole moments for RR"""

        dipoles={}

        for f in D.keys():
            dipoles[f]=D[f][axis]

        return dipoles

    def SetBaseField(self,field):
        if ((self.fstep == 0 and abs(max(field)) > 0) or
            (abs(max(field)) > 0 and abs(max(field)) < self.fstep)):
            self.fstep = abs(max(field))

    def SortFields(self,fields):

        sfields = sorted( sorted( sorted( fields,
            lambda a,b: cmp(a[2], b[2]) ),
            lambda a,b: cmp(a[0], b[0]) ),
            lambda a,b: cmp(a[1], b[1]) )

        nfields=(len(sfields)-1)/3/2

        return sfields

    def PrintData(self):
        '''Punch extracted energies and dipoles'''

        for e in self.energies.keys():
            print "\n%s energies:" % (e)
            print 55*"-"

            #for f in sorted( self.energies[e].keys(), key=self.energies[e].get ):
            for f in sorted( self.energies[e].keys() ):
                row='E(%7.4f,%7.4f,%7.4f) = ' % f
                row+='%26.15f' % self.energies[e][f]
                print row

#----------------------------------------------------------------------------
# Common input template routines
#----------------------------------------------------------------------------
class INPUT_TEMPLATE(Template):
    delimiter = '@'

class INPUTS:
    """Common input routines"""
    def __init__(self, data, fstep, frange, grr):
        self.data = data
        self.fstep = fstep
        self.grr = grr

        if self.grr[0] > 0:
            self.frange = range(2*self.grr[0]+1)
        else:
            self.frange = frange

        self.ReadTemplate()

    def ReadTemplate(self):
        """Read or punch standard template"""
        try:
            self.tmpl = open(self.pkg+'.tmpl','r').read()
            self.tmpl = INPUT_TEMPLATE(self.tmpl)
            self.WriteInputs()
        except IOError:
            print "There's no " + self.pkg + " template. I'm punching one - please check"
            open(self.pkg+'.tmpl','w').write(self.tmpl)
            sys.exit()
        except AttributeError:
            pass

    def WriteInputs(self):
        pass

    def RRfield(self,p,a,fstep,axis):
        field = [ 0.0, 0.0, 0.0 ]
        field[axis] = a**(p)*fstep
        return tuple(field)

    def Fields(self,i,f):
        '''Set fields for GR or FF calculations'''

        # Generalized Romberg formulas: grr(iter,a,axis)
        if self.grr[0] > 0:

            np = self.grr[0]
            a = self.grr[1]
            axis = self.grr[2]

            fields = [ '%7.4f %7.4f %7.4f' % ( 0.0, 0.0, 0.0 ) ]

            for p in range(np):
                fields.append( '%7.4f %7.4f %7.4f' % self.RRfield(p,a, f,axis) )
            for p in range(np):
                fields.append( '%7.4f %7.4f %7.4f' % self.RRfield(p,a,-f,axis) )

        # 7 point formulas
        elif len(self.frange) == 55 or len(self.frange) == 19 or len(self.frange) == 7 :
            fields=['%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % (-2*f,  0.0,- 0.0),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0, -2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0, -2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  2*f),
                    '%7.4f %7.4f %7.4f'  % (-3*f,  0.0,- 0.0),
                    '%7.4f %7.4f %7.4f'  % ( 3*f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0, -3*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  3*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0, -3*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  3*f),
                    '%7.4f %7.4f %7.4f'  % (  -f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,    f),
                    '%7.4f %7.4f %7.4f'  % (-2*f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (-2*f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (-2*f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % (-2*f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0, -2*f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  2*f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0, -2*f,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  2*f,    f),
                    '%7.4f %7.4f %7.4f'  % (  -f, -2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f, -2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,  2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0, -2*f),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,  2*f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0, -2*f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,  2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f, -2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,  2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f, -2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,  2*f)]

        # 5 point formulas
        else:
            fields=['%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % (-2*f,  0.0,- 0.0),
                    '%7.4f %7.4f %7.4f'  % ( 2*f,  0.0,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0, -2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  2*f,  0.0),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0, -2*f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  2*f),
                    '%7.4f %7.4f %7.4f'  % (  -f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,   -f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (   f,    f,  0.0),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,   -f),
                    '%7.4f %7.4f %7.4f'  % (  -f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % (   f,  0.0,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,   -f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,    f),
                    '%7.4f %7.4f %7.4f'  % ( 0.0,    f,    f)]

        return  fields[i]

    def PrintFields(self):
        '''Prepare selected fields for printout'''

        log = "Selected fields:\n"

        for f in self.frange:
            log += "F %4d  %s\n" % (f, self.Fields(f,self.fstep))

        return log

#----------------------------------------------------------------------------
# Gamess (US) routines
#----------------------------------------------------------------------------

class GAMESS_INPUTS(INPUTS):
    """Gamess US input routines"""

    def __init__(self, data, fstep, frange, grr):
        # template name
        self.pkg = "gamess"
        # template content
        self.tmpl="""\
 $system mwords=10 memddi=0 parall=.t. $end
 $contrl scftyp=rhf runtyp=ffield icharg=0 mult=1 units=angs
         maxit=100 nprint=8 exetyp=run mplevl=0 ispher=-1
         icut=20 itol=30 aimpac=.f. cctyp=none $end
!$eds    mch(1)=0,0 mmul(1)=1,1 mnr(1)=1 ffeds=.t. $end
!$mp2    mp2prp=.t. code=ddi cutoff=1d-20 $end
 $ffcalc offdia=.t. estep=@fstep onefld=@onefld
         efield(1)=@field $end
!$efield evec(1)=@field sym=.f. $end
 $ccinp  iconv=14 maxcc=100 $end
 $trans  cuttrf=1d-15 $end
 $scf    dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
         conv=1d-11 swdiis=0.0001 $end
 $basis  gbasis=sto ngauss=3 $end
@data
"""
        INPUTS.__init__(self, data, fstep, frange, grr)

    def WriteInputs(self):

        # initialize periodic table
        p=Periodic(0)
        
        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)
        
        for i in range(len(xyz)):
            xyz[i]=xyz[i].split()
            xyz[i].insert(1, str( Atomn(xyz[i][0], p) ))
            xyz[i]='%-5s %5s %15s %15s %15s\n' % tuple([xyz[i][j] for j in range(5)])
        
        xyz.insert(0, ' $data\nFinite Field\nc1 0\n')
        xyz.append(' $end')
        xyz=''.join(xyz)
        
        # write input files
        for f in self.frange:

            # prepare onefld inputs
            filename = self.data.replace('.xyz','')+'_%.4f_F%2d_' % (self.fstep, f)
            filename = filename.replace(' ','0')
            ffield = ','.join(self.Fields(f,self.fstep).split())
            finput = self.tmpl.substitute(data=xyz, fstep=self.fstep, field=ffield, onefld='.t.')
            if f==0:
                finput = finput.replace(' $efield','!$efield')

            # write onefld inputs
            open(filename+'.inp','w').write(finput)

            # write kurtz's ffield input
            if f==0:
                finput = self.tmpl.substitute(data=xyz, fstep=self.fstep, field=ffield, onefld='.f.')
                filename = self.data.replace('.xyz','')+'_ffield_%.4f' % self.fstep
                filename = filename.replace(' ','_')
                open(filename+'.inp','w').write(finput)

class GAMESS(PARSER):
    """A GAMESS log parser."""

    def ParseFile(self,filename):
        """Read gamess(us) energies for this system."""

        self.log = open(filename)

        # Gamess conversion factor
        self.au2d=2.541766

        if FindLine(self.log,'$CONTRL OPTIONS') !=-1:

            # Read control options
            line = SkipLines(self.log,2)
            self.scftyp = re.split('\s+|=',line)[2].lower()
            self.runtyp = re.split('\s+|=',line)[4].lower()
            line = SkipLines(self.log,1)
            self.mplevl = int(re.split('\s+|=',line)[3])
            self.cctyp = re.split('\s+|=',line)[9].upper()
            if self.mplevl == 2:
                line = FindLine(self.log,'MP2 CONTROL INFORMATION')
                line = SkipLines(self.log,5)
                self.mp2prp = re.split('\s+|=',line)[4].lower()

            # Try parsing this file
            try:

                if self.runtyp == 'eds':
                    term_ok = self.ParseEds()
                else:
                    term_ok = self.ParseGms()

                if not term_ok:
                    print "Warning: file ", filename, " did not finished properly!"

            except AttributeError:
                print "Attribute error. Check for junk or redundand log files."
                sys.exit(1)
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise

        else:
            print "File: ", filename, " is probably not a Gamess US file or has an unexpected format."

    def ParseGms(self):
        """Read gamess(us) energies for this system."""

        termination_code = False

        # Find field, energies and dipoles
        while 1:
            line = self.log.readline()

            if line.find('EXECUTION OF GAMESS TERMINATED NORMALLY') !=-1:
                termination_code = True
                break
            if line == '':
                break
    
            # Case of FF calculations
            if line.find('APPLIED FIELD') !=-1:
    
                # ... determine best energy ...
                if self.mplevl == 2:
                    EnKey = 'MP2'
                elif self.cctyp != 'NONE':
                    EnKey = self.cctyp
                else:
                    EnKey = 'SCF'
                if EnKey not in self.energies:
                    self.energies[EnKey]={}
    
                line = line.split()

                # ... set field label ...
                field = ( round(float(line[-3]),9),
                          round(float(line[-2]),9),
                          round(float(line[-1]),9) )

                # ... update base field ...
                self.SetBaseField(field)

                if field not in self.energies[EnKey]:
                    self.energies[EnKey][field]= float64(SkipLines(self.log,1).split()[-1])
    
            # Case of single field calculations
            if line.find('ELECTRIC FIELD') !=-1:
                line = line.split()
                # ... set field label ...
                field = ( round(float(line[-3]),9),
                          round(float(line[-2]),9),
                          round(float(line[-1]),9) )
                # ... update base field ...
                self.SetBaseField(field)
                # ... read scf energy ...
                line = FindLine(self.log,'SCF CALCULATION')
                if 'SCF' not in self.energies:
                    self.energies['SCF']={}
                if 'SCF' not in self.dipoles:
                    self.dipoles['SCF']={}
                if field not in self.energies['SCF']:
                    self.energies['SCF'][field]= float64(FindLine(self.log,' FINAL ').split()[4])
                # ... read scf dipole ...
                if field not in self.dipoles['SCF']:
                    line = FindLine(self.log,'ELECTROSTATIC MOMENTS')
                    line = SkipLines(self.log,6)
                    self.dipoles['SCF'][field] = array(line.split()[:3],dtype=float64)/self.au2d
                # ... read mp2 energy ...
                if self.mplevl == 2:
                    if 'MP2' not in self.energies:
                        self.energies['MP2']={}
                    if field not in self.energies['MP2']:
                        self.energies['MP2'][field]= float64(FindLine(self.log,'E(MP2)=').split()[-1])
                if self.mplevl == 2 and self.mp2prp == 't':
                    if 'MP2' not in self.dipoles:
                        self.dipoles['MP2']={}
                    if field not in self.dipoles['MP2']:
                        line = FindLine(self.log,'MP2 PROPERTIES')
                        line = FindLine(self.log,'ELECTROSTATIC MOMENTS')
                        line = SkipLines(self.log,6)
                        self.dipoles['MP2'][field] = array(line.split()[:3],dtype=float64)/self.au2d
                # ... read cc energy ...
                if self.cctyp != 'NONE':
                    if 'MP2' not in self.energies:
                        self.energies['MP2']={}
                    if field not in self.energies['MP2']:
                        self.energies['MP2'][field]= float64(FindLine(self.log,'MBPT(2) ENERGY:').split()[-3])
                    if 'CCSD' not in self.energies:
                        self.energies['CCSD']={}
                    if field not in self.energies['CCSD']:
                        self.energies['CCSD'][field]= float64(FindLine(self.log,'CCSD    ENERGY:').split()[-3])
                if self.cctyp == 'CCSD(T)':
                    if 'CCSD(T)' not in self.energies:
                        self.energies['CCSD(T)']={}
                    if field not in self.energies['CCSD(T)']:
                        self.energies['CCSD(T)'][field]= float64(FindLine(self.log,'CCSD(T) ENERGY:').split()[-3])

        return termination_code

    def ParseEds(self):
        """Read gamess(us) eds energies for this system."""

        termination_code = False

        # Check onefld value
        line = FindLine(self.log,'FFCALC INPUT OPTIONS')
        line = SkipLines(self.log,4)
        onefld = re.split('\s*=\s*|\s*',line)[2].lower()

        if onefld == 't':

            # Case of single field calculations
            line = FindLine(self.log,'APPLIED FIELD')
            if line !=-1:
                l = line.split()
                # ... set field label ...
                field = ( round(float(l[-3]),9),
                          round(float(l[-2]),9),
                          round(float(l[-1]),9) )
                # ... update base field ...
                self.SetBaseField(field)
            else:
                return termination_code

            # Find field, interaction energies and total energies and dipoles
            while 1:
                line = self.log.readline()

                if line == '':
                    break
                if line.find('GAMESS TERMINATED NORMALLY') !=-1:
                    termination_code = True
                    break
            
                if line.find('MONOMER:') !=-1:
                    line = SkipLines(self.log,2)
                    conf_no = re.split('\D+',line)[1]
                    conf_id = ''.join(re.split('\D+',line)[2:-1])

                    # ... read scf energy ...
                    line = FindLine(self.log,'SCF CALCULATION')
                    vale = float64(FindLine(self.log,' FINAL ').split()[4])
                    self.SaveEnergy('SCF',conf_id,conf_no,field,vale)

                    # ... read scf interaction energy terms ...
                    if self.Nmer(conf_id) >= 2:
                        self.ReadScfTotTerms(conf_id,conf_no,field)
                    if self.Nmer(conf_id) == len(conf_id) and self.Nmer(conf_id) >= 3:
                        self.ReadScfMnbTerms('n-body',conf_no,field)

                    # ... read scf dipole ...
                    line = FindLine(self.log,'ELECTROSTATIC MOMENTS')
                    line = SkipLines(self.log,6)
                    vald = array(line.split()[:3],dtype=float64)/self.au2d
                    self.SaveDipole('SCF',conf_id,conf_no,field,vald)

                    # ... read mp2 energy ...
                    if self.mplevl == 2:
                        line = FindLine(self.log,'E(MP2)=')
                        vale = float64(line.split()[-1])
                        self.SaveEnergy('MP2',conf_id,conf_no,field,vale)

                    # ... read mp2 interaction energy terms ...
                    if self.mplevl == 2 and self.Nmer(conf_id) >= 2:
                        self.ReadMp2TotTerms(conf_id,conf_no,field)
                    if ( self.mplevl == 2 and self.Nmer(conf_id) == len(conf_id) 
                         and self.Nmer(conf_id) >= 3 ): 
                        self.ReadMp2MnbTerms('n-body',conf_no,field)

                    if self.mplevl == 2 and self.mp2prp == 't':
                        line = FindLine(self.log,'ELECTROSTATIC MOMENTS')
                        line = SkipLines(self.log,6)
                        vald = array(line.split()[:3],dtype=float64)/self.au2d
                        self.SaveDipole('MP2',conf_id,conf_no,field,vald)

        return termination_code

    def Nmer(self,conf_id):
        """Return n for n-mer having conf_id."""
        n = array(list(conf_id),dtype=int).sum()
        return n

    def ReadScfTotTerms(self,conf_id,conf_no,field):
        """Read interaction energies for this system."""

        line   = FindLine(self.log,'  INTERACTION ENERGY TERMS')
        line   = SkipLines(self.log,4)

        e = {}

        while 1:

            line = self.log.readline()

            if line == '':
                break

            # ... all energies are read ...
            if line.find(20*'-') !=-1:
                # ... save ...
                self.SaveEnergy('DG(HL)',conf_id,conf_no,field,e['DE(HL)']+e['DG(HL,F)'])
                self.SaveEnergy('DG(DEL,HF)',conf_id,conf_no,field,e['DE(DEL,HF)']+e['DG(DEL,F)'])
                self.SaveEnergy('DG(HF)',conf_id,conf_no,field,e['DG(HF,F)'])
                if self.Nmer(conf_id) == 2:
                    self.SaveEnergy('E(EL,10)',conf_id,conf_no,field,e['E(EL,10)'])
                    self.SaveEnergy('G(EX,HL)',conf_id,conf_no,field,e['E(EX,HL)']+e['DG(HL,F)'])

                # ... and exit loop ...
                break

            # ... stuck all energies in e{} ...
            if line.find('(') !=-1:
                label = line.split()[0]
                value = line.split()[1]

                e[label]=float64(value)

    def ReadMp2TotTerms(self,conf_id,conf_no,field):
        """Read interaction energies for this system."""

        line   = FindLine(self.log,'  INTERACTION ENERGY TERMS')
        line   = SkipLines(self.log,4)

        e = {}

        while 1:

            line = self.log.readline()

            if line == '':
                break

            # ... all energies are read ...
            if line.find(20*'-') !=-1:

                # ... save ...
                for k in e.keys():
                    if k != 'DE(MP2)':
                        self.SaveEnergy(k,conf_id,conf_no,field,e[k])

                # ... and exit loop ...
                break

            # ... stuck all energies in e{} ...
            if line.find('(') !=-1:
                label = line.split()[0]
                value = line.split()[1]

                e[label]=float64(value)

    def ReadScfMnbTerms(self,conf_id,conf_no,field):
        """Read SCF many-body partitioning for this system."""

        line   = FindLine(self.log,'MANY BODY INTERACTION ENERGY TERMS')
        line   = SkipLines(self.log,4)

        e = {}

        while 1:

            line = self.log.readline()

            if line == '':
                break

            # ... all energies are read ...
            if line.find(20*'-') !=-1:

                # ... save ...
                for k in e.keys():
                    if k.find('E(EL,10)') !=-1:
                        self.SaveEnergy(k,conf_id,conf_no,field,e[k])
                    if k.find('E(EX,HL)') !=-1:
                        label = 'G(EX,HL),'+k.split(',')[-1]
                        value = e[k]+e['DG(HL,F),'+k.split(',')[-1]]
                        self.SaveEnergy(label,conf_id,conf_no,field,value)
                    if k.find('DE(DEL,HF)') !=-1:
                        label = 'DG(DEL,HF),'+k.split(',')[-1]
                        value = e[k]+e['DG(DEL,F),'+k.split(',')[-1]]
                        self.SaveEnergy(label,conf_id,conf_no,field,value)
                    if k.find('DG(HF)') !=-1:
                        self.SaveEnergy(k,conf_id,conf_no,field,e[k])

                # ... and exit loop ...
                break

            # ... stuck all energies in e{} ...
            if line.find('(') !=-1:

                l = line.split()

                if len(l) == 3:
                    label = l[0]
                    value = l[1]
                if len(l) == 4:
                    label = ''.join(l[:2])
                    value = l[2]

                e[label]=float64(value)

    def ReadMp2MnbTerms(self,conf_id,conf_no,field):
        """Read MP2 many-body partitioning for this system."""

        line   = FindLine(self.log,'MANY BODY INTERACTION ENERGY TERMS')
        line   = SkipLines(self.log,4)

        e = {}

        while 1:

            line = self.log.readline()

            if line == '':
                break

            # ... all energies are read ...
            if line.find(20*'-') !=-1:

                # ... save ...
                for k in e.keys():
                    if k != 'DE(MP2)':
                        self.SaveEnergy(k,conf_id,conf_no,field,e[k])

                # ... and exit loop ...
                break

            # ... stuck all energies in e{} ...
            if line.find('(') !=-1:

                l = line.split()

                if len(l) == 3:
                    label = l[0]
                    value = l[1]
                if len(l) == 4:
                    label = ''.join(l[:2])
                    value = l[2]

                e[label]=float64(value)

#----------------------------------------------------------------------------
# Molcas routines
#----------------------------------------------------------------------------

class MOLCAS_INPUTS(INPUTS):
    """Molcas input routines"""

    def __init__(self, data, fstep, frange, grr):
        # template name
        self.pkg = "molcas"
        # template content
        self.tmpl="""\
 &GATEWAY
Coord
@data
Basis=STO-3G
Title=Finite field calculations
Group=NoSym

 &SEWARD

@field

 &SCF
thre     = 1.0e-15, 1.0e-10, 1.0e-08, 1.0e-07

 &MBPT2

 &RASSCF
outorb   = canonical

 &MOTRA
frozen = 0

 &CCSDT
CCT
iterations
200
accuracy = 1.0e-15

 &RASSCF
nactel   = ? 0 0
inactive = ?
ras2     = ?
thrs     = 1.0e-15, 1.0e-06, 1.0e-06

 &CASPT2

"""
        INPUTS.__init__(self, data, fstep, frange, grr)

    def MakeCoords(self):
        '''read xyz file'''
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[0:int(xyz[0])+2]
            xyz=''.join(xyz).rstrip()
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)

        return xyz

    def WriteInputs(self):

        # read xyz file
        xyz=self.MakeCoords()

        # write input files
        for f in self.frange:
            filename = self.data.replace('.xyz','')+'_%.4f_F%2d_' % (self.fstep, f)
            filename = filename.replace(' ','0')

            # Prepare FFPT
            FFPT=' &FFPT\nDIPO\n'
            F=self.Fields(f,self.fstep).split()
            FFPT+='X ' + F[0] + '\n'
            FFPT+='Y ' + F[1] + '\n'
            FFPT+='Z ' + F[2] + '\n'

            finput = self.tmpl.substitute(data=xyz, field=FFPT)
            open(filename+'.inp','w').write(finput)

class MOLCAS_NATIVE_INPUTS(INPUTS):
    """Molcas input routines"""

    def __init__(self, data, fstep, frange, grr):
        # template content
        self.tmpl="""\
 &GATEWAY
@data
Title=Finite field calculations
Group=NoSym

 &SEWARD

@field

 &SCF
thre     = 1.0e-15, 1.0e-10, 1.0e-08, 1.0e-07

 &MBPT2

 &RASSCF
outorb   = canonical

 &MOTRA
frozen = 0

 &CCSDT
CCT
iterations
200
accuracy = 1.0e-15

 &RASSCF
nactel   = ? 0 0
inactive = ?
ras2     = ?
thrs     = 1.0e-15, 1.0e-06, 1.0e-06

 &CASPT2

"""
        MOLCAS_INPUTS.__init__(self, data, fstep, frange, grr)

    def MakeCoords(self):
        '''read native xyz file'''
        try:
            xyz=open(self.data,'r').readlines()
            xyz=''.join(xyz).rstrip()
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)

        return xyz

class MOLCAS(PARSER):
    """A MOLCAS log parser."""

    def ParseFile(self,filename):
        '''Read molcas energies for this system.'''

        self.log=open(filename)

        # Molcas conversion factor
        au2d=2.5417463077102775

        # Find field, energies and dipoles
        while 1:
            line = self.log.readline()
            if line == '': break

            if line.find('FFPT    DIPO    COMP') !=-1:
                line = line.split()
                # ... set field label ...
                field = ( self.sign*round(float(line[-5]),9),
                          self.sign*round(float(line[-3]),9),
                          self.sign*round(float(line[-1]),9) )
                # ... update base field ...
                self.SetBaseField(field)
                # ... read all energies for this field ...
                while 1:
                    line = self.log.readline()
                    if line == '': break
                    #if line.find('MOLCAS executing module FFPT') !=-1: break
                    if line.find('FFPT    DIPO    COMP') !=-1: break
                    # ... read scf energy and dipole ...
                    if line.find('MOLCAS executing module SCF') !=-1:
                        if ('SCF' not in self.energies):
                            self.energies['SCF']={}
                        if (field not in self.energies['SCF']):
                            self.energies['SCF'][field]= float64(ReFindLine(self.log,'Total .+ energy').split()[-1])
                    if line.find('Dipole Moment') !=-1:
                            if ('SCF' not in self.dipoles):
                                self.dipoles['SCF']={}
                            if (field not in self.dipoles['SCF']):
                                dipl = SkipLines(self.log,2).split()
                                self.dipoles['SCF'][field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    if line.find('1-st cartesian moments') !=-1:
                            if ('SCF' not in self.dipoles):
                                self.dipoles['SCF']={}
                            if (field not in self.dipoles['SCF']):
                                dipl = SkipLines(self.log,6).split()
                                self.dipoles['SCF'][field] = array([dipl[1],dipl[2],dipl[3]],dtype=float64)/au2d
                    # ... read mp2 energy and dipole ...
                    if line.find('MOLCAS executing module MBPT2') !=-1:
                        if ('MP2' not in self.energies):
                            self.energies['MP2']={}
                        if (field not in self.energies['MP2']):
                            self.energies['MP2'][field]= float64(FindLine(self.log,'Total energy').split()[-2])
                    # ... read rasscf energy and dipole ...
                    if line.find('MOLCAS executing module RASSCF') !=-1:
                        if ('RASSCF' not in self.energies):
                            self.energies['RASSCF']={}
                        if (field not in self.energies['RASSCF']):
                            line = FindLine(self.log,'Final state energy')
                            line = SkipLines(self.log,3)
                            self.energies['RASSCF'][field] = float64(self.re_number.findall(line)[-1])
                        if ('RASSCF' not in self.dipoles):
                            self.dipoles['RASSCF']={}
                        if (field not in self.dipoles['RASSCF']):
                            line = FindLine(self.log,'Dipole Moment')
                            dipl = SkipLines(self.log,2).split()
                            self.dipoles['RASSCF'][field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    # ... read caspt2 energy and dipole ...
                    if line.find('MOLCAS executing module CASPT2') !=-1:
                        if ('CASPT2' not in self.energies):
                            self.energies['CASPT2']={}
                        if (field not in self.energies['CASPT2']):
                            line = FindLine(self.log,'FINAL CASPT2 RESULT')
                            self.energies['CASPT2'][field]= float64(SkipLines(self.log,5).split()[-1])
                        if ('CASPT2' not in self.dipoles):
                            self.dipoles['CASPT2']={}
                        if (field not in self.dipoles['CASPT2']):
                            line = FindLine(self.log,'Dipole Moment')
                            dipl = SkipLines(self.log,2).split()
                            self.dipoles['CASPT2'][field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    # ... read ccsdt energy ...
                    if line.find('MOLCAS executing module CCSD(T)') !=-1:
                        if ('CCSD' not in self.energies):
                            self.energies['CCSD']={}
                        if (field not in self.energies['CCSD']):
                            line = FindLine(self.log,'Total energy (diff)')
                            self.energies['CCSD'][field]= float64(line.split()[-2])
                        while 1:
                            line = self.log.readline()
                            if line == '':
                                break
                            if line.find('Stop Module: ccsdt') !=-1:
                                break
                            if line.find('CCSD + T3') !=-1:
                                if ('CCSD(T)' not in self.energies):
                                    self.energies['CCSD(T)']={}
                                self.energies['CCSD(T)'][field]= float64(line.split()[-1])

#----------------------------------------------------------------------------
# Gaussian routines
#----------------------------------------------------------------------------

class GAUSSIAN_INPUTS(INPUTS):
    """Gaussian 09 input routines"""

    def __init__(self, data, fstep, frange, grr):
        # template name
        self.pkg = "gaussian"
        # template content
        self.tmpl="""\
%chk=@chk
%mem=512mb
%nproc=2
#p hf/sto-3g scf(conver=11,novaracc,xqc) nosymm field=read
maxdisk=20gb iop(5/7=512) iop(3/27=20) iop(11/27=20)
integral(grid=ultrafine)

gaussian ffield

0 1
@data
@field
"""
        INPUTS.__init__(self, data, fstep, frange, grr)

    def WriteInputs(self):

        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            xyz=''.join(xyz)
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)

        fieldtail=''
        #for i in range(10):
        #    fieldtail+='%7.4f\n' % 0.0

        # write input files
        for f in self.frange:
            filename = self.data.replace('.xyz','')+'_%.4f_F%2d_' % (self.fstep, f)
            filename = filename.replace(' ','0')
            ffield = self.Fields(f,self.fstep)+fieldtail
            finput = self.tmpl.substitute(data=xyz, field=ffield, chk=filename+'.chk')
            open(filename+'.inp','w').write(finput)

class GAUSSIAN(PARSER):
    """A gaussian 09 log parser."""

    def ParseFiles(self):
        """read logfiles in logpath"""

        self.fstep = 0

        if 'g09root' not in os.environ:
            #os.putenv('g09root','/opt/gaussian')
            os.environ['g09root']='/opt/gaussian'

        logfiles=[]
        for f in os.listdir(self.logpath):
            if re.search("\.chk", f):
                logfiles.append(self.logpath+'/'+f)
        logfiles.sort()

        for chk in logfiles:
            # prepare formated checkpoints
            fchk = chk.replace('chk','fchk')
            if not os.path.exists(fchk):
                print "Using %s to format checkpoints" % (os.environ['g09root']+'/g09/formchk')
                Run(os.environ['g09root']+'/g09/formchk', chk + ' ' + fchk)

            # parse current logfile
            self.ParseFile(fchk)

        if self.density == 1:

            for chk in logfiles:

                # ... read applied external field ...
                fchk = chk.replace('chk','fchk')
                if os.path.exists(fchk):
                    log=open(fchk)
                    line = FindLine(log,'External E-field')
                    line = SkipLines(log,1).split()
                    field = ( self.sign*round(float(line[1]),9), 
                              self.sign*round(float(line[2]),9),
                              self.sign*round(float(line[3]),9) )

                # parse current cube file
                cubefile = chk.replace('chk','cube')
                if os.path.exists(cubefile):
                    t = CUBE(cubefile)
                    self.densities[field]=t.data
                    #t.writeCubeFile(chk.replace('chk','cube2'))
                if field == (0,0,0):
                    t.data = None
                    self.cubedata=t

    def ParseFile(self,filename):
        """Read g09 energies and dipoles for this system."""

        self.log=open(filename)

        # Read control options
        puhf=0
        mp2=0
        pmp2=0
        mp2dip=0
        mp3=0
        pmp3=0
        mp4=0
        ccsd=0
        ccsdt=0

        # available energies
        line=SkipLines(self.log,2)
        if line.find('MP2') !=-1:
            mp2=1
        if line.find('UMP2') !=-1:
            puhf=1
            pmp2=1
        if line.find('CCSD') !=-1:
            mp2=1
            mp2dip=0
            mp3=1
            mp4=1
            ccsd=1
        if line.find('CCSD(T)') !=-1:
            ccsdt=1
        if line.find('UCCSD') !=-1:
            puhf=1
            pmp2=1
            pmp3=1

        # available dipoles
        if mp2==1:
            if FindLine(self.log,'Total MP2 Density') !=-1:
                mp2dip=1
            else:
                mp2dip=0

        # setup dictionaries
        if 'SCF' not in self.energies:
            self.energies['SCF']={}
        if puhf==1 and 'PUHF' not in self.energies:
            self.energies['PUHF']={}
        if mp2==1 and 'MP2' not in self.energies:
            self.energies['MP2']={}
        if pmp2==1 and 'PMP2' not in self.energies:
            self.energies['PMP2']={}
        if mp3==1 and 'MP3' not in self.energies:
            self.energies['MP3']={}
        if pmp3==1 and 'PMP3' not in self.energies:
            self.energies['PMP3']={}
        if mp4==1 and 'MP4(SDQ)' not in self.energies:
            self.energies['MP4(SDQ)']={}
        if ccsd==1 and 'CCSD' not in self.energies:
            self.energies['CCSD']={}
        if ccsdt==1 and 'CCSD(T)' not in self.energies:
            self.energies['CCSD(T)']={}

        self.log.seek(0)

        # ... read SCF energy ...
        line = FindLine(self.log,'SCF Energy')
        ESCF = float64(line.split()[-1])

        # ... read MP2 energy ...
        if mp2==1:
            line = FindLine(self.log,'MP2 Energy')
            EMP2 = float64(line.split()[-1])

        # ... read MP3 energy ...
        if mp3==1:
            line = FindLine(self.log,'MP3 Energy')
            EMP3 = float64(line.split()[-1])

        # ... read MP4(SDQ) energy ...
        if mp4==1:
            line = FindLine(self.log,'MP4SDQ Energy')
            EMP4 = float64(line.split()[-1])

        # ... read CCSD energy ...
        if ccsd==1:
            line = FindLine(self.log,'Cluster Energy  ')
            ECCSD = float64(line.split()[-1])

        # ... read CCSD(T) energy ...
        if ccsdt==1:
            line = FindLine(self.log,'Cluster Energy with triples')
            ECCSDT = float64(line.split()[-1])

        # ... read PUHF energy ...
        if puhf==1:
            line = FindLine(self.log,'PUHF Energy')
            EPUHF = float64(line.split()[-1])

        # ... read PMP2 energy ...
        if pmp2==1:
            line = FindLine(self.log,'PMP2-0 Energy')
            EPMP2 = float64(line.split()[-1])

        # ... read PMP3 energy ...
        if pmp3==1:
            line = FindLine(self.log,'PMP3-0 Energy')
            EPMP3 = float64(line.split()[-1])

        # ... read applied external field ...
        line = FindLine(self.log,'External E-field')
        line = SkipLines(self.log,1).split()
        field = ( self.sign*round(float(line[1]),9), 
                  self.sign*round(float(line[2]),9),
                  self.sign*round(float(line[3]),9) )

        # ... update base field ...
        self.SetBaseField(field)

        if field not in self.energies['SCF']:
            self.energies['SCF'][field]=ESCF
        if puhf==1 and field not in self.energies['PUHF']:
            self.energies['PUHF'][field]=EPUHF
        if mp2==1 and field not in self.energies['MP2']:
            self.energies['MP2'][field]=EMP2
        if pmp2==1 and field not in self.energies['PMP2']:
            self.energies['PMP2'][field]=EPMP2
        if mp3==1 and field not in self.energies['MP3']:
            self.energies['MP3'][field]=EMP3
        if pmp3==1 and field not in self.energies['PMP3']:
            self.energies['PMP3'][field]=EPMP3
        if mp4==1 and field not in self.energies['MP4(SDQ)']:
            self.energies['MP4(SDQ)'][field]=EMP4
        if ccsd==1 and field not in self.energies['CCSD']:
            self.energies['CCSD'][field]=ECCSD
        if ccsdt==1 and field not in self.energies['CCSD(T)']:
            self.energies['CCSD(T)'][field]=ECCSDT

        # ... read dipoles ...
        if FindLine(self.log,'Dipole Moment') != -1:

            if mp2dip == 0:
                if 'SCF' not in self.dipoles:
                    self.dipoles['SCF']={}
            else:
                if 'MP2' not in self.dipoles:
                    self.dipoles['MP2']={}

            line = SkipLines(self.log,1)

            D = array(line.split(),dtype=float64)

            if mp2dip == 0 and field not in self.dipoles['SCF']:
                self.dipoles['SCF'][field]=D

            if mp2dip == 1 and field not in self.dipoles['MP2']:
                self.dipoles['MP2'][field]=D

#----------------------------------------------------------------------------
# Utilities
#----------------------------------------------------------------------------

def Run(program, args):
    '''Run shell command'''
    os.system(program + ' ' + args)

def SkipLines(open_file,n):
    '''Read n lines from file f.'''

    for i in range(n):
        line = open_file.readline()
        if line == '' :
            break

    return line

def FindLine(open_file,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = open_file.readline()
        if line.find(pattern) !=-1 :
            break
        if line == '' :
            line = -1
            break

    return line

def ReFindLine(open_file,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = open_file.readline()
        if re.search(pattern,line) !=None : break
        if line == '' :
            linie = -1
            break

    return line

def PunchEn(out,i,field,energies):
    '''Punch energies to out file'''
    row='%26.15f ' % energies[F]
    row+='%7.4f' % field[i]
    row+='\t\tF=%s' % str(field)
    out.write( row+'\n' )

def FF25Input(energies, fstep):
    '''Write input for RZ's ff25.f'''
    out=open('ff25.inp','w')
    out.write('%f\n' % fstep)
    for i in range(25):
        field=fields(i,fstep).split()
        for j in range(len(field)):
            field[j]=float64(field[j])
        out.write('%25.16e _%d_\n' % (energies[tuple(field)], i+1))

#----------------------------------------------------------------------------
# Generalized Romberg analysis
#----------------------------------------------------------------------------

class ROMBERG:
    """Generalized romberg analysis"""

    def __init__(self, energies, sign, label):
        """Expects energies cast in a following dictionary:
           energies{ (fx, fy, fz): value, ... }
        """
        self.energies = energies
        self.sign = sign # -1 for energies and +1 for dipoles
        self.label = label
        self.fields = array(energies.keys(),dtype=float64)
        self.max_field = self.fields.max(axis=0)
        self.rr_log = {}
        self.rr = {}
        self.enforce = False

        if self.max_field.max() == self.max_field.sum() > 0 or self.enforce:

            self.axis = self.max_field.argmax()
            self.romberg = abs(self.fields.sum(axis=0)[self.axis]) <= 1.0e-10

            if self.romberg:

                try:
                    # setup fields
                    self.fields = self.fields[:,self.axis]
                    self.fields.sort(axis=0)

                    # n: number steps; h: step size; a: quotient;
                    self.n = (self.fields.size-1)/2
                    self.h = self.fields[self.n+1:][0]
                    self.a = self.fields[0] / self.fields[1]

                    # recast energies
                    energies = {}
                    for f,e in self.energies.items():
                        energies[ f[self.axis] ] = float64(e)

                    self.energies = energies

                    # perform the analysis
                    for o in [1, 2, 3, 4]:
                        self.GR(o)
                        print self.rr_log[o]

                except KeyError:
                    pass

    def GR(self,order):
        """Generalized Romberg estimates"""

        # Set number of GR iterations and output format
        if order == 1:
            n = self.n
            line = '%10.4f'
            sufix='st'
        if order == 2:
            n = self.n
            line = '%10.3f'
            sufix='nd'
        if order == 3:
            n = self.n - 1
            line = '%10.2f'
            sufix='rd'
        if order == 4:
            n = self.n - 1
            line = '%10.1f'
            sufix='th'

        # Romberg iterations forming P[k,p]
        P = []
        dP = {}

        for p in range(n):

            if p == 0:
                P.append( [ self.D(order,k) for k in range(n) ] )
            else:
                P.append( [ ( self.a**(2*p)*P[p-1][k] - P[p-1][k+1] ) / 
                            ( self.a**(2*p)-1 ) for k in range(n-p) ] )

        # Format results for final printout and storage
        log  = 'GR scheme %1d%s derivative (a=%3.1f) %s\n' % (order, sufix, self.a, self.label)
        log += (9+10*n)*'-'
        log += '\n%9s' % ( 'field / k' )
        log += '%10d'*n % tuple( [k for k in range(n) ] )
        log += '\n' + (9+10*n)*'-'

        for i in range(len(P)):
            log += '\n%9.4f' % ( self.h*self.a**i )

            for j in range(len(P[i])):
                log += line % P[j][i]
                if j > 0:
                    dP[(j,i)] = ( 2*P[j][i] - P[j-1][i] - P[j-1][i+1] )/2

        try:
            r_i = abs( array(dP.values()) ).argmin()
            r_k = dP.keys()[r_i] 

            log += '\n' + (9+10*n)*'-'
            log += '\n P[p=%d,k=%d] = ' % ( r_k[1]+1, r_k[0] )
            self.rr[order] = P[r_k[0]][r_k[1]]
            log += line % self.rr[order]
            log += ' +/- ' + line % dP.values()[r_i]
            log += '\n' + (9+10*n)*'-'
        except ValueError:
            self.rr[order] = NaN

        log += '\n'

        self.rr_log[order] = log

        # Wrap up and exit
        return 0

        # Just for testing
        log  = 'GR scheme %1d derivative (a=%3.1f)\n' % (order, self.a)
        log += (9+10*n)*'-'
        log += '\n%9s' % ( 'field / k' )
        log += '%10d'*n % tuple( [k for k in range(n) ] )
        log += '\n' + (9+10*n)*'-'

        for i in range(len(P)):
            log += '\n%9.4f' % ( self.h*self.a**i )

            for j in range(len(P[i])):
                if j > 0:
                    dP[(j,i)] = ( 2*P[j][i] - P[j-1][i] - P[j-1][i+1] )/2
                    log += line % dP[(j,i)]
                else:
                    log += line % NaN

        log += '\n'
        print log

    # Zeroth-order central-difference derivatives
    def D(self,order,k):
        """Finite-difference derivatives"""

        # ... 1st derivative ...
        if order == 1:
            D = ( self.energies[ self.h*self.a**k] - 
                  self.energies[-self.h*self.a**k] 
                ) / ( 2*self.h*self.a**k )

        # ... 2nd derivative ...
        if order == 2:
            D = ( self.energies[-self.h*self.a**k] +
                  self.energies[ self.h*self.a**k] - 
                  2*self.energies[0.0]
                ) / (self.h*self.a**k)**2

        # ... 3rd derivative ...
        if order == 3:
            D = 3*( self.a*self.energies[-self.h*self.a**k] -
                    self.a*self.energies[ self.h*self.a**k] +
                    self.energies[ self.h*self.a**(k+1)] -
                    self.energies[-self.h*self.a**(k+1)]
                  ) / ( self.a*(self.a**2-1)*(self.h*self.a**k)**3 )

        # ... 4th derivative ...
        if order == 4:
            D = 12*( ( self.energies[-self.h*self.a**(k+1)] +
                       self.energies[ self.h*self.a**(k+1)] ) -
                     ( (self.a**2)*self.energies[-self.h*self.a**k] +
                       (self.a**2)*self.energies[ self.h*self.a**k] ) +
                       2*(self.a**2-1)*self.energies[0.0]
                   ) / ( (self.a**2)*(self.a**2-1)*(self.h*self.a**k)**4 )

        # ... return ...
        return self.sign*D

#----------------------------------------------------------------------------
# Property calculation routines
#----------------------------------------------------------------------------

class FDIFF:

    def __init__(self, base_property, fstep, units, method):
        self.base_property = base_property
        self.fstep = fstep
        self.units = units
        self.method = method

        self.properties={}

        self.Calculate()
        self.Log()
        print self.log

    def Calculate(self):
        '''Calculate properties using finite field method (dipole expansion).

           self.properties dictionary stores properties as tensors:
           {M, A, B, G} and isotropic properties {Mv, Av, Bx, By, Bz, Bm, Gv}
           
           individual tensor elements are also stored as:
           self.label   x,   xx, xy, xz,   xxx, yxx, zxx,   xxxx, xxyy, xxzz,
                        y,   yx, yy, yz,   xyy, yyy, zyy,   yyxx, yyyy, yyzz,
                        z,   zx, zy, zz,   xzz, yzz, zzz,   zzxx, zzyy, zzzz,
                       Mv,           Av,    Bx,  By,  Bz,                 Gv, 
                                                      Bm,
        '''

        # set base property and field step
        P = self.base_property
        f = self.fstep  

        # calculate vector/tensor elements
        self.x = x = self.Mu(0, P, f)
        self.y = y = self.Mu(1, P, f)
        self.z = z = self.Mu(2, P, f)
        # pack tensor elements into an array
        self.M = M = array([x,y,z],dtype=float64)
        # vector norm
        self.Mv = Mv = sqrt(dot(M,M))

        self.properties['M'] = M
        self.properties['Mv'] = Mv

        # polarizability
        self.xx = xx = self.Alpha(0,0, P, f)
        self.xy = xy = self.Alpha(0,1, P, f)
        self.xz = xz = self.Alpha(0,2, P, f)
        self.yx = yx = self.Alpha(1,0, P, f)
        self.yy = yy = self.Alpha(1,1, P, f)
        self.yz = yz = self.Alpha(1,2, P, f)
        self.zx = zx = self.Alpha(2,0, P, f)
        self.zy = zy = self.Alpha(2,1, P, f)
        self.zz = zz = self.Alpha(2,2, P, f)
        # pack tensor elements into an array
        self.A = A = array([xx,xy,xz,
                            yx,yy,yz,
                            zx,zy,zz],dtype=float64).reshape(3,3)
        # isotropic average
        self.Av = Av = trace(A)/3.0

        self.properties['A'] = A
        self.properties['Av'] = Av

        # first hyperpolarizability
        self.xxx = xxx = self.Beta(0,0, P, f)
        self.xyy = xyy = self.Beta(0,1, P, f)
        self.xzz = xzz = self.Beta(0,2, P, f)
        self.yxx = yxx = self.Beta(1,0, P, f)
        self.yyy = yyy = self.Beta(1,1, P, f)
        self.yzz = yzz = self.Beta(1,2, P, f)
        self.zxx = zxx = self.Beta(2,0, P, f)
        self.zyy = zyy = self.Beta(2,1, P, f)
        self.zzz = zzz = self.Beta(2,2, P, f)
        # pack tensor elements into an array
        B = array([xxx,yxx,zxx,
                   xyy,yyy,zyy,
                   xzz,yzz,zzz],dtype=float64).reshape(3,3)
        # vector component
        self.Bx = Bx = (3.0/5.0)*B.sum(axis=0)[0]
        self.By = By = (3.0/5.0)*B.sum(axis=0)[1]
        self.Bz = Bz = (3.0/5.0)*B.sum(axis=0)[2]
        # projection to the dipole moment vector
        self.Bm = Bm = (3.0/5.0)*dot(M,B.sum(axis=0))/Mv
        self.Bv = Bv = (1.0/5.0)*sqrt(dot(B.sum(axis=0),B.sum(axis=0)))

        self.properties['B'] = B
        self.properties['Bx'] = Bx
        self.properties['By'] = By
        self.properties['Bz'] = Bz
        self.properties['Bm'] = Bm
        self.properties['Bv'] = Bv

        # second hyperpolarizability
        self.xxxx = xxxx = self.Gamma(0,0, P, f)
        self.xxyy = xxyy = self.Gamma(0,1, P, f)
        self.xxzz = xxzz = self.Gamma(0,2, P, f)
        self.yyxx = yyxx = self.Gamma(1,0, P, f)
        self.yyyy = yyyy = self.Gamma(1,1, P, f)
        self.yyzz = yyzz = self.Gamma(1,2, P, f)
        self.zzxx = zzxx = self.Gamma(2,0, P, f)
        self.zzyy = zzyy = self.Gamma(2,1, P, f)
        self.zzzz = zzzz = self.Gamma(2,2, P, f)
        # pack tensor elements into an array
        self.G = G = array([xxxx,xxyy,xxzz,
                            yyxx,yyyy,yyzz,
                            zzxx,zzyy,zzzz],dtype=float64).reshape(3,3)
        # isotropic average
        self.Gv = Gv = (G.trace()+2.0*(G[0,1]+G[0,2]+G[1,2]))/5.0

        self.properties['G'] = G
        self.properties['Gv'] = Gv

    def Log(self):
        '''Prepare log'''

        # Prepare data for printout
        self.log  = 'Static electric dipole properties (finite field %s, f=%.4f) \n\n' % (self.method, self.fstep)

        self.log += 'Dipole moment [%s]\n\n' % self.units['u']
        line = 4*self.units['m']['f'] + '\n\n'
        self.log += '%15s %15s %15s %15s\n' % ( 'x'.rjust(15), 'y'.rjust(15), 'z'.rjust(15), '<mu>'.rjust(15) )
        self.log += line % ( self.x, self.y, self.z, self.Mv )

        self.log += 'Polarizability [%s]\n\n' % self.units['u']
        line = 4*self.units['a']['f'] + '\n\n'
        self.log += '%15s %15s %15s %15s\n' % ( 'xx'.rjust(15), 'xy'.rjust(15), 'xz'.rjust(15), '<alpha>'.rjust(15) )
        self.log += line % ( self.xx, self.xy, self.xz, self.Av )
        line = 3*self.units['a']['f'] + '\n\n'
        self.log += '%15s %15s %15s\n' % ( 'yx'.rjust(15), 'yy'.rjust(15), 'yz'.rjust(15) )
        self.log += line % ( self.yx, self.yy, self.yz )
        self.log += '%15s %15s %15s\n' % ( 'zx'.rjust(15), 'zy'.rjust(15), 'zz'.rjust(15) )
        self.log += line % ( self.zx, self.zy, self.zz )

        self.log += 'First Hyperpolarizability [%s]\n\n' % self.units['u']
        line = 4*self.units['b']['f'] + '\n\n'
        self.log += '%15s %15s %15s %15s\n' % ( 'xxx'.rjust(15), 'yxx'.rjust(15), 'zxx'.rjust(15), '<beta_mu>'.rjust(15) )
        self.log += line % ( self.xxx, self.yxx, self.zxx, self.Bm )
        line = 3*self.units['b']['f'] + '\n\n'
        self.log += '%15s %15s %15s\n' % ( 'xyy'.rjust(15), 'yyy'.rjust(15), 'zyy'.rjust(15) )
        self.log += line % ( self.xyy, self.yyy, self.zyy )
        self.log += '%15s %15s %15s\n' % ( 'xzz'.rjust(15), 'yzz'.rjust(15), 'zzz'.rjust(15) )
        self.log += line % ( self.xzz, self.yzz, self.zzz )
        line = 4*self.units['b']['f'] + '\n\n'
        self.log += '%15s %15s %15s %15s\n' % ( '<beta_x>'.rjust(15), '<beta_y>'.rjust(15), '<beta_z>'.rjust(15), '<beta>'.rjust(15) )
        self.log += line % ( self.Bx, self.By, self.Bz, self.Bv )

        self.log += 'Second Hyperpolarizability [%s]\n\n' % self.units['u']
        line = 4*self.units['g']['f'] + '\n\n'
        self.log += '%15s %15s %15s %15s\n' % ( 'xxxx'.rjust(15), 'xxyy'.rjust(15), 'xxzz'.rjust(15), '<gamma>'.rjust(15) )
        self.log += line % ( self.xxxx, self.xxyy, self.xxzz, self.Gv )
        line = 3*self.units['g']['f'] + '\n\n'
        self.log += '%15s %15s %15s\n' % ( 'yyxx'.rjust(15), 'yyyy'.rjust(15), 'yyzz'.rjust(15) )
        self.log += line % ( self.yyxx, self.yyyy, self.yyzz )
        self.log += '%15s %15s %15s\n' % ( 'zzxx'.rjust(15), 'zzyy'.rjust(15), 'zzzz'.rjust(15) )
        self.log += line % ( self.zzxx, self.zzyy, self.zzzz )

    def Fi(self,i,bi):
        '''Select diagonal field.'''
        Fi = [0,0,0]
        Fi[i] = bi
        return tuple(Fi)

    def Fij(self,i,j,bi,bj):
        '''Select off-diagonal field.'''
        Fij = [0,0,0]
        Fij[i] = bi
        Fij[j] = bj
        return tuple(Fij)

    def Mu(self,i):
        pass

    def Alpha(self,i,j):
        pass

    def Beta(self,i,j):
        pass

    def Gamma(self,i,j):
        pass

    def D_i(self,i,E,f,np):
        '''Components of dipole moment (energy expansion).'''
        if np == 7:
            try:
                D_i= ( E[Fi(i,3*f)] - E[Fi(i,-3*f)] -
                       9.0*(E[Fi(i,2*f)] - E[Fi(i,-2*f)] +
                       5.0*(E[Fi(i,-f)] - E[Fi(i,f)]) ) ) / (60.0*abs(f))
            except KeyError:
                D_i=NaN

        return D_i

    def D_ij(self,i,j,E,f,np):
        '''Components of polarizability (energy expansion).'''
        if np == 7:
            try:
                if i==j:
                    D_ij=  ( 2.0*( E[Fi(i,3*f)] + E[Fi(i,-3*f)] ) -
                             27.0*( E[Fi(i,2*f)] + E[Fi(i,-2*f)] -
                                    10.0*( E[Fi(i,f)] + E[Fi(i,-f)] )) -
                             490.0*E[Fi(0,0)] ) / (360.0*abs(f*f))
                if i!=j:
                    D_ij=  ( 12.0*( E[Fij(i,j,3*f,f)] + E[Fij(i,j,f,3*f)] -
                                    E[Fij(i,j,3*f,0)] - E[Fij(i,j,0,3*f)] +
                                    E[Fij(i,j,-2*f,0)] + E[Fij(i,j,0,-2*f)] ) +
                             480.0*E[Fij(i,j,f,f)] +
                             330.0*( E[Fi(0,0)] - E[Fij(i,j,0,f)] - E[Fij(i,j,f,0)] ) +
                             160.0*E[Fij(i,j,-f,-f)] -
                             105.0*( E[Fij(i,j,f,2*f)] + E[Fij(i,j,2*f,f)] ) -
                             15.0*( E[Fij(i,j,-2*f,-f)] + E[Fij(i,j,-f,-2*f)] ) +
                             90.0*( E[Fij(i,j,2*f,0)] + E[Fij(i,j,0,2*f)] -
                                    E[Fij(i,j,-f,0)] - E[Fij(i,j,0,-f)] ) -
                             60.0*( E[Fij(i,j,f,-f)] + E[Fij(i,j,-f,f)] ) +
                              5.0*( E[Fij(i,j,2*f,-f)] + E[Fij(i,j,-f,2*f)] ) +
                              3.0*( E[Fij(i,j,-2*f,f)] + E[Fij(i,j,f,-2*f)] ) +
                              10.0*E[Fij(i,j,2*f,2*f)] ) / (360.0*abs(f*f))
            except KeyError:
                D_ij=NaN

        return D_ij

    def D_ijj(self,i,j,E,f,np):
        '''Components of hyperpolarizability (energy expansion).'''
        if np == 7:
            try:
                # Beta_iii
                if i==j:
                    D_ijj=  ( E[Fi(i,-3*f)] - E[Fi(i,3*f)] +
                              8.0*( E[Fi(i,2*f)] - E[Fi(i,-2*f)] ) -
                              13.0*( E[Fi(i,f)] - E[Fi(i,-f)] ) ) / (48.0*abs(f*f*f))
                # Beta_ijj
                if i!=j:
                    D_ijj=  ( E[Fij(i,j,-f,2*f)] + E[Fij(i,j,-f,-2*f)] -
                              E[Fij(i,j,f,2*f)] - E[Fij(i,j,f,-2*f)] +
                              2.0*( E[Fij(i,j,-2*f,f)] + E[Fij(i,j,-2*f,-f)] -
                                    E[Fij(i,j,2*f,f)] - E[Fij(i,j,2*f,-f)] ) +
                              4.0*( E[Fij(i,j,2*f,0)] - E[Fij(i,j,-2*f,0)] ) +
                              20.0*( E[Fij(i,j,f,f)] + E[Fij(i,j,f,-f)] -
                                     E[Fij(i,j,-f,f)] - E[Fij(i,j,-f,-f)] ) -
                              38.0*( E[Fij(i,j,f,0)] - E[Fij(i,j,-f,0)] ) ) / (144.0*abs(f*f*f))
            except KeyError:
                D_ijj=NaN

        return D_ijj

    def D_iijj(self,i,j,E,f,np):
        '''Components of second hyperpolarizability (energy expansion).'''
        if np == 7:
            try:
                # Gamma_iiii
                if i==j:
                    D_iijj=  ( 56.0*E[Fi(0,0)] - E[Fi(i,3*f)] - E[Fi(i,-3*f)] +
                               12.0*( E[Fi(i,2*f)] + E[Fi(i,-2*f)] ) -
                               39.0*( E[Fi(i,f)] + E[Fi(i,-f)] ) ) / (144.0*abs(f*f*f*f))
                # Gamma_iijj
                if i!=j:
                    D_iijj=  ( 72.0*E[Fi(0,0)] - E[Fij(i,j,-2*f,-f)] - E[Fij(i,j,-f,-2*f)] -
                               E[Fij(i,j,-f,2*f)] - E[Fij(i,j,f,-2*f)] -
                               E[Fij(i,j,f,2*f)] - E[Fij(i,j,2*f,-f)] -
                               E[Fij(i,j,2*f,f)] - E[Fij(i,j,-2*f,f)] +
                               2.0*( E[Fij(i,j,0,2*f)] + E[Fij(i,j,2*f,0)] +
                                     E[Fij(i,j,-2*f,0)] + E[Fij(i,j,0,-2*f)] ) +
                               20.0*( E[Fij(i,j,f,-f)] + E[Fij(i,j,f,f)] +
                                      E[Fij(i,j,-f,-f)] + E[Fij(i,j,-f,f)] ) -
                               38.0*( E[Fij(i,j,-f,0)] + E[Fij(i,j,0,-f)] +
                                      E[Fij(i,j,0,f)] + E[Fij(i,j,f,0)] ) ) / (288.0*abs(f*f*f*f))
            except KeyError:
                D_iijj=NaN

        return D_iijj


class FDIFF_E_V(FDIFF):

    def __init__(self, E, fstep, units, method):
        '''Initialize for energy expansion'''
        self.base_name = 'energy'

        FDIFF.__init__(self, E, fstep, units, method)

    def Mu(self,i,E,f):
        '''Components of dipole moment (energy expansion).'''
        try:
            mu_i=( (-2.0/3.0)*( E[Fi(i,f)] - E[Fi(i,-f)] ) + 
                   ( E[Fi(i,2*f)] - E[Fi(i,-2*f)] )/12.0 )/abs(f)
        except KeyError:
            mu_i=NaN

        return mu_i

    def Alpha(self,i,j,E,f):
        '''Components of polarizability (energy expansion).'''
        try:
            if i==j:
                alpha_ij=( (5.0/2.0)*E[Fi(0,0)] - (4.0/3.0)*( E[Fi(i,f)] + E[Fi(i,-f)] ) + 
                           ( E[Fi(i,2*f)] + E[Fi(i,-2*f)] )/12.0 )/abs(f*f)
            if i!=j:
                alpha_ij=( ( E[Fij(i,j,2*f,2*f)] - E[Fij(i,j,2*f,-2*f)] - 
                             E[Fij(i,j,-2*f,2*f)] + E[Fij(i,j,-2*f,-2*f)] )/48.0 -
                           ( E[Fij(i,j,f,f)] - E[Fij(i,j,f,-f)] -
                             E[Fij(i,j,-f,f)] + E[Fij(i,j,-f,-f)] )/3.0 )/abs(f*f)
        except KeyError:
            alpha_ij=NaN

        return alpha_ij

    def Beta(self,i,j,E,f):
        '''Components of hyperpolarizability (energy expansion).'''
        try:
            # Beta_iii
            if i==j:
                beta_ijj=( ( E[Fi(i,f)] - E[Fi(i,-f)] ) -
                           ( E[Fi(i,2*f)] - E[Fi(i,-2*f)] )/2.0 )/abs(f*f*f)
            # Beta_ijj
            if i!=j:
                beta_ijj=( ( E[Fij(i,j,-f,-f)] - E[Fij(i,j,f,f)] +
                             E[Fij(i,j,-f,f)] - E[Fij(i,j,f,-f)] )/2.0 +
                             E[Fi(i,f)] - E[Fi(i,-f)] )/abs(f*f*f)
        except KeyError:
            beta_ijj=NaN

        return beta_ijj

    def Gamma(self,i,j,E,f):
        '''Components of second hyperpolarizability (energy expansion).'''
        try:
            # Gamma_iiii
            if i==j:
                gamma_iijj=( 4.0*( E[Fi(i,f)] + E[Fi(i,-f)] ) - 6.0*E[Fi(0,0)] -
                             ( E[Fi(i,2*f)] + E[Fi(i,-2*f)] ) )/abs(f*f*f*f)
            # Gamma_iijj
            if i!=j:
                gamma_iijj=( -4.0*E[Fi(0,0)] - ( E[Fij(i,j,f,f)] + E[Fij(i,j,-f,-f)] +
                             E[Fij(i,j,f,-f)] + E[Fij(i,j,-f,f)] ) +
                             2.0*( E[Fi(i,f)] + E[Fi(i,-f)] + E[Fi(j,f)] + E[Fi(j,-f)] ) )/abs(f*f*f*f)
        except KeyError:
            gamma_iijj=NaN

        return gamma_iijj

class FDIFF_E_VII(FDIFF):

    def __init__(self, E, fstep, units, method):
        '''Initialize for energy expansion'''
        self.base_name = 'energy'

        FDIFF.__init__(self, E, fstep, units, method)

    def Mu(self,i,E,f):
        '''Components of dipole moment (energy expansion).'''
        mu_i= -1.0*self.D_i(i,E,f,7)
        return mu_i

    def Alpha(self,i,j,E,f):
        '''Components of polarizability (energy expansion).'''
        alpha_ij= -2.0*self.D_ij(i,j,E,f,7)
        return alpha_ij

    def Beta(self,i,j,E,f):
        '''Components of hyperpolarizability (energy expansion).'''
        beta_ijj = -6.0*self.D_ijj(i,j,E,f,7)
        return beta_ijj

    def Gamma(self,i,j,E,f):
        '''Components of second hyperpolarizability (energy expansion).'''
        gamma_iijj = -24.0*self.D_iijj(i,j,E,f,7)
        return gamma_iijj

class FDIFF_DD():

    def __init__(self, base_property, fstep, units, cube):
        self.base_property = base_property
        self.fstep = fstep
        self.units = units
        self.cube = cube
        self.Calculate()

    def Calculate(self):
        P = self.base_property
        f = self.fstep  
        indexes = ['x','y','z']

        # polarizability density
        for i in range(3):
            t = self.D_i(i,P,f,7)
            if isinstance(t, ndarray):
                self.cube.data = t
                self.cube.writeCubeFile('p_alpha_'+indexes[i]+'.cube')

        # first hyperpolarizability density
        for i in range(3):
            t = self.D_ij(i,i,P,f,7)
            if isinstance(t, ndarray):
                self.cube.data = t
                self.cube.writeCubeFile('p_beta_'+2*indexes[i]+'.cube')

        # second hyperpolarizability density
        for i in range(3):
            t = self.D_iijj(i,i,P,f,7)
            if isinstance(t, ndarray):
                self.cube.data = t
                self.cube.writeCubeFile('p_gamma_'+3*indexes[i]+'.cube')

    def Fi(self,i,bi):
        '''Select diagonal field.'''
        Fi = [0,0,0]
        Fi[i] = bi
        return tuple(Fi)

    def Fij(self,i,j,bi,bj):
        '''Select off-diagonal field.'''
        Fij = [0,0,0]
        Fij[i] = bi
        Fij[j] = bj
        return tuple(Fij)

    def D_i(self,i,E,f,np):
        '''Components of dipole moment (energy expansion).'''
        if np == 7:
            try:
                D_i= ( E[Fi(i,3*f)] - E[Fi(i,-3*f)] -
                       9.0*(E[Fi(i,2*f)] - E[Fi(i,-2*f)] +
                       5.0*(E[Fi(i,-f)] - E[Fi(i,f)]) ) ) / (60.0*abs(f))
            except KeyError:
                D_i=NaN

        return D_i

    def D_ij(self,i,j,E,f,np):
        '''Components of polarizability (energy expansion).'''
        if np == 7:
            try:
                if i==j:
                    D_ij=  ( 2.0*( E[Fi(i,3*f)] + E[Fi(i,-3*f)] ) -
                             27.0*( E[Fi(i,2*f)] + E[Fi(i,-2*f)] -
                                    10.0*( E[Fi(i,f)] + E[Fi(i,-f)] )) -
                             490.0*E[Fi(0,0)] ) / (360.0*abs(f*f))
                if i!=j:
                    D_ij=  ( 12.0*( E[Fij(i,j,3*f,f)] + E[Fij(i,j,f,3*f)] -
                                    E[Fij(i,j,3*f,0)] - E[Fij(i,j,0,3*f)] +
                                    E[Fij(i,j,-2*f,0)] + E[Fij(i,j,0,-2*f)] ) +
                             480.0*E[Fij(i,j,f,f)] +
                             330.0*( E[Fi(0,0)] - E[Fij(i,j,0,f)] - E[Fij(i,j,f,0)] ) +
                             160.0*E[Fij(i,j,-f,-f)] -
                             105.0*( E[Fij(i,j,f,2*f)] + E[Fij(i,j,2*f,f)] ) -
                             15.0*( E[Fij(i,j,-2*f,-f)] + E[Fij(i,j,-f,-2*f)] ) +
                             90.0*( E[Fij(i,j,2*f,0)] + E[Fij(i,j,0,2*f)] -
                                    E[Fij(i,j,-f,0)] - E[Fij(i,j,0,-f)] ) -
                             60.0*( E[Fij(i,j,f,-f)] + E[Fij(i,j,-f,f)] ) +
                              5.0*( E[Fij(i,j,2*f,-f)] + E[Fij(i,j,-f,2*f)] ) +
                              3.0*( E[Fij(i,j,-2*f,f)] + E[Fij(i,j,f,-2*f)] ) +
                              10.0*E[Fij(i,j,2*f,2*f)] ) / (360.0*abs(f*f))
            except KeyError:
                D_ij=NaN

        return D_ij

    def D_ijj(self,i,j,E,f,np):
        '''Components of hyperpolarizability (energy expansion).'''
        if np == 7:
            try:
                # Beta_iii
                if i==j:
                    D_ijj=  ( E[Fi(i,-3*f)] - E[Fi(i,3*f)] +
                              8.0*( E[Fi(i,2*f)] - E[Fi(i,-2*f)] ) -
                              13.0*( E[Fi(i,f)] - E[Fi(i,-f)] ) ) / (48.0*abs(f*f*f))
                # Beta_ijj
                if i!=j:
                    D_ijj=  ( E[Fij(i,j,-f,2*f)] + E[Fij(i,j,-f,-2*f)] -
                              E[Fij(i,j,f,2*f)] - E[Fij(i,j,f,-2*f)] +
                              2.0*( E[Fij(i,j,-2*f,f)] + E[Fij(i,j,-2*f,-f)] -
                                    E[Fij(i,j,2*f,f)] - E[Fij(i,j,2*f,-f)] ) +
                              4.0*( E[Fij(i,j,2*f,0)] - E[Fij(i,j,-2*f,0)] ) +
                              20.0*( E[Fij(i,j,f,f)] + E[Fij(i,j,f,-f)] -
                                     E[Fij(i,j,-f,f)] - E[Fij(i,j,-f,-f)] ) -
                              38.0*( E[Fij(i,j,f,0)] - E[Fij(i,j,-f,0)] ) ) / (144.0*abs(f*f*f))
            except KeyError:
                D_ijj=NaN

        return D_ijj

    def D_iijj(self,i,j,E,f,np):
        '''Components of second hyperpolarizability (energy expansion).'''
        if np == 7:
            try:
                # Gamma_iiii
                if i==j:
                    D_iijj=  ( 56.0*E[Fi(0,0)] - E[Fi(i,3*f)] - E[Fi(i,-3*f)] +
                               12.0*( E[Fi(i,2*f)] + E[Fi(i,-2*f)] ) -
                               39.0*( E[Fi(i,f)] + E[Fi(i,-f)] ) ) / (144.0*abs(f*f*f*f))
                # Gamma_iijj
                if i!=j:
                    D_iijj=  ( 72.0*E[Fi(0,0)] - E[Fij(i,j,-2*f,-f)] - E[Fij(i,j,-f,-2*f)] -
                               E[Fij(i,j,-f,2*f)] - E[Fij(i,j,f,-2*f)] -
                               E[Fij(i,j,f,2*f)] - E[Fij(i,j,2*f,-f)] -
                               E[Fij(i,j,2*f,f)] - E[Fij(i,j,-2*f,f)] +
                               2.0*( E[Fij(i,j,0,2*f)] + E[Fij(i,j,2*f,0)] +
                                     E[Fij(i,j,-2*f,0)] + E[Fij(i,j,0,-2*f)] ) +
                               20.0*( E[Fij(i,j,f,-f)] + E[Fij(i,j,f,f)] +
                                      E[Fij(i,j,-f,-f)] + E[Fij(i,j,-f,f)] ) -
                               38.0*( E[Fij(i,j,-f,0)] + E[Fij(i,j,0,-f)] +
                                      E[Fij(i,j,0,f)] + E[Fij(i,j,f,0)] ) ) / (288.0*abs(f*f*f*f))
            except KeyError:
                D_iijj=NaN

        return D_iijj



class FDIFF_D_V(FDIFF):

    def __init__(self, D, fstep, units, method):
        '''Initialize for dipole expansion'''
        self.base_name = 'dipole'

        FDIFF.__init__(self, D, fstep, units, method)

    def Mu(self,i,D,f):
        '''Components of dipole moment (dipole expansion).'''
        try:
            mu_i=( (2.0/3.0)*( D[Fi(i,f)][i] + D[Fi(i,-f)][i] ) - 
                   (1.0/6.0)*( D[Fi(i,2*f)][i] + D[Fi(i,-2*f)][i] ) )
        except KeyError:
            mu_i=NaN

        return mu_i

    def Alpha(self,i,j,D,f):
        '''Components of polarizability (dipole expansion).'''
        try:
            if i==j:
                alpha_ij=(  (2.0/3.0)*( D[Fi(i,f)][i] - D[Fi(i,-f)][i] ) -
                            (1.0/12.0)*( D[Fi(i,2*f)][i] - D[Fi(i,-2*f)][i] ) )/abs(f)
            if i!=j:
                alpha_ij=(  (2.0/3.0)*( D[Fi(j,f)][i] - D[Fi(j,-f)][i] ) -
                            (1.0/12.0)*( D[Fi(j,2*f)][i] - D[Fi(j,-2*f)][i] ) )/abs(f)
        except KeyError:
            alpha_ij=NaN

        return alpha_ij

    def Beta(self,i,j,D,f):
        '''Components of hyperpolarizability (dipole expansion).'''
        try:
            # Beta_iii
            if i==j:
                beta_ijj=( (1.0/3.0)*( D[Fi(i,2*f)][i] + D[Fi(i,-2*f)][i] -
                           D[Fi(i,f)][i] - D[Fi(i,-f)][i] ) )/abs(f*f)
            # Beta_ijj
            if i!=j:
                beta_ijj=( (1.0/3.0)*( D[Fi(j,2*f)][i] + D[Fi(j,-2*f)][i] -
                           D[Fi(j,f)][i] - D[Fi(j,-f)][i] ) )/abs(f*f)
        except KeyError:
            beta_ijj=NaN

        return beta_ijj

    def Gamma(self,i,j,D,f):
        '''Components of second hyperpolarizability (dipole expansion).'''
        try:
            # Gamma_iiii
            if i==j:
                gamma_iijj=( (1.0/2.0)*( D[Fi(i,2*f)][i] - D[Fi(i,-2*f)][i] ) -
                             D[Fi(i,f)][i] + D[Fi(i,-f)][i] )/abs(f*f*f)
            # Gamma_iijj
            if i!=j:
                gamma_iijj=( (1.0/2.0)*( D[Fij(i,j,f,f)][i] - D[Fij(i,j,-f,f)][i] +
                             D[Fij(i,j,f,-f)][i] - D[Fij(i,j,-f,-f)][i] ) -
                             D[Fi(i,f)][i] + D[Fi(i,-f)][i] )/abs(f*f*f)
        except KeyError:
            gamma_iijj=NaN

        return gamma_iijj

#----------------------------------------------------------------------------
# Property calculation routines (old)
#----------------------------------------------------------------------------

def fields(i,f):
    '''25 fields for average properties'''

    fields=['%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  0.0),
            '%7.4f %7.4f %7.4f'  % (  -f,  0.0,  0.0),
            '%7.4f %7.4f %7.4f'  % (   f,  0.0,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0,    f,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,   -f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,    f),
            '%7.4f %7.4f %7.4f'  % (-2*f,  0.0,- 0.0),
            '%7.4f %7.4f %7.4f'  % ( 2*f,  0.0,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0, -2*f,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0,  2*f,  0.0),
            '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0, -2*f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,  0.0,  2*f),
            '%7.4f %7.4f %7.4f'  % (  -f,   -f,  0.0),
            '%7.4f %7.4f %7.4f'  % (   f,   -f,  0.0),
            '%7.4f %7.4f %7.4f'  % (  -f,    f,  0.0),
            '%7.4f %7.4f %7.4f'  % (   f,    f,  0.0),
            '%7.4f %7.4f %7.4f'  % (  -f,  0.0,   -f),
            '%7.4f %7.4f %7.4f'  % (   f,  0.0,   -f),
            '%7.4f %7.4f %7.4f'  % (  -f,  0.0,    f),
            '%7.4f %7.4f %7.4f'  % (   f,  0.0,    f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,   -f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,    f,   -f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,   -f,    f),
            '%7.4f %7.4f %7.4f'  % ( 0.0,    f,    f)]

    field=fields[i]
    
    return field

def Fi(i,bi):
    '''Select diagonal field.'''
    Fi = [0,0,0]
    Fi[i] = bi
    return tuple(Fi)

def Fij(i,j,bi,bj):
    '''Select off-diagonal field.'''
    Fij = [0,0,0]
    Fij[i] = bi
    Fij[j] = bj
    return tuple(Fij)

def Mu_E(i,E,f):
    '''Components of dipole moment (energy expansion).'''
    try:
        mu_i=( (-2.0/3.0)*( E[Fi(i,f)] - E[Fi(i,-f)] ) + 
               ( E[Fi(i,2*f)] - E[Fi(i,-2*f)] )/12.0 )/abs(f)
    except KeyError:
        mu_i=NaN

    return mu_i

def Alpha_E(i,j,E,f):
    '''Components of polarizability (energy expansion).'''
    try:
        if i==j:
            alpha_ij=( (5.0/2.0)*E[Fi(0,0)] - (4.0/3.0)*( E[Fi(i,f)] + E[Fi(i,-f)] ) + 
                       ( E[Fi(i,2*f)] + E[Fi(i,-2*f)] )/12.0 )/abs(f*f)
        if i!=j:
            alpha_ij=( ( E[Fij(i,j,2*f,2*f)] - E[Fij(i,j,2*f,-2*f)] - 
                         E[Fij(i,j,-2*f,2*f)] + E[Fij(i,j,-2*f,-2*f)] )/48.0 -
                       ( E[Fij(i,j,f,f)] - E[Fij(i,j,f,-f)] -
                         E[Fij(i,j,-f,f)] + E[Fij(i,j,-f,-f)] )/3.0 )/abs(f*f)
    except KeyError:
        alpha_ij=NaN

    return alpha_ij

def Beta_E(i,j,E,f):
    '''Components of hyperpolarizability (energy expansion).'''
    try:
        # Beta_iii
        if i==j:
            beta_ijj=( ( E[Fi(i,f)] - E[Fi(i,-f)] ) -
                       ( E[Fi(i,2*f)] - E[Fi(i,-2*f)] )/2.0 )/abs(f*f*f)
        # Beta_ijj
        if i!=j:
            beta_ijj=( ( E[Fij(i,j,-f,-f)] - E[Fij(i,j,f,f)] +
                         E[Fij(i,j,-f,f)] - E[Fij(i,j,f,-f)] )/2.0 +
                         E[Fi(i,f)] - E[Fi(i,-f)] )/abs(f*f*f)
    except KeyError:
        beta_ijj=NaN

    return beta_ijj

def Gamma_E(i,j,E,f):
    '''Components of second hyperpolarizability (energy expansion).'''
    try:
        # Gamma_iiii
        if i==j:
            gamma_iijj=( 4.0*( E[Fi(i,f)] + E[Fi(i,-f)] ) - 6.0*E[Fi(0,0)] -
                         ( E[Fi(i,2*f)] + E[Fi(i,-2*f)] ) )/abs(f*f*f*f)
        # Gamma_iijj
        if i!=j:
            gamma_iijj=( -4.0*E[Fi(0,0)] - ( E[Fij(i,j,f,f)] + E[Fij(i,j,-f,-f)] +
                         E[Fij(i,j,f,-f)] + E[Fij(i,j,-f,f)] ) +
                         2.0*( E[Fi(i,f)] + E[Fi(i,-f)] + E[Fi(j,f)] + E[Fi(j,-f)] ) )/abs(f*f*f*f)
    except KeyError:
        gamma_iijj=NaN

    return gamma_iijj

def Mu_D(i,D,f):
    '''Components of dipole moment (dipole expansion).'''
    try:
        mu_i=( (2.0/3.0)*( D[Fi(i,f)][i] + D[Fi(i,-f)][i] ) - 
               (1.0/6.0)*( D[Fi(i,2*f)][i] + D[Fi(i,-2*f)][i] ) )
    except KeyError:
        mu_i=NaN

    return mu_i

def Alpha_D(i,j,D,f):
    '''Components of polarizability (dipole expansion).'''
    try:
        if i==j:
            alpha_ij=(  (2.0/3.0)*( D[Fi(i,f)][i] - D[Fi(i,-f)][i] ) -
                        (1.0/12.0)*( D[Fi(i,2*f)][i] - D[Fi(i,-2*f)][i] ) )/abs(f)
        if i!=j:
            alpha_ij=(  (2.0/3.0)*( D[Fi(j,f)][i] - D[Fi(j,-f)][i] ) -
                        (1.0/12.0)*( D[Fi(j,2*f)][i] - D[Fi(j,-2*f)][i] ) )/abs(f)
    except KeyError:
        alpha_ij=NaN

    return alpha_ij

def Beta_D(i,j,D,f):
    '''Components of hyperpolarizability (dipole expansion).'''
    try:
        # Beta_iii
        if i==j:
            beta_ijj=( (1.0/3.0)*( D[Fi(i,2*f)][i] + D[Fi(i,-2*f)][i] -
                       D[Fi(i,f)][i] - D[Fi(i,-f)][i] ) )/abs(f*f)
        # Beta_ijj
        if i!=j:
            beta_ijj=( (1.0/3.0)*( D[Fi(j,2*f)][i] + D[Fi(j,-2*f)][i] -
                       D[Fi(j,f)][i] - D[Fi(j,-f)][i] ) )/abs(f*f)
    except KeyError:
        beta_ijj=NaN

    return beta_ijj

def Gamma_D(i,j,D,f):
    '''Components of second hyperpolarizability (dipole expansion).'''
    try:
        # Gamma_iiii
        if i==j:
            gamma_iijj=( (1.0/2.0)*( D[Fi(i,2*f)][i] - D[Fi(i,-2*f)][i] ) -
                         D[Fi(i,f)][i] + D[Fi(i,-f)][i] )/abs(f*f*f)
        # Gamma_iijj
        if i!=j:
            gamma_iijj=( (1.0/2.0)*( D[Fij(i,j,f,f)][i] - D[Fij(i,j,-f,f)][i] +
                         D[Fij(i,j,f,-f)][i] - D[Fij(i,j,-f,-f)][i] ) -
                         D[Fi(i,f)][i] + D[Fi(i,-f)][i] )/abs(f*f*f)
    except KeyError:
        gamma_iijj=NaN

    return gamma_iijj

def CalcKurtzD(D, f, U, Properties, Label):
    '''Calculate properties using finite field method (dipole expansion)'''

    # calculate vector/tensor elements
    x = Mu_D(0,D,f)
    y = Mu_D(1,D,f)
    z = Mu_D(2,D,f)
    # pack tensor elements into an array
    M = array([x,y,z],dtype=float64)
    # vector norm
    Mv =sqrt(dot(M,M))

    Properties['M']=M
    Properties['Mv']=Mv

    # polarizability
    xx =Alpha_D(0,0,D,f)
    xy =Alpha_D(0,1,D,f)
    xz =Alpha_D(0,2,D,f)
    yx =Alpha_D(1,0,D,f)
    yy =Alpha_D(1,1,D,f)
    yz =Alpha_D(1,2,D,f)
    zx =Alpha_D(2,0,D,f)
    zy =Alpha_D(2,1,D,f)
    zz =Alpha_D(2,2,D,f)
    # pack tensor elements into an array
    A = array([xx,xy,xz,
               yx,yy,yz,
               zx,zy,zz],dtype=float64).reshape(3,3)
    # isotropic average
    Av = trace(A)/3.0

    Properties['A']=A
    Properties['Av']=Av

    # first hyperpolarizability
    xxx=Beta_D(0,0,D,f)
    xyy=Beta_D(0,1,D,f)
    xzz=Beta_D(0,2,D,f)
    yxx=Beta_D(1,0,D,f)
    yyy=Beta_D(1,1,D,f)
    yzz=Beta_D(1,2,D,f)
    zxx=Beta_D(2,0,D,f)
    zyy=Beta_D(2,1,D,f)
    zzz=Beta_D(2,2,D,f)
    # pack tensor elements into an array
    B = array([xxx,yxx,zxx,
               xyy,yyy,zyy,
               xzz,yzz,zzz],dtype=float64).reshape(3,3)
    # vector component
    Bx = (3.0/5.0)*B.sum(axis=0)[0]
    By = (3.0/5.0)*B.sum(axis=0)[1]
    Bz = (3.0/5.0)*B.sum(axis=0)[2]
    # projection to the dipole moment vector
    Bm = (3.0/5.0)*dot(M,B.sum(axis=0))/Mv

    Properties['B']=B
    Properties['Bx']=Bx
    Properties['By']=By
    Properties['Bz']=Bz
    Properties['Bm']=Bm

    # second hyperpolarizability
    xxxx=Gamma_D(0,0,D,f)
    xxyy=Gamma_D(0,1,D,f)
    xxzz=Gamma_D(0,2,D,f)
    yyxx=Gamma_D(1,0,D,f)
    yyyy=Gamma_D(1,1,D,f)
    yyzz=Gamma_D(1,2,D,f)
    zzxx=Gamma_D(2,0,D,f)
    zzyy=Gamma_D(2,1,D,f)
    zzzz=Gamma_D(2,2,D,f)
    # pack tensor elements into an array
    G = array([xxxx,xxyy,xxzz,
               yyxx,yyyy,yyzz,
               zzxx,zzyy,zzzz],dtype=float64).reshape(3,3)
    # isotropic average
    Gv = (G.trace()+2.0*(G[0,1]+G[0,2]+G[1,2]))/5.0

    Properties['G']=G
    Properties['Gv']=Gv

    # Prepare data for printout
    log  = 'Static electric dipole properties (finite field '+Label+' dipole, f=%.4f) \n\n' % f

    log += 'Dipole moment [%s]\n\n' % U['u']
    line = 4*U['m']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'x'.rjust(15), 'y'.rjust(15), 'z'.rjust(15), '<mu>'.rjust(15) )
    log += line % ( x, y, z, Mv )

    log += 'Polarizability [%s]\n\n' % U['u']
    line = 4*U['a']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xx'.rjust(15), 'xy'.rjust(15), 'xz'.rjust(15), '<alpha>'.rjust(15) )
    log += line % ( xx, xy, xz, Av )
    line = 3*U['a']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'yx'.rjust(15), 'yy'.rjust(15), 'yz'.rjust(15) )
    log += line % ( yx, yy, yz )
    log += '%15s %15s %15s\n' % ( 'zx'.rjust(15), 'zy'.rjust(15), 'zz'.rjust(15) )
    log += line % ( zx, zy, zz )

    log += 'First Hyperpolarizability [%s]\n\n' % U['u']
    line = 4*U['b']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xxx'.rjust(15), 'yxx'.rjust(15), 'zxx'.rjust(15), '<beta_mu>'.rjust(15) )
    log += line % ( xxx, yxx, zxx, Bm )
    line = 3*U['b']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'xyy'.rjust(15), 'yyy'.rjust(15), 'zyy'.rjust(15) )
    log += line % ( xyy, yyy, zyy )
    log += '%15s %15s %15s\n' % ( 'xzz'.rjust(15), 'yzz'.rjust(15), 'zzz'.rjust(15) )
    log += line % ( xzz, yzz, zzz )
    log += '%15s %15s %15s\n' % ( '<beta_x>'.rjust(15), '<beta_y>'.rjust(15), '<beta_z>'.rjust(15) )
    log += line % ( Bx, By, Bz )

    log += 'Second Hyperpolarizability [%s]\n\n' % U['u']
    line = 4*U['g']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xxxx'.rjust(15), 'xxyy'.rjust(15), 'xxzz'.rjust(15), '<gamma>'.rjust(15) )
    log += line % ( xxxx, xxyy, xxzz, Gv )
    line = 3*U['g']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'yyxx'.rjust(15), 'yyyy'.rjust(15), 'yyzz'.rjust(15) )
    log += line % ( yyxx, yyyy, yyzz )
    log += '%15s %15s %15s\n' % ( 'zzxx'.rjust(15), 'zzyy'.rjust(15), 'zzzz'.rjust(15) )
    log += line % ( zzxx, zzyy, zzzz )

    print log

def CalcKurtzE(E, f, U, Properties, Label):
    '''Calculate properties using finite field method (energy expansion)'''

    # calculate vector/tensor elements
    x = Mu_E(0,E,f)
    y = Mu_E(1,E,f)
    z = Mu_E(2,E,f)
    # pack tensor elements into an array
    M = array([x,y,z],dtype=float64)
    # vector norm
    Mv =sqrt(dot(M,M))

    Properties['M']=M
    Properties['Mv']=Mv

    # polarizability
    xx =Alpha_E(0,0,E,f)
    xy =Alpha_E(0,1,E,f)
    xz =Alpha_E(0,2,E,f)
    yx =Alpha_E(1,0,E,f)
    yy =Alpha_E(1,1,E,f)
    yz =Alpha_E(1,2,E,f)
    zx =Alpha_E(2,0,E,f)
    zy =Alpha_E(2,1,E,f)
    zz =Alpha_E(2,2,E,f)
    # pack tensor elements into an array
    A = array([xx,xy,xz,
               yx,yy,yz,
               zx,zy,zz],dtype=float64).reshape(3,3)
    # isotropic average
    Av = trace(A)/3.0

    Properties['A']=A
    Properties['Av']=Av

    # first hyperpolarizability
    xxx=Beta_E(0,0,E,f)
    xyy=Beta_E(0,1,E,f)
    xzz=Beta_E(0,2,E,f)
    yxx=Beta_E(1,0,E,f)
    yyy=Beta_E(1,1,E,f)
    yzz=Beta_E(1,2,E,f)
    zxx=Beta_E(2,0,E,f)
    zyy=Beta_E(2,1,E,f)
    zzz=Beta_E(2,2,E,f)
    # pack tensor elements into an array
    B = array([xxx,yxx,zxx,
               xyy,yyy,zyy,
               xzz,yzz,zzz],dtype=float64).reshape(3,3)
    # vector component
    Bx = (3.0/5.0)*B.sum(axis=0)[0]
    By = (3.0/5.0)*B.sum(axis=0)[1]
    Bz = (3.0/5.0)*B.sum(axis=0)[2]
    # projection to the dipole moment vector
    Bm = (3.0/5.0)*dot(M,B.sum(axis=0))/Mv

    Properties['B']=B
    Properties['Bx']=Bx
    Properties['By']=By
    Properties['Bz']=Bz
    Properties['Bm']=Bm

    # second hyperpolarizability
    xxxx=Gamma_E(0,0,E,f)
    xxyy=Gamma_E(0,1,E,f)
    xxzz=Gamma_E(0,2,E,f)
    yyxx=Gamma_E(1,0,E,f)
    yyyy=Gamma_E(1,1,E,f)
    yyzz=Gamma_E(1,2,E,f)
    zzxx=Gamma_E(2,0,E,f)
    zzyy=Gamma_E(2,1,E,f)
    zzzz=Gamma_E(2,2,E,f)
    # pack tensor elements into an array
    G = array([xxxx,xxyy,xxzz,
               yyxx,yyyy,yyzz,
               zzxx,zzyy,zzzz],dtype=float64).reshape(3,3)
    # isotropic average
    Gv = (G.trace()+2.0*(G[0,1]+G[0,2]+G[1,2]))/5.0

    Properties['G']=G
    Properties['Gv']=Gv

    # Prepare data for printout
    log  = 'Static electric dipole properties (finite field '+Label+' energy, f=%.4f) \n\n' % f

    log += 'Dipole moment [%s]\n\n' % U['u']
    line = 4*U['m']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'x'.rjust(15), 'y'.rjust(15), 'z'.rjust(15), '<mu>'.rjust(15) )
    log += line % ( x, y, z, Mv )

    log += 'Polarizability [%s]\n\n' % U['u']
    line = 4*U['a']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xx'.rjust(15), 'xy'.rjust(15), 'xz'.rjust(15), '<alpha>'.rjust(15) )
    log += line % ( xx, xy, xz, Av )
    line = 3*U['a']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'yx'.rjust(15), 'yy'.rjust(15), 'yz'.rjust(15) )
    log += line % ( yx, yy, yz )
    log += '%15s %15s %15s\n' % ( 'zx'.rjust(15), 'zy'.rjust(15), 'zz'.rjust(15) )
    log += line % ( zx, zy, zz )

    log += 'First Hyperpolarizability [%s]\n\n' % U['u']
    line = 4*U['b']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xxx'.rjust(15), 'yxx'.rjust(15), 'zxx'.rjust(15), '<beta_mu>'.rjust(15) )
    log += line % ( xxx, yxx, zxx, Bm )
    line = 3*U['b']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'xyy'.rjust(15), 'yyy'.rjust(15), 'zyy'.rjust(15) )
    log += line % ( xyy, yyy, zyy )
    log += '%15s %15s %15s\n' % ( 'xzz'.rjust(15), 'yzz'.rjust(15), 'zzz'.rjust(15) )
    log += line % ( xzz, yzz, zzz )
    log += '%15s %15s %15s\n' % ( '<beta_x>'.rjust(15), '<beta_y>'.rjust(15), '<beta_z>'.rjust(15) )
    log += line % ( Bx, By, Bz )

    log += 'Second Hyperpolarizability [%s]\n\n' % U['u']
    line = 4*U['g']['f'] + '\n\n'
    log += '%15s %15s %15s %15s\n' % ( 'xxxx'.rjust(15), 'xxyy'.rjust(15), 'xxzz'.rjust(15), '<gamma>'.rjust(15) )
    log += line % ( xxxx, xxyy, xxzz, Gv )
    line = 3*U['g']['f'] + '\n\n'
    log += '%15s %15s %15s\n' % ( 'yyxx'.rjust(15), 'yyyy'.rjust(15), 'yyzz'.rjust(15) )
    log += line % ( yyxx, yyyy, yyzz )
    log += '%15s %15s %15s\n' % ( 'zzxx'.rjust(15), 'zzyy'.rjust(15), 'zzzz'.rjust(15) )
    log += line % ( zzxx, zzyy, zzzz )

    print log

#----------------------------------------------------------------------------
# Units, conversion factors and formats
#----------------------------------------------------------------------------

def SetUnits(u):
    '''Set units of electric properties and the respective formats.'''

    units = {}

    if u.lower() == 'au':
        units['u'] = 'a.u.'
        units['m'] = {'x':1.0,          'f':'%15.4f ', 'r':4}
        units['a'] = {'x':1.0,          'f':'%15.3f ', 'r':3}
        units['b'] = {'x':1.0,          'f':'%15.2f ', 'r':2}
        units['g'] = {'x':1.0,          'f':'%15.1f ', 'r':1}
    elif u.lower() == 'si':
        units['u'] = 'si'
        units['m'] = {'x':8.478358e-30, 'f':'%15.5e ', 'r':5} # C m
        units['a'] = {'x':1.648778e-41, 'f':'%15.5e ', 'r':5} # C^2 m^2 J^-1
        units['b'] = {'x':3.206361e-53, 'f':'%15.5e ', 'r':5} # C^3 m^3 J^-2
        units['g'] = {'x':6.235377e-65, 'f':'%15.5e ', 'r':5} # C^4 m^4 J^-3
    elif u.lower() == 'asi':
        units['u'] = 'asi'
        units['m'] = {'x':8.4784e-30,   'f':'%15.5e ', 'r':5} # C m
        units['a'] = {'x':1.8621e-30,   'f':'%15.5e ', 'r':5} # m^3
        units['b'] = {'x':3.6213e-42,   'f':'%15.5e ', 'r':5} # m^4 V^-1
        units['g'] = {'x':7.0423e-54,   'f':'%15.5e ', 'r':5} # m^5 V^-2
    elif u.lower() == 'esu':
        units['u'] = 'esu'
        units['m'] = {'x':2.5418e-18,   'f':'%15.5e ', 'r':5} # statvolt cm^2
        units['a'] = {'x':1.4817e-25,   'f':'%15.5e ', 'r':5} # cm^3
        units['b'] = {'x':8.6392e-33,   'f':'%15.5e ', 'r':5} # statvolt^-1 cm^4
        units['g'] = {'x':5.0367e-40,   'f':'%15.5e ', 'r':5} # statvolt^-2 cm^5
    else:
        print "Unknown units requested!"
        sys.exit(1)

    # Selected tensor elements
    #                  Mu:     Alpha:     Beta:         Gamma:
    # [[00, 01, 02],   x y z   xx xy xz   xxx yxx zxx   xxxx xxyy xxzz
    #  [10, 11, 12],           yx yy yz   xyy yyy zyy   yyxx yyyy yyzz
    #  [20, 21, 22]]           zx zy zz   xzz yzz zzz   zzxx zzyy zzzz
    # 
    property_index = {
        'Mu'    :  [    'x',    'y',    'z' ],
        'Alpha' : [[   'xx',   'xy',   'xz' ],
                   [   'yx',   'yy',   'yz' ],
                   [   'zx',   'zy',   'zz' ]],
        'Beta'  : [[  'xxx',  'yxx',  'zxx' ],
                   [  'xyy',  'yyy',  'zyy' ],
                   [  'xzz',  'yzz',  'zzz' ]],
        'Gamma' : [[ 'xxxx', 'xxyy', 'xxzz' ],
                   [ 'yyxx', 'yyyy', 'yyzz' ],
                   [ 'zzxx', 'zzyy', 'zzzz' ]] }

    # Descriptions of properties
    property_description = {
        'Mu'    : '# Dipole Moment Vector i=' + \
                     property_index['Mu'][2],
        '|D|'   : '# Total Dipole Moment',
        'Alpha' : '# Polarizability Tensor i=' + \
                     property_index['Alpha'][2][2],
        '<A>'   : '# Isotropic Polatizability',
        '<B>'   : '# Anisotropy of Polarizability (Z-axis is the rotation axis)',
        'Beta'  : '# First hyperpolarizability tensor i=' + \
                     property_index['Beta'][2][2],
        'B(Z)'  : '# Vector component of hyperpolarizability tensor (Z is the permament dipole moment direction)',
        'Gamma' : '# Second hyperpolarizability tensor i=' + \
                     property_index['Gamma'][2][2],
        '<G>'   : '# Scalar component of second hyperpolarizability tensor given by the isotropic average' }

    return units

def Periodic(mendeleiev):
    '''Returns the mendeleiev table as a python list of tuples. Each cell
    contains either None or a tuple (symbol, atomic number), or a list of pairs
    for the cells * and **. Requires: "import re". Source: Gribouillis at
    www.daniweb.com - 2008 '''

    # L is a consecutive list of tuples ('Symbol', atomic number)
    L = [ (e,i+1) for (i,e) in enumerate( re.compile ("[A-Z][a-z]*").findall('''
    HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKr
    RbSrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoEr
    TmYbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFm
    MdNoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'''))]

    # The following fills the void with nones and returns the list of lists
    mendeleiev = 0

    if mendeleiev:
        for i,j in ( (88,103), (56,71) ):
            L[i] = L[i:j]
            L[i+1:] = L[j:]
        for i,j in ( (12,10 ), (4,10), (1,16) ):
            L[i:i]=[None]*j 

        return [ L[18*i:18*(i+1)] for i in range (7) ]

    # Return a plain list of tuples
    else:
        return L

def Atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''

    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)

class CUBE:
    def __init__(self,filename):
        """Represents a Gaussian03 cubefile object.
        Lifted from Noel O'Boyle Cubefile class

        Attributes:
            data -- a scipy array containing the data
            numAtoms -- the number of atoms
            origin -- the origin
            numPoints -- a list of the number of points in x,y,z
            spacing -- a list of the spacing in x,y,z
            dv -- the volume of chunk
        """

        self.readCubeFile(filename)
        self.CubeHeader()
        #print self.header
        #print self.numPoints
        #print self.spacing
        #print self.dv
        #print self.origin

    def I_P_dv(self):
        self.P = self.dv*self.data.sum()
        return self.P

    def I_rP_dv(self):

        self.rP = zeros(3)

        for i in range(self.numPoints[0]):
            for j in range(self.numPoints[1]):
                for k in range(self.numPoints[2]):
                    x,y,z = self.xyz(i,j,k)
                    self.rP[0] += x*self.data[j,i,k]*self.dv
                    self.rP[1] += y*self.data[j,i,k]*self.dv
                    self.rP[2] += z*self.data[j,i,k]*self.dv

        return self.rP

    def xyz(self,ix,iy,iz):
        x = self.origin[0]+(ix)*self.spacing[0]
        y = self.origin[1]+(iy)*self.spacing[1]
        z = self.origin[2]+(iz)*self.spacing[2]
        return array([x,y,z])

    def readCubeFile(self,filename):
        """Read in a cube file."""

        inputfile = open(filename,"r")
        self.head = [inputfile.readline(),inputfile.readline()]

        temp = inputfile.readline().strip().split()
        self.numAtoms = int(temp[0])
        self.origin = map(float,temp[1:])
        
        self.numPoints = [0]*3
        self.spacing = [0]*3
        self.dv=1.0
        for i in range(3):
            line = inputfile.readline().strip().split()
            self.numPoints[i] = int(line[0])
            temp = map(float,line[1:])
            self.spacing[i] = temp[i]
            self.dv *= temp[i]
            assert sum(temp[:i]+temp[i+1:])==0
    
        # Read in the lines with atom data
        self.atoms = ''
        for i in range(self.numAtoms):
            self.atoms += inputfile.readline()
        
        self.data = zeros( (self.numPoints[1],self.numPoints[0],self.numPoints[2]),"float" )
        i = j = k =0
        while i<self.numPoints[0]:
            line = inputfile.next()
            temp = map(float,line.strip().split())
            for x in range(0,len(temp)):
                self.data[j,i,x+k] = temp[x]

            k+=len(temp)
            if k==self.numPoints[2]:
                j+=1
                k = 0
                if j==self.numPoints[1]:
                    i+=1
                    j = 0
                    
        inputfile.close()

    def CubeHeader(self):
        """Read in a cube file."""

        log=self.head[0]+self.head[1]

        line = 3*'%12.6f'+'\n'
        log += '%5d' % self.numAtoms
        log += line % tuple(self.origin)

        for i in range(3):
            temp = [0.0, 0.0, 0.0]
            temp[i] = self.spacing[i]
            log += '%5d' % self.numPoints[i]
            log += line % tuple(temp)

        log += self.atoms

        self.header = log

    def writeCubeFile(self,filename):
        """Read in a cube file."""

        outputfile = open(filename,"w")

        log=self.head[0]+self.head[1]

        line = 3*'%12.6f'+'\n'
        log += '%5d' % self.numAtoms
        log += line % tuple(self.origin)

        for i in range(3):
            temp = [0.0, 0.0, 0.0]
            temp[i] = self.spacing[i]
            log += '%5d' % self.numPoints[i]
            log += line % tuple(temp)

        log += self.atoms

        #savetxt(filename, self.data, fmt='%13.5e')
        #self.data.tofile(outputfile, sep='\n', format='%13.5e')
        #self.data = zeros( (self.numPoints[1],self.numPoints[0],self.numPoints[2]),"float" )

        i = j = k =0
        nlines = int(self.numPoints[2]/6)
        ntails = self.numPoints[2] - 6*nlines

        while i < self.numPoints[0]:

            line = 6*'%13.5e'+'\n'
            for l in range(nlines):
                log += line % tuple(self.data[j,i,k:k+6])
                k+=6

            if ntails > 0:
                line = ntails*'%13.5e'+'\n'
                log += line % tuple(self.data[j,i,k:l+k])
                k+=ntails

            if k==self.numPoints[2]:
               j+=1
               k = 0
               if j==self.numPoints[1]:
                   i+=1
                   j = 0
                    
        outputfile.write(log)
        outputfile.close()

#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__":
    Main(sys.argv[1:])

