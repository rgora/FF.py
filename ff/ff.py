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
  -f, --base-field set base field
  -s, --gaussian   prepare gaussian inputs (expects gaussian.tmpl)
  -g, --gamess     prepare gamess inputs (expects gamess.tmpl)
  -m, --molcas     prepare molcas inputs (expects molcas.tmpl)
  -u, --units      set units; chose from: esu, si, asi 
                   (si multiplied by electric permittivity of free space),
                   or au which is the default
  --ff25           prepares ff25.inp file for RZ's ff25.f
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

from numpy import *
from string import Template

import os, sys, getopt, re

# Regular expressions
reflags = re.DOTALL

#----------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------
def usage():
    '''Print usage information.'''
    print __doc__
    print "Machine epsilon is: ",finfo(float64).eps,"for float64 type\n"


#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------
def Main(argv):
    '''Parse commandline and loop throught the logs'''

    Data = {}

    # Set up defaults
    fstep = 0.001
    calculate = 0
    gamess = 0
    gaussian = 0
    molcas = 0
    native = 0
    ff25 = 0
    units=SetUnits('au')

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hu:f:cgsmn",
                                        ["help",
                                         "units=",
                                         "base-field=",
                                         "calculate",
                                         "gamess",
                                         "gaussian",
                                         "molcas",
                                         "molcas-native",
                                         "ff25"
                                         ])
    except getopt.GetoptError, error:
        print(error)
        usage()
        sys.exit(2)
    if not argv:
        usage()
        sys.exit(0)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-u", "--units="):
            Units=SetUnits(arg)
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
        elif opt in ("--ff25"):
            ff25=1

    DataFiles = args
    # Parse each data file/dir
    for f in DataFiles:

        # Analyze outputs
        if calculate:

            if gamess:
                Data[f]=GAMESS(f, fstep,units)
            if molcas:
                Data[f]=MOLCAS(f,-fstep,units)
            if gaussian:
                Data[f]=GAUSSIAN(f,-fstep,units)

        # Prepare inputs
        else:
            if gaussian:
                GaussianInputs(f,fstep)
            if gamess:
                GamessInputs(f,fstep)
            if molcas:
                MolcasInputs(f,fstep)
            if native:
                MolcasNativeInputs(f,fstep)

#----------------------------------------------------------------------------
# Common parser routines
#----------------------------------------------------------------------------
class Parser:
    """A common parser routines"""

    def __init__(self, logpath, fstep, units):
        # initialize
        self.logpath = logpath
        self.fstep = fstep
        self.units = units
        self.ext = "\.log"
        # input data
        self.energies = {}
        self.dipoles = {}
        # results
        self.properties = {'en':{}, 'dm':{}}
        # useful regexp's
        self.re_number=re.compile(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?")
        # parse files and calculate properties
        self.parsefiles()
        self.calculate()

    def parsefiles(self):
        '''read logfiles in logpath'''
    
        self.fstep = 0

        logfiles=[]
        for f in os.listdir(self.logpath):
            if re.search(self.ext, f):
                logfiles.append(self.logpath+'/'+f)
        logfiles.sort()
    
        for log in logfiles:
            self.parsefile(log)
    
    def parsefile(self,filename):
        pass

    def calculate(self):
        for e in self.energies.keys():
            self.properties['en'][e]={}
            CalcKurtzE(self.energies[e], self.fstep, self.units, self.properties['en'][e], e)
        for d in self.dipoles.keys():
            self.properties['dm'][d]={}
            CalcKurtzD(self.dipoles[d], self.fstep, self.units, self.properties['dm'][d], d)

    def sortfields(self):

        self.fields = sorted( sorted( sorted( self.energies.keys(),
            lambda a,b: cmp(abs(a[0]), abs(b[0])) ),
            lambda a,b: cmp(abs(a[1]), abs(b[1])) ),
            lambda a,b: cmp(abs(a[2]), abs(b[2])) )

        self.nfields=(len(SFields)-1)/3/2
        self.fstep=abs(SFields[1][0])

        return self.fields, self.nfields, self.fstep

    def path(self):
        return self.logpath

class InputTemplate(Template):
    delimiter = '@'

class Inputs:
    """Common input routines"""
    def __init__(self, data, fstep):
        self.data = data
        self.fstep = fstep
        self.ReadTemplate()

    def ReadTemplate(self):
        """Read or punch standard template"""
        try:
            self.tmpl = open(self.pkg+'.tmpl','r').read()
            self.tmpl = InputTemplate(self.tmpl)
            self.WriteInputs()
        except IOError:
            print "There's no " + self.pkg + " template. I'm punching one - please check"
            open(self.pkg+'.tmpl','w').write(self.tmpl)
            sys.exit()

    def WriteInputs(self):
        pass

#----------------------------------------------------------------------------
# EDS Gamess (US) routines
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Gamess (US) routines
#----------------------------------------------------------------------------

class GamessInputs(Inputs):
    """Gamess US input routines"""

    def __init__(self, data, fstep):
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
        Inputs.__init__(self, data, fstep)

    def WriteInputs(self):

        # initialize periodic table
        p=periodic(0)
        
        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)
        
        for i in range(len(xyz)):
            xyz[i]=xyz[i].split()
            xyz[i].insert(1, str( atomn(xyz[i][0], p) ))
            xyz[i]='%-5s %5s %15s %15s %15s\n' % tuple([xyz[i][j] for j in range(5)])
        
        xyz.insert(0, ' $data\nFinite Field\nc1 0\n')
        xyz.append(' $end')
        xyz=''.join(xyz)
        
        # write input files
        for f in range(25):

            # prepare onefld inputs
            filename = self.data.replace('.xyz','')+'_F%2d_' % f
            filename = filename.replace(' ','0')
            ffield = ','.join(fields(f,self.fstep).split())
            finput = self.tmpl.substitute(data=xyz, fstep=self.fstep, field=ffield, onefld='.t.')
            if f==0:
                finput = finput.replace(' $efield','!$efield')

            # write onefld inputs
            open(filename+'.inp','w').write(finput)

            # write kurtz's ffield input
            if f==0:
                finput = self.tmpl.substitute(data=xyz, fstep=self.fstep, field=ffield, onefld='.f.')
                filename = self.data.replace('.xyz','')+'_%.4f' % self.fstep
                filename = filename.replace(' ','_')
                open(filename+'.inp','w').write(finput)

class GAMESS(Parser):
    """A GAMESS log parser."""

    def parsefile(self,filename):
        """Read gamess(us) energies for this system."""

        self.File = open(filename)

        # Molcas conversion factor
        self.au2d=2.541766

        # Read control options
        try:
            line = FindLine(self.File,'$CONTRL OPTIONS')
            line = SkipLines(self.File,2)
            self.scftyp = re.split('\s+|=',line)[2].lower()
            self.runtyp = re.split('\s+|=',line)[4].lower()
            line = SkipLines(self.File,1)
            self.mplevl = int(re.split('\s+|=',line)[3])
            self.cctyp = re.split('\s+|=',line)[9].upper()
            if self.mplevl == 2:
                line = FindLine(self.File,'MP2 CONTROL INFORMATION')
                line = SkipLines(self.File,5)
                self.mp2prp = re.split('\s+|=',line)[4].lower()
            
            if self.runtyp == 'eds':
                self.parse_eds()
            else:
                self.parse_gms()
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise
            print "File: ", filename, " is not a Gamess US file"

    def parse_gms(self):
        """Read gamess(us) energies for this system."""

        File = self.File

        # Find field, energies and dipoles
        while 1:
            line = File.readline()
            if line == '': break
    
            # Case of FF calculations
            if line.find('APPLIED FIELD') !=-1:
    
                # ... determine best energy ...
                if self.mplevl == 2:
                    EnKey = 'MP2'
                elif self.cctyp != 'NONE':
                    EnKey = cctyp
                else:
                    EnKey = 'SCF'
                if EnKey not in self.energies:
                    self.energies[EnKey]={}
    
                line = line.split()

                # ... set field label ...
                Field = ( round(float(line[-3]),4), round(float(line[-2]),4), round(float(line[-1]),4) )

                # ... set base field ...
                if ( (self.fstep == 0 and abs(Field[0]) > 0) or
                     (abs(Field[0]) > 0 and abs(Field[0]) < self.fstep) ):
                    self.fstep = abs(Field[0])

                if Field not in self.energies[EnKey]:
                    self.energies[EnKey][Field]= float64(SkipLines(File,1).split()[-1])
    
            # Case of single field calculations
            if line.find('ELECTRIC FIELD') !=-1:
                line = line.split()
                # ... set field label ...
                Field = ( round(float(line[-3]),4), round(float(line[-2]),4), round(float(line[-1]),4) )
                # ... set base field ...
                if ( (self.fstep == 0 and abs(Field[0]) > 0) or
                     (abs(Field[0]) > 0 and abs(Field[0]) < self.fstep) ):
                    self.fstep = abs(Field[0])
                # ... read scf energy ...
                line = FindLine(File,'SCF CALCULATION')
                if 'SCF' not in self.energies:
                    self.energies['SCF']={}
                if 'SCF' not in self.dipoles:
                    self.dipoles['SCF']={}
                if Field not in self.energies['SCF']:
                    self.energies['SCF'][Field]= float64(FindLine(File,' FINAL ').split()[4])
                # ... read scf dipole ...
                if Field not in self.dipoles['SCF']:
                    line = FindLine(File,'ELECTROSTATIC MOMENTS')
                    line = SkipLines(File,6)
                    self.dipoles['SCF'][Field] = array(line.split()[:3],dtype=float64)/self.au2d
                # ... read scf energy ...
                if self.mplevl == 2:
                    if 'MP2' not in self.energies:
                        self.energies['MP2']={}
                    if Field not in self.energies['MP2']:
                        self.energies['MP2'][Field]= float64(FindLine(File,'E(MP2)=').split()[-1])
                if self.mplevl == 2 and self.mp2prp == 't':
                    if 'MP2' not in self.dipoles:
                        self.dipoles['MP2']={}
                    if Field not in self.dipoles['MP2']:
                        line = FindLine(File,'MP2 PROPERTIES')
                        line = FindLine(File,'ELECTROSTATIC MOMENTS')
                        line = SkipLines(File,6)
                        self.dipoles['MP2'][Field] = array(line.split()[:3],dtype=float64)/self.au2d

    def parse_eds(self):
        """Read gamess(us) eds energies for this system."""

        File = self.File

        # Check onefld value
        line = FindLine('FFCALC INPUT OPTIONS') !=-1
        line = SkipLines(File,4)
        onefld = re.split('\s+|=',line)[2].lower()

        if onefld == 't':

            # Find field, interaction energies and total energies and dipoles
            while 1:
                line = File.readline()
                if line == '': break
            
                # Case of single field calculations
                if line.find('APPLIED FIELD') !=-1:
                    line = line.split()
                    # ... set field label ...
                    Field = ( round(float(line[-3]),4), round(float(line[-2]),4), round(float(line[-1]),4) )
                    # ... set base field ...
                    if ( (self.fstep == 0 and abs(Field[0]) > 0) or
                         (abs(Field[0]) > 0 and abs(Field[0]) < self.fstep) ):
                        self.fstep = abs(Field[0])
                    # ... read scf energy ...
                    line = FindLine(File,'SCF CALCULATION')
                    if 'SCF' not in self.energies:
                        self.energies['SCF']={}
                    if 'SCF' not in self.dipoles:
                        self.dipoles['SCF']={}
                    if Field not in self.energies['SCF']:
                        self.energies['SCF'][Field]= float64(FindLine(File,' FINAL ').split()[4])
                    # ... read scf dipole ...
                    if Field not in self.dipoles['SCF']:
                        line = FindLine(File,'ELECTROSTATIC MOMENTS')
                        line = SkipLines(File,6)
                        self.dipoles['SCF'][Field] = array(line.split()[:3],dtype=float64)/self.au2d
                    # ... read scf energy ...
                    if self.mplevl == 2:
                        if 'MP2' not in self.energies:
                            self.energies['MP2']={}
                        if Field not in self.energies['MP2']:
                            self.energies['MP2'][Field]= float64(FindLine(File,'E(MP2)=').split()[-1])
                    if self.mplevl == 2 and self.mp2prp == 't':
                        if 'MP2' not in self.dipoles:
                            self.dipoles['MP2']={}
                        if Field not in self.dipoles['MP2']:
                            line = FindLine(File,'MP2 PROPERTIES')
                            line = FindLine(File,'ELECTROSTATIC MOMENTS')
                            line = SkipLines(File,6)
                            self.dipoles['MP2'][Field] = array(line.split()[:3],dtype=float64)/self.au2d

#----------------------------------------------------------------------------
# Molcas routines
#----------------------------------------------------------------------------

class MolcasInputs(Inputs):
    """Molcas input routines"""

    def __init__(self, data, fstep):
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
        Inputs.__init__(self, data, fstep)

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
        for f in range(25):
            filename = self.data.replace('.xyz','')+'_F%2d_' % f
            filename = filename.replace(' ','0')

            # Prepare FFPT
            FFPT=' &FFPT\nDIPO\n'
            F=fields(f,self.fstep).split()
            FFPT+='X ' + F[0] + '\n'
            FFPT+='Y ' + F[1] + '\n'
            FFPT+='Z ' + F[2] + '\n'

            finput = self.tmpl.substitute(data=xyz, field=FFPT)
            open(filename+'.inp','w').write(finput)

class MolcasNativeInputs(Inputs):
    """Molcas input routines"""

    def __init__(self, data, fstep):
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
        MolcasInputs.__init__(self, data, fstep)

    def MakeCoords(self):
        '''read native xyz file'''
        try:
            xyz=open(self.data,'r').readlines()
            xyz=''.join(xyz).rstrip()
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)

        return xyz

class MOLCAS(Parser):
    """A MOLCAS log parser."""

    def parsefile(self,filename):
        '''Read molcas energies for this system.'''

        File=open(filename)

        # Molcas conversion factor
        au2d=2.5417463077102775

        # Find field, energies and dipoles
        while 1:
            line = File.readline()
            if line == '': break

            if line.find('FFPT    DIPO    COMP') !=-1:
                line = line.split()
                # ... set field label ...
                Field = ( round(float(line[-5]),4), round(float(line[-3]),4), round(float(line[-1]),4) )
                # ... set base field ...
                if ( (self.fstep == 0 and abs(Field[0]) > 0) or
                     (abs(Field[0]) > 0 and abs(Field[0]) < self.fstep) ):
                    self.fstep = abs(Field[0])
                # ... read all energies for this field ...
                while 1:
                    line = File.readline()
                    if line == '': break
                    if line.find('MOLCAS executing module FFPT') !=-1: break
                    # ... read scf energy and dipole ...
                    if line.find('MOLCAS executing module SCF') !=-1:
                        if ('SCF' not in self.energies):
                            self.energies['SCF']={}
                        if (Field not in self.energies['SCF']):
                            self.energies['SCF'][Field]= float64(ReFindLine(File,'Total .+ energy').split()[-1])
                        if ('SCF' not in self.dipoles):
                            self.dipoles['SCF']={}
                        if (Field not in self.dipoles['SCF']):
                            line = FindLine(File,'Dipole Moment')
                            dipl = SkipLines(File,2).split()
                            self.dipoles['SCF'][Field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    # ... read rasscf energy and dipole ...
                    if line.find('MOLCAS executing module RASSCF') !=-1:
                        if ('RASSCF' not in self.energies):
                            self.energies['RASSCF']={}
                        if (Field not in self.energies['RASSCF']):
                            line = FindLine(File,'Final state energy')
                            line = SkipLines(File,3)
                            self.energies['RASSCF'][Field] = float64(self.re_number.findall(line)[-1])
                        if ('RASSCF' not in self.dipoles):
                            self.dipoles['RASSCF']={}
                        if (Field not in self.dipoles['RASSCF']):
                            line = FindLine(File,'Dipole Moment')
                            dipl = SkipLines(File,2).split()
                            self.dipoles['RASSCF'][Field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    # ... read caspt2 energy and dipole ...
                    if line.find('MOLCAS executing module CASPT2') !=-1:
                        if ('CASPT2' not in self.energies):
                            self.energies['CASPT2']={}
                        if (Field not in self.energies['CASPT2']):
                            line = FindLine(File,'FINAL CASPT2 RESULT')
                            self.energies['CASPT2'][Field]= float64(SkipLines(File,5).split()[-1])
                        if ('CASPT2' not in self.dipoles):
                            self.dipoles['CASPT2']={}
                        if (Field not in self.dipoles['CASPT2']):
                            line = FindLine(File,'Dipole Moment')
                            dipl = SkipLines(File,2).split()
                            self.dipoles['CASPT2'][Field] = array([dipl[1],dipl[3],dipl[5]],dtype=float64)/au2d
                    # ... read ccsdt energy ...
                    if line.find('MOLCAS executing module CCSD(T)') !=-1:
                        if ('CCSD' not in self.energies):
                            self.energies['CCSD']={}
                        if (Field not in self.energies['CCSD']):
                            line = FindLine(File,'Total energy (diff)')
                            self.energies['CCSD'][Field]= float64(line.split()[-2])
                        while 1:
                            line = File.readline()
                            if line == '':
                                break
                            if line.find('Stop Module: ccsdt') !=-1:
                                break
                            if line.find('CCSD + T3') !=-1:
                                if ('CCSD(T)' not in self.energies):
                                    self.energies['CCSD(T)']={}
                                self.energies['CCSD(T)'][Field]= float64(line.split()[-1])

#----------------------------------------------------------------------------
# Gaussian routines
#----------------------------------------------------------------------------

class GaussianInputs(Inputs):
    """Gaussian 09 input routines"""

    def __init__(self, data, fstep):
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
        Inputs.__init__(self, data, fstep)

    def WriteInputs(self):

        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            xyz=''.join(xyz)
        except ValueError:
            print "Problem with *.xyz file?"
            sys.exit(1)

        fieldtail='\n'
        for i in range(10):
            fieldtail+='%7.4f\n' % 0.0

        # write input files
        for f in range(25):
            filename = self.data.replace('.xyz','')+'_F%2d_' % f
            filename = filename.replace(' ','0')
            ffield = fields(f,self.fstep)+fieldtail
            finput = self.tmpl.substitute(data=xyz, field=ffield, chk=filename+'.chk')
            open(filename+'.inp','w').write(finput)

class GAUSSIAN(Parser):
    """A gaussian 09 log parser."""

    def parsefiles(self):
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
                run(os.environ['g09root']+'/g09/formchk', chk)

            # parse current logfile
            self.parsefile(fchk)

    def parsefile(self,filename):
        """Read g09 energies and dipoles for this system."""

        File=open(filename)

        # Read control options
        mp2=0
        mp2dip=0
        mp3=0
        mp4=0
        ccsd=0
        ccsdt=0

        # available energies
        line=SkipLines(File,2)
        if line.find('MP2') !=-1:
            mp2=1
        if line.find('CCSD') !=-1:
            mp2=1
            mp2dip=0
            mp3=1
            mp4=1
            ccsd=1
            if line.find('CCSD(T)') !=-1:
                ccsdt=1
            else:
                ccsdt=0

        # available dipoles
        if mp2==1:
            if FindLine(File,'Total MP2 Density') !=-1:
                mp2dip=1
            else:
                mp2dip=0

        # setup dictionaries
        if 'SCF' not in self.energies:
            self.energies['SCF']={}
        if mp2==1 and 'MP2' not in self.energies:
            self.energies['MP2']={}
        if mp3==1 and 'MP3' not in self.energies:
            self.energies['MP3']={}
        if mp4==1 and 'MP4(SDQ)' not in self.energies:
            self.energies['MP4(SDQ)']={}
        if ccsd==1 and 'CCSD' not in self.energies:
            self.energies['CCSD']={}
        if ccsdt==1 and 'CCSD(T)' not in self.energies:
            self.energies['CCSD(T)']={}

        File.seek(0)

        # ... read SCF energy ...
        line = FindLine(File,'SCF Energy')
        ESCF = float64(line.split()[-1])

        # ... read MP2 energy ...
        if mp2==1:
            line = FindLine(File,'MP2 Energy')
            EMP2 = float64(line.split()[-1])

        # ... read MP3 energy ...
        if mp3==1:
            line = FindLine(File,'MP3 Energy')
            EMP3 = float64(line.split()[-1])

        # ... read MP4(SDQ) energy ...
        if mp4==1:
            line = FindLine(File,'MP4SDQ Energy')
            EMP4 = float64(line.split()[-1])

        # ... read CCSD energy ...
        if ccsd==1:
            line = FindLine(File,'Cluster Energy  ')
            ECCSD = float64(line.split()[-1])

        # ... read CCSD(T) energy ...
        if ccsdt==1:
            line = FindLine(File,'Cluster Energy with triples')
            ECCSDT = float64(line.split()[-1])

        # ... read applied external field ...
        line = FindLine(File,'External E-field')
        line = SkipLines(File,1).split()
        Field = ( round(float(line[1]),4), round(float(line[2]),4), round(float(line[3]),4) )

        # ... set base field ...
        if ( (self.fstep == 0 and abs(Field[0]) > 0) or
             (abs(Field[0]) > 0 and abs(Field[0]) < self.fstep) ):
            self.fstep = abs(Field[0])


        if Field not in self.energies['SCF']:
            self.energies['SCF'][Field]=ESCF
        if mp2==1 and Field not in self.energies['MP2']:
            self.energies['MP2'][Field]=EMP2
        if mp3==1 and Field not in self.energies['MP3']:
            self.energies['MP3'][Field]=EMP3
        if mp4==1 and Field not in self.energies['MP4(SDQ)']:
            self.energies['MP4(SDQ)'][Field]=EMP4
        if ccsd==1 and Field not in self.energies['CCSD']:
            self.energies['CCSD'][Field]=ECCSD
        if ccsdt==1 and Field not in self.energies['CCSD(T)']:
            self.energies['CCSD(T)'][Field]=ECCSDT

        # ... read dipoles ...
        if FindLine(File,'Dipole Moment') != -1:

            if mp2dip == 0:
                if 'SCF' not in self.dipoles:
                    self.dipoles['SCF']={}
            else:
                if 'MP2' not in self.dipoles:
                    self.dipoles['MP2']={}

            line = SkipLines(File,1)

            D = array(line.split(),dtype=float64)

            if mp2dip == 0 and Field not in self.dipoles['SCF']:
                self.dipoles['SCF'][Field]=D

            if mp2dip == 1 and Field not in self.dipoles['MP2']:
                self.dipoles['MP2'][Field]=D

#----------------------------------------------------------------------------
# Utilities
#----------------------------------------------------------------------------

def run(program, args):
    '''Run shell command'''
    os.system(program + ' ' + args)

def SkipLines(File,n):
    '''Read n lines from file f.'''

    for i in range(n):
        line = File.readline()
        if line == '' :
            break

    return line

def FindLine(File,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = File.readline()
        if line.find(pattern) !=-1 :
            break
        if line == '' :
            line = -1
            break

    return line

def ReFindLine(File,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = File.readline()
        if re.search(pattern,line) !=None : break
        if line == '' :
            linie = -1
            break

    return line

def PunchEn(out,i,F,Energies):
    '''Punch energies to out file'''
    row='%26.15f ' % Energies[F]
    row+='%7.4f' % F[i]
    row+='\t\tF=%s' % str(F)
    out.write( row+'\n' )

def ff25input(Energies, fstep):
    '''Write input for RZ's ff25.f'''
    out=open('ff25.inp','w')
    out.write('%f\n' % fstep)
    for i in range(25):
        field=fields(i,fstep).split()
        for j in range(len(field)):
            field[j]=float64(field[j])
        out.write('%25.16e _%d_\n' % (Energies[tuple(field)], i+1))

#----------------------------------------------------------------------------
# Property calculation routines
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
    mu_i=( (-2.0/3.0)*( E[Fi(i,f)] - E[Fi(i,-f)] ) + 
           ( E[Fi(i,2*f)] - E[Fi(i,-2*f)] )/12.0 )/abs(f)
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
    mu_i=( (2.0/3.0)*( D[Fi(i,f)][i] + D[Fi(i,-f)][i] ) - 
           (1.0/6.0)*( D[Fi(i,2*f)][i] + D[Fi(i,-2*f)][i] ) )
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

    Properties['B']=B
    Properties['Bx']=Bx
    Properties['By']=By
    Properties['Bz']=Bz

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
    log  = 'Static electric dipole properties (finite field '+Label+' dipole) \n\n'

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
    line = 3*U['b']['f'] + '\n\n'
    B = array([xxx,yxx,zxx,
               xyy,yyy,zyy,
               xzz,yzz,zzz],dtype=float64).reshape(3,3)
    log += '%15s %15s %15s\n' % ( 'xxx'.rjust(15), 'yxx'.rjust(15), 'zxx'.rjust(15) )
    log += line % ( xxx, yxx, zxx )
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

    Properties['B']=B
    Properties['Bx']=Bx
    Properties['By']=By
    Properties['Bz']=Bz

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
    log  = 'Static electric dipole properties (finite field '+Label+' energy) \n\n'

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
    line = 3*U['b']['f'] + '\n\n'
    B = array([xxx,yxx,zxx,
               xyy,yyy,zyy,
               xzz,yzz,zzz],dtype=float64).reshape(3,3)
    log += '%15s %15s %15s\n' % ( 'xxx'.rjust(15), 'yxx'.rjust(15), 'zxx'.rjust(15) )
    log += line % ( xxx, yxx, zxx )
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

    Units = {}

    if u.lower() == 'au':
        Units['u'] = 'a.u.'
        Units['m'] = {'x':1.0,          'f':'%15.4f ', 'r':4}
        Units['a'] = {'x':1.0,          'f':'%15.3f ', 'r':3}
        Units['b'] = {'x':1.0,          'f':'%15.2f ', 'r':2}
        Units['g'] = {'x':1.0,          'f':'%15.2f ', 'r':2}
    elif u.lower() == 'si':
        Units['u'] = 'si'
        Units['m'] = {'x':8.478358e-30, 'f':'%15.5e ', 'r':5} # C m
        Units['a'] = {'x':1.648778e-41, 'f':'%15.5e ', 'r':5} # C^2 m^2 J^-1
        Units['b'] = {'x':3.206361e-53, 'f':'%15.5e ', 'r':5} # C^3 m^3 J^-2
        Units['g'] = {'x':6.235377e-65, 'f':'%15.5e ', 'r':5} # C^4 m^4 J^-3
    elif u.lower() == 'asi':
        Units['u'] = 'asi'
        Units['m'] = {'x':8.4784e-30,   'f':'%15.5e ', 'r':5} # C m
        Units['a'] = {'x':1.8621e-30,   'f':'%15.5e ', 'r':5} # m^3
        Units['b'] = {'x':3.6213e-42,   'f':'%15.5e ', 'r':5} # m^4 V^-1
        Units['g'] = {'x':7.0423e-54,   'f':'%15.5e ', 'r':5} # m^5 V^-2
    elif u.lower() == 'esu':
        Units['u'] = 'esu'
        Units['m'] = {'x':2.5418e-18,   'f':'%15.5e ', 'r':5} # statvolt cm^2
        Units['a'] = {'x':1.4817e-25,   'f':'%15.5e ', 'r':5} # cm^3
        Units['b'] = {'x':8.6392e-33,   'f':'%15.5e ', 'r':5} # statvolt^-1 cm^4
        Units['g'] = {'x':5.0367e-40,   'f':'%15.5e ', 'r':5} # statvolt^-2 cm^5
    else:
        print "Unknown units requested!"
        sys.exit(1)

    # Selected tensor elements
    #               Mu:     Alpha:     Beta:         Gamma:
    # [[0, 1, 2],   x y z   xx xy xz   xxx xxy xxz   xxxx xxyy xxzz
    #  [3, 4, 5],           yx yy yz   yyx yyy yyz   yyxx yyyy yyzz
    #  [6, 7, 8]]           zx zy zz   zzx zzy zzz   zzxx zzyy zzzz
    # 
    PropertyIndex = {
        'Mu'    : [2, 'x'   ,'y'   ,'z'  ],
        'Alpha' : [8, 'xx'  ,'xy'  ,'xz'  , 'yx'  ,'yy'  ,'yz'  , 'zx'  ,'zy'  ,'zz'  ],
        'Beta'  : [8, 'xxx' ,'xxy' ,'xxz' , 'yyx' ,'yyy' ,'yyz' , 'zzx' ,'zzy' ,'zzz' ],
        'Gamma' : [8, 'xxxx','xxyy','xxzz', 'yyxx','yyyy','yyzz', 'zzxx','zzyy','zzzz'] }

    # Descriptions of properties
    PropertyDescription = {
        'Mu'    : '# Dipole Moment Vector i=' + \
                     PropertyIndex['Mu'][PropertyIndex['Mu'][0]+1],
        '|D|'   : '# Total Dipole Moment',
        'Alpha' : '# Polarizability Tensor i=' + \
                     PropertyIndex['Alpha'][PropertyIndex['Alpha'][0]+1],
        '<A>'   : '# Isotropic Polatizability',
        '<B>'   : '# Anisotropy of Polarizability (Z-axis is the rotation axis)',
        'Beta'  : '# First hyperpolarizability tensor i=' + \
                     PropertyIndex['Beta'][PropertyIndex['Beta'][0]+1],
        'B(Z)'  : '# Vector component of hyperpolarizability tensor (Z is the permament dipole moment direction)',
        'Gamma' : '# Second hyperpolarizability tensor i=' + \
                     PropertyIndex['Gamma'][PropertyIndex['Gamma'][0]+1],
        '<G>'   : '# Scalar component of second hyperpolarizability tensor given by the isotropic average' }

    return Units

def periodic(mendeleiev):
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

def atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''

    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)

#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__": Main(sys.argv[1:])

