#!/usr/bin/env python
import os,sys,re,subprocess
from os.path import join as jpth

rpathFlag = '-Wl,-rpath,'

AVHLIBpath = '/root/physics/montecarlo/avhlib'
KATIEpath = '/root/physics/katieWithPDFxTMD'
KATIEpath = '/root/physics/katieWithPDFxTMD'
srcDir = jpth(KATIEpath,'src')
buildDir = jpth(KATIEpath,'build')

sys.path.append(jpth(AVHLIBpath,'python'))
from avh import *
import cuts_and_pars
from cuts_and_pars import pdf_particle as pdf_particle
from cuts_and_pars import particle as particle
from cuts_and_pars import frmt as frmt
from cuts_and_pars import frmt_char as frmt_char
from cuts_and_pars import frmt_list as frmt_list
from cuts_and_pars import Sdouble as Sdouble

LINKMETHOD = 'dynamic'
LINKER = 'ld'
FC = 'gfortran -fcheck=bounds'
LHAPDFpath = '/usr/local/lib'
TMDLIBpath = '/usr/local/lib'
TMDLIBmulti = ''
PDFXTMDpath = '/usr/local/lib'
GSLpath = '/usr/local/lib'
MINCASpath = ''
MERGERpath = ''

FCF = FC.split()

def print_array(array):
    rslt = '['+array[0]
    for jj in range(1,len(array)): rslt = rslt+','+array[jj]
    rslt = rslt+']'
    return rslt

def print_process(array):
    rslt = '['+particle(array[0])+','+particle(array[1])+']'+',['+particle(array[3])
    for jj in range(4,len(array)): rslt = rslt+','+particle(array[jj])
    rslt = rslt+']'
    return rslt

def print_inst(array):
    rslt = '['+particle(array[0])+','+particle(array[1])+']'
    return rslt


def cmplcmd_libs(srcList,target,options=[]):
    outlist = []
    if ('withMINCAS' in options) and (not '-fno-underscoring' in FCF): 
        FCF.append('-fno-underscoring')
    outlist.append(FCF)
    outlist.append(['-I',buildDir])
    outlist.append(srcList)
    outlist.append(['-L'+buildDir,rpathFlag+buildDir,'-lhead'])
    if 'withMINCAS' in options: outlist.append(['-lmincas'])
    outlist.append(['-L'+LHAPDFpath,rpathFlag+LHAPDFpath,'-lLHAPDF'])
    if TMDLIBpath!='':
        if os.path.isfile(jpth(TMDLIBpath,'TMDlib-config')):
            line = subprocess.check_output([jpth(TMDLIBpath,'TMDlib-config') \
                                           ,'--ldflags']).decode(sys.stdout.encoding)
            words = line.rstrip().split(' ')
            outlist.append([words[0],rpathFlag+words[0][2:]])
            outlist.append(words[1:])
        else:
            outlist.append(['-L'+TMDLIBpath,rpathFlag+TMDLIBpath,'-lTMDlib'])
        outlist.append(['-L'+GSLpath,rpathFlag+GSLpath,'-lgsl','-lgslcblas','-lm'])
    if PDFXTMDpath!='':
        outlist.append(['-L'+PDFXTMDpath,rpathFlag+PDFXTMDpath,'-lPDFxTMDLib','-lstdc++'])
    outlist.append([cr.moduleDir(FCF[0]),buildDir])
    outlist.append(['-o',target+'.out'])
    return outlist

def cmplcmd_file(srcFile,target,options=[]):
    mainName,extension = os.path.splitext(srcFile)
    if extension=='.f' or extension=='.f90' or extension=='f03': 
        srcfile = srcFile
    else:
        srcfile = mainName+'.f90'
    if target=='':
        return cmplcmd_libs([srcfile],mainName,options)
    else:
        return cmplcmd_libs([srcfile],target,options)

def cmplcmd_list(srcList,target,options=[]):
    mainName,extension = os.path.splitext(srcList[-1])
    if target=='':
        return cmplcmd_libs(srcList,mainName,options)
    else:
        return cmplcmd_libs(srcList,target,options)

def flatlist(inlist):
    outlist = []
    for sublist in inlist:
        for item in sublist:
            outlist.append(item)
    return outlist

def cmplcmd_file_as_str(srcFile,target='',options=[]):
    return ' '.join(flatlist(cmplcmd_file(srcFile,target,options)))

def cmplcmd_file_as_strlist(srcFile,target='',options=[]):
    outlist = []
    for item in cmplcmd_file(srcFile,target,options):
        outlist.append(' '.join(item)+' \\\n')
    return outlist

def compile_file(srcFile,target='',options=[]):
    cr.execute(flatlist(cmplcmd_file(srcFile,target,options)),delPattern='',options=options)

def compile_list(srcList,target='',options=[]):
    cr.execute(flatlist(cmplcmd_list(srcList,target,options)),delPattern='',options=options)


def myFloat(xx):
    return float(re.sub('D','E',re.sub('d','e',xx)))

def tmdExists(tmpPath,tmpFile):
    if not os.path.exists(jpth(tmpPath[0],tmpFile[1])):
        cr.prnt('')
        cr.prnt('ERROR: the file')
        cr.prnt('           '+jpth(tmpPath[0],tmpFile[1]))
        cr.prnt('       does not seem to exist. KaTie does not provide TMDs.')
        cr.prnt('       It needs rectangular grids with one of the two formats:')
        cr.prnt('           ln(x) ln(kT^2) x*F(x,kT)')
        cr.prnt('           ln(x) ln(kT^2) ln(mu^2) x*F(x,kT,mu)')
        cr.prnt('       You can find such grids, for example the set MRW-CT10nlo, at')
        cr.prnt('           https://tmdlib.hepforge.org/downloads?f=datafiles/')
        cr.prnt('       Alternatively, you can use TMDlib itself to evaluate TMDs,')
        cr.prnt('       by setting,')
        cr.prnt('           TMDlibSet = MRW-CT10nlo')
        cr.prnt('       in the input file.')
        cr.prnt('')
        sys.exit()
