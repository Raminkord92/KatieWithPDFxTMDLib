#!/usr/bin/env python

import os,sys,subprocess,platform
from os.path import join as jpth

AVHLIBpath = ''
LHAPDFpath = ''
TMDLIBpath = ''
TMDLIBmulti = ''
PDFXTMDpath = ''
GSLpath = ''
MINCASpath = ''
MERGERpath = ''
FC = ''
LINKMETHOD = 'dynamic'
if platform.system()=='Darwin':
    LINKER = 'gcc'
else:
    LINKER = 'ld'
from settings import *

subprocess.call(['python',jpth(AVHLIBpath,'config.py')])

KATIEpath = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(jpth(AVHLIBpath,'python'))
from avh import *

if not os.path.exists(jpth(KATIEpath,'build')): pfd.mkdirp(jpth(KATIEpath,'build'))

edpy.set_value( 'KATIEpath'  ,KATIEpath  ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'AVHLIBpath' ,AVHLIBpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'LHAPDFpath' ,LHAPDFpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'TMDLIBpath' ,TMDLIBpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'TMDLIBmulti',TMDLIBmulti,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'PDFXTMDpath' ,PDFXTMDpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'GSLpath'    ,GSLpath    ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'MINCASpath' ,MINCASpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'MERGERpath' ,MERGERpath ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'FC'         ,FC         ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'LINKER'     ,LINKER     ,jpth(KATIEpath,'create_def.py') )
edpy.set_value( 'LINKMETHOD' ,LINKMETHOD ,jpth(KATIEpath,'create_def.py') )
ed.set_value( ['^ *','','=','$'] ,'KATIEpath' ,KATIEpath ,jpth(KATIEpath,'run.sh') )

if len(sys.argv)>1:
    if sys.argv[1]=='lib':
        pfd.rmdirstar(jpth(KATIEpath,'build'))
        cr.execute(['python',jpth(KATIEpath,'create_lib.py'),'lib'])
