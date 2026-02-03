#!/usr/bin/env python
import platform
from create_def import *

def createSource(lines):
    cr.addavh(AVHLIBpath,'iounits',[],lines)
    cr.addavh(AVHLIBpath,'kindpars',[],lines)
    cr.addavh(AVHLIBpath,'prnt',[],lines)
    cr.addavh(AVHLIBpath,'mathcnst',[],lines)
    cr.addavh(AVHLIBpath,'iotools',[],lines)
    cr.addavh(AVHLIBpath,'bivinp',[],lines)
    cr.addavh(AVHLIBpath,'doloops',[],lines)
    cr.addavh(AVHLIBpath,'statdist',[],lines)
    cr.addavh(AVHLIBpath,'random',[],lines)
    cr.addavh(AVHLIBpath,'ranvar',[],lines)
    cr.addavh(AVHLIBpath,'grid',[],lines)
    cr.addavh(AVHLIBpath,'mctools',[],lines)
    cr.addavh(AVHLIBpath,'lorentz',[],lines)
    cr.addavh(AVHLIBpath,'kskeleton',[],lines)
    cr.addavh(AVHLIBpath,'skeletools',[],lines)
    cr.addavh(AVHLIBpath,'kinematics',[],lines)
    cr.addavh(AVHLIBpath,'colorflow',[],lines)
    cr.addavh(AVHLIBpath,'qftmodels',[],lines)
    cr.addavh(AVHLIBpath,'kaleu',[],lines)
    cr.addfile(jpth(srcDir,'version.f90'),lines)
    cr.addfile(jpth(srcDir,'model_kaleu_qcd.f90'),lines)
    cr.addfile(jpth(srcDir,'ranInst.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_DISinst.f90'),lines)
    cr.addfile(jpth(srcDir,'events.f90'),lines)
    cr.addfile(jpth(srcDir,'partlumi.f90'),lines)
    cr.addfile(jpth(srcDir,'pdfs.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_quadpack.f90'),lines)
    cr.addfile(jpth(srcDir,'sudakov.f90'),lines)
    cr.addfile(jpth(srcDir,'histogramtools.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_eventfile.f90'),lines)
    cr.addfile(jpth(srcDir,'read_events.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_amplitudes.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_matrixelement.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_itmds.f90'),lines)
    fpp.incl(srcDir,lines)
    fpp.xpnd('path_tbldir',buildDir+os.sep,lines)
    fpp.xpnd('rngenerator','call rangen(#1)',lines)
    fpp.xpnd('rngmodule','avh_random',lines)
    fpp.xpnd('realknd2','real(kind(1d0))',lines)
    fpp.xpnd('realknd1','real(kind(1e0))',lines)
    fpp.xpnd('complex2','complex(kind(1d0))',lines)
    fpp.xpnd('integer1','integer(selected_int_kind(4))',lines)
    fpp.xpnd('integer0','integer(selected_int_kind(2))',lines)
    fpp.blck('sudakovTableDir',frmt_char(buildDir+os.sep,cont=False),lines)
    fpp.blck('itmdColorTableFile',frmt_char(jpth(srcDir,'katie_ITMDmatrix.tbl'),cont=False),lines)
    if TMDLIBpath=='':
        TMDLIBversion = '0'
    elif TMDLIBmulti!='yes':
        TMDLIBversion = '1'
    else:
        TMDLIBversion = '2'
    fpp.case('TMDlibVersion',TMDLIBversion,lines)
    cr.putdate(lines)
    return lines


if len(sys.argv)<2:
    cr.prnt('================================================================')
    cr.prnt(' Examples:')
    cr.prnt(' $ ./create.py clean')
    cr.prnt(' $ ./create.py lib')
    cr.prnt(' $ ./create.py input=<filename> dir=<dirname>')
    cr.prnt(' $ ./create.py compile sourcefile.f90')
    cr.prnt('================================================================')

elif sys.argv[1]=='clean' or sys.argv[1]=='task=clean':
    pfd.rmdirstar(buildDir)

elif sys.argv[1]=='source' or sys.argv[1]=='task=source':
    lines = []
    createSource(lines)
    cr.wfile(jpth(buildDir,'head.f90'),lines)

elif sys.argv[1]=='head' or sys.argv[1]=='lib' or sys.argv[1]=='task=lib':
    lines = []
    createSource(lines)

    options = []
    delpat = KATIEpath+'/'
    if len(sys.argv)>2:
        if sys.argv[2]=='twostep':
            options.append('onlyPrint')
            delpat=''

    cr.modlib('head',FCF+['-I '+srcDir],lines,buildDir \
                    ,linkmethod=LINKMETHOD,linker=LINKER \
                    ,delPattern=delpat,options=options)

    if MINCASpath!='' \
       and not os.path.exists(jpth(buildDir,'libmincas.dylib'))\
       and not os.path.exists(jpth(buildDir,'libmincas.so')):
        CCF = ['g++','-std=c++17','-lstdc++','-fPIC','-c','-I'+jpth(MINCASpath,'include'),'-I'+MERGERpath]
        objects = [jpth(buildDir,'mincas0.o'),jpth(buildDir,'mincas1.o')\
                  ,jpth(buildDir,'mincas2.o'),jpth(buildDir,'mincas3.o')\
                  ,jpth(buildDir,'mincas4.o')]
        cr.execute(CCF+['-o',objects[0],jpth(MINCASpath,'src','mixmax.c')],delPattern='')
        cr.execute(CCF+['-o',objects[1],jpth(MINCASpath,'src','mt19937ar.c')],delPattern='')
        cr.execute(CCF+['-o',objects[2],jpth(MINCASpath,'src','mincas.c')],delPattern='')
        cr.execute(CCF+['-o',objects[3],jpth(MERGERpath,'rlib.cpp')],delPattern='')
        cr.execute(CCF+['-o',objects[4],jpth(MERGERpath,'launchmincas.cpp')],delPattern='')
        if platform.system()=='Darwin':
            cr.execute([LINKER,'-dynamiclib','-o',jpth(buildDir,'libmincas.dylib')]+objects,delPattern='')
        else:
            cr.execute([LINKER,'-shared','-o',jpth(buildDir,'libmincas.so')]+objects,delPattern='')
    
    compile_file(jpth(srcDir,'merge_raw.f90'),target=jpth(buildDir,'merge_raw'),options=options)

elif sys.argv[1]=='katamp':
    if len(sys.argv)<3: targetFile=jpth(buildDir,'katamp.f90')
    else: targetFile= sys.argv[2]
    lines = []
    cr.addavh(AVHLIBpath,'iounits',[],lines)
    cr.addavh(AVHLIBpath,'kindpars',[],lines)
    cr.addavh(AVHLIBpath,'prnt',[],lines)
    cr.addavh(AVHLIBpath,'mathcnst',[],lines)
    cr.addavh(AVHLIBpath,'doloops',[],lines)
    cr.addavh(AVHLIBpath,'lorentz',[],lines)
    cr.addavh(AVHLIBpath,'kskeleton',[],lines)
    cr.addavh(AVHLIBpath,'skeletools',[],lines)
    cr.addavh(AVHLIBpath,'kinematics',[],lines)
    cr.addavh(AVHLIBpath,'colorflow',[],lines)
    cr.addavh(AVHLIBpath,'qftmodels',[],lines)
    cr.addfile(jpth(srcDir,'version.f90'),lines)
    cr.addfile(jpth(srcDir,'katie_amplitudes.f90'),lines)
    fpp.incl(srcDir,lines)
    fpp.xpnd('path_tbldir',buildDir+os.sep,lines)
    fpp.xpnd('realknd2','real(kind(1d0))',lines)
    fpp.xpnd('realknd1','real(kind(1e0))',lines)
    fpp.xpnd('complex2','complex(kind(1d0))',lines)
    fpp.xpnd('integer1','integer(selected_int_kind(4))',lines)
    fpp.xpnd('integer0','integer(selected_int_kind(2))',lines)
    fpp.blck('itmdColorTableFile',frmt_char(jpth(srcDir,'katie_ITMDmatrix.tbl'),cont=False),lines)
    cr.putdate(lines)
    cr.wfile(targetFile,lines,delPattern='')
    
elif sys.argv[1]=='compile' or sys.argv[1]=='task=compile':
    if len(sys.argv)==3: compile_file(sys.argv[2])
    else               : compile_list(sys.argv[2:])

elif sys.argv[1]=='help':
    if len(sys.argv)<3:
        print('Please include a help topic, for example: help compile')
    elif sys.argv[2]=='compile':
        cmndList = cmplcmd_libs('DUMMYsrc','DUMMYtarget')
        print('To compile source requiring KaTie, include the following lines:')
        for item in cmndList:
            if item[0][:2]=='-I' or item[0][:2]=='-L':
                print(' ',' '.join(item))
    else:
        print('No help available for topic: ',sys.argv[2])

else:
    cr.prnt('ERROR in create_lib: argument(s) not defined.')
