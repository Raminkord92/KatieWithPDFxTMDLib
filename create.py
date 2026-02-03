#!/usr/bin/env python
import os,sys,re
from os.path import join as jpth

rpathFlag = '-Wl,-rpath,'

AVHLIBpath = '/home/user0/repos/avhlib'
KATIEpath = '/home/user0/repos/katie'
srcDir = jpth(KATIEpath,'src')
buildDir = jpth(KATIEpath,'build')

sys.path.append(jpth(AVHLIBpath,'python'))
from avh import *
import cuts_and_pars
from cuts_and_pars import particle as particle
from cuts_and_pars import Sdouble as Sdouble

FC = 'gfortran -fcheck=bounds'
LHAPDFpath = '/usr/lib'
TMDLIBpath = '/usr/local/lib'
PDFXTMDpath = '/usr/local/lib'
GSLpath = '/usr/lib'

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

def compile_main_basic(FCFIsrc,mainName):
    args2 = ['-L'+LHAPDFpath,rpathFlag+LHAPDFpath,'-lLHAPDF']
    if TMDLIBpath!='':
        args2 = args2 + ['-L'+GSLpath,rpathFlag+GSLpath,'-lgsl','-lgslcblas','-lm']
        args2 = args2 + ['-L'+TMDLIBpath,rpathFlag+TMDLIBpath,'-lTMDlib']
    args2 = args2 + ['-L'+buildDir]
    args3 = ['-lhead',cr.moduleDir(FCF[0]),buildDir,'-o',mainName+'.out']
    return FCFIsrc+args2+args3

def compile_main_arglist(name,label=''):
    mainName,extension = os.path.splitext(name)
    if extension=='.f' or extension=='.f90' or extension=='f03': 
        srcfile = name
    else:
        srcfile = mainName+'.f90'
    return compile_main_basic(FCF+['-I',buildDir,srcfile],mainName+label)

def compile_main_list_arglist(srcList,label=''):
    mainName,extension = os.path.splitext(srcList[-1])
    return compile_main_basic(FCF+['-I',buildDir]+srcList,mainName+label)

def compile_main(name,label=''):
    cr.execute(compile_main_arglist(name,label),delPattern='')

def compile_main_list(srcList,label=''):
    cr.execute(compile_main_list_arglist(srcList,label),delPattern='')

def frmt_char(string,linelen=84,cont=True):
    stringlen = len(string)
    words = []
    offset = 0
    while offset<stringlen:
        words.append('//\''+string[offset:min(offset+linelen-5,stringlen)]+'\'&\n')
        offset = offset+linelen-5
    words[0] = re.sub(r'^//','  ',words[0])
    if cont==False: words[-1] = re.sub(r'&','',words[-1])
    return words

def myFloat(xx):
    return float(re.sub('D','E',re.sub('d','e',xx)))

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

elif sys.argv[1]=='head' or sys.argv[1]=='lib' or sys.argv[1]=='task=lib':
    lines = []
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
    cr.addfile(jpth(srcDir,'events.f90'),lines)
    cr.addfile(jpth(srcDir,'read_events.f90'),lines)
    cr.addfile(jpth(srcDir,'partlumi.f90'),lines)
    cr.addfile(jpth(srcDir,'pdfs.f90'),lines)
    cr.addfile(jpth(srcDir,'amplitude.f90'),lines)
    fpp.xpnd('path_tbldir',buildDir+os.sep,lines)
    fpp.xpnd('rngenerator','call rangen(#1)',lines)
    fpp.xpnd('rngmodule','avh_random',lines)
    fpp.xpnd('realknd2','real(kind(1d0))',lines)
    fpp.xpnd('realknd1','real(kind(1e0))',lines)
    fpp.xpnd('complex2','complex(kind(1d0))',lines)
    fpp.xpnd('integer1','integer(selected_int_kind(4))',lines)
    fpp.xpnd('integer0','integer(selected_int_kind(2))',lines)
    if TMDLIBpath!='':
        fpp.case('withTMDlib','yes',lines)
    else:
        fpp.case('withTMDlib','no',lines)
    if PDFXTMDpath!='':
        fpp.case('withPDFxTMD','yes',lines)
    else:
        fpp.case('withPDFxTMD','no',lines)
    cr.putdate(lines)
    cr.modlib('head',FCF+['-I '+srcDir],lines,buildDir,delPattern=KATIEpath+'/')
    cr.execute(FCF+[jpth(srcDir,'merge_raw.f90') \
                   ,'-o',jpth(buildDir,'merge_raw.out')],delPattern=KATIEpath+'/')
    
elif sys.argv[1]=='compile' or sys.argv[1]=='task=compile':
    if len(sys.argv)==3: compile_main(sys.argv[2])
    else               : compile_main_list(sys.argv[2:])

else:
    for word in sys.argv:
        if re.match('^input=',word):
            inputFile = os.path.abspath(re.sub('^input=','',word))
        if re.match('^dir=',word):
            targetDir = os.path.abspath(re.sub('^dir=','',word))

    inputLines = ed.lines_from_file(inputFile)
    ii = 0
    for line in inputLines:
        if cr.lhs('include file','=',line):
            insertFile = cr.rhs('=',line)
            del inputLines[ii]
            inputLines[ii:ii] = ed.lines_from_file(insertFile)
            ii = ii-1
        ii = ii+1

    offshell = ['0','0']
    Ngroup = '1'
    Nfinst = ['2']
    Nflavors = '0'
    lhaSet = 'MSTW2008lo68'
    lhaSetSet = False
    withLHAPDF = 'yes'
    TMDlibSet = []
    PDFxTMDSet = []
    PDFxTMDType = 'tmd'
    PDFxTMDCoupling = 'no'
    PDFxTMDCouplingSet = ''
    tmdTableDir = '/home/user0/kTfac/krzysztof/'
    Qsqueez = '-1d0'
    seed = '987654321'
    helicity = 'sampling'
    optimHel = 'sampling'
    Noptim = '100000' #
    NbatchGrid = '0'  # Noptim is overruled by NbatchGrid,
    Nbatch = '0'      # NbatchGrid is overruled by Nbatch.
    NstepGrid = '100'
    Nstep = '10'
    thrs = '1d-3'
    instPartons = 'g_g'
    instLabel = 'g_g'
    sigma_eff = '1d0'
    partlumi = 'individual'
    fluxFactor = 'collinear'
    Esoft = '20'
    Ebeam1 = '0'
    Ebeam2 = '0'
    pars = cuts_and_pars.parameters()

    for line in inputLines:
        if cr.lhs(        'Ngroup','=',line): Ngroup = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(        'Nfinst','=',line): Nfinst = cr.rhs('=',line).split()
        if cr.lhs(          'seed','=',line): seed = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(        'Noptim','=',line): Noptim = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(        'Nbatch','=',line): Nbatch = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(         'Nstep','=',line): Nstep = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(          'thrs','=',line): thrs = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(    'NbatchGrid','=',line): NbatchGrid = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(     'NstepGrid','=',line): NstepGrid = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(     'lha[sS]et','=',line):
            lhaSet = cr.rhs('=',line).split(' ',1)[0]
            lhaSetSet = True
        if cr.lhs('useLHAPDF','=',line) or cr.lhs('withLHAPDF','=',line):
            rhs = cr.rhs('=',line).split(' ',1)[0].lower()
            if rhs in ('yes','y','1','true'):
                withLHAPDF = 'yes'
            elif rhs in ('no','n','0','false'):
                withLHAPDF = 'no'
        if cr.lhs(     'TMDlibSet','=',line): TMDlibSet = cr.rhs('=',line).split()
        if cr.lhs(   'PDFxTMDSet','=',line): PDFxTMDSet = cr.rhs('=',line).split()
        if cr.lhs( 'PDFxTMD[Tt]ype','=',line) or cr.lhs( 'PDFxTMD[Mm]ode','=',line):
            PDFxTMDType = cr.rhs('=',line).split(' ',1)[0].lower()
        if cr.lhs('PDFxTMDCouplingSet','=',line) or cr.lhs('PDFxTMDAlphaSet','=',line):
            PDFxTMDCouplingSet = cr.rhs('=',line).split(' ',1)[0]
            PDFxTMDCoupling = 'yes'
        if cr.lhs('PDFxTMDCoupling','=',line) or cr.lhs('PDFxTMDAlpha','=',line):
            rhs = cr.rhs('=',line).split(' ',1)[0]
            if rhs.lower() in ('yes','y','1','true'):
                PDFxTMDCoupling = 'yes'
            elif rhs.lower() in ('no','n','0','false'):
                PDFxTMDCoupling = 'no'
            else:
                PDFxTMDCouplingSet = rhs
                PDFxTMDCoupling = 'yes'
        if cr.lhs(       'Qsqueez','=',line): Qsqueez = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs(      'Nflavors','=',line): Nflavors = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(      'offshell','=',line): offshell = cr.rhs('=',line).split()
        if cr.lhs(     'sigma_eff','=',line): sigma_eff = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs(         'Esoft','=',line): Esoft = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs(        'Ebeam1','=',line): Ebeam1 = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs(        'Ebeam2','=',line): Ebeam2 = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs(      'partlumi','=',line): partlumi = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs('flux *[fF]actor','=',line): fluxFactor = cr.rhs('=',line).split(' ',)[0]
        if cr.lhs('print skeleton','=',line): printskel = cr.rhs('=',line).split(' ',1)[0]
        if cr.lhs(      'helicity','=',line): 
            helicity = cr.rhs('=',line).split(' ',1)[0]
            optimHel = helicity
        if cr.lhs('optimization helicity','=',line):
            optimHel = cr.rhs('=',line).split(' ',1)[0]

        if cr.lhs(    'mass','=',line): pars.mass(    cr.rhs('=',line))
        if cr.lhs(  'switch','=',line): pars.switch(  cr.rhs('=',line))
        if cr.lhs('coupling','=',line): pars.coupling(cr.rhs('=',line))
        if cr.lhs(     'Ecm','=',line): pars.energy( cr.rhs('=',line))

    pars.esoft(Esoft)
    if Ebeam1=='0': pars.eNegRap('Ecm/2')
    else:           pars.eNegRap(Sdouble(Ebeam1))
    if Ebeam2=='0': pars.ePosRap('Ecm/2')
    else:           pars.ePosRap(Sdouble(Ebeam2))
    cuts = []
    for ii in range(0,1+int(Ngroup)): cuts.append(cuts_and_pars.cuts())
    for line in inputLines:
        if cr.lhs('cut','=',line):
             rhsline = cr.rhs('=',line)
             rhs = rhsline.split()
             if rhs[0]=='group':
                 cuts[int(rhs[1])].add_cut(re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
             else:
                 cuts[0].add_cut(rhsline)
        if cr.lhs('scale','=',line): 
             rhsline = cr.rhs('=',line)
             rhs = rhsline.split()
             if rhs[0]=='group':
                 cuts[int(rhs[1])].set_scale(re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
             elif rhs[0]=='entry':
                 cuts[0].set_scale(re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline),group=rhs[1])
             else:
                 cuts[0].set_scale(rhsline)

    Noptim = re.sub(',','',Noptim)
    NbatchGrid = re.sub(',','',NbatchGrid)
    Nbatch = re.sub(',','',Nbatch)
    if NbatchGrid!='0': Noptim=str(int(int(NbatchGrid)*int(NstepGrid)))
    if Nbatch!='0': Noptim=str(int(int(Nbatch)*int(Nstep)))

    procNum = []
    procDir = []
    procList = []
    instPartons = []
    ii = 0
    for line in inputLines:
        if cr.lhs('process','=',line):
            procList.append(cr.rhs('=',line))
            procNum.append(str(ii+1).zfill(3))
            procDir.append('proc'+procNum[ii])
            word = procList[-1].split()
            instPartons.append(re.sub('~','b',word[0]+'_'+word[1]))
            if partlumi!='combined':
                if word[0][0]=='q' or word[0][0]=='r' \
                or word[1][0]=='q' or word[1][0]=='r':
                    cr.prnt('')
                    cr.prnt('WARNING from '+KATIEpath+'/create.py:')
                    cr.prnt('  if you indicate quarks with q/r, then you should actually set')
                    cr.prnt('  partlumi = combined')
                    cr.prnt('  (Yes, the default used to be different, sorry.)')
                    cr.prnt('  Setting  partlumi = combined  myself...')
                    cr.prnt('')
                    partlumi = 'combined'
            if instPartons[-1]=='g_qb':  instPartons[-1]='g_q' 
            if instPartons[-1]=='qb_g':  instPartons[-1]='q_g' 
            if instPartons[-1]=='qb_qb': instPartons[-1]='q_q' 
            if instPartons[-1]=='qb_q':  instPartons[-1]='q_qb' 
            if instPartons[-1]=='r_q':   instPartons[-1]='q_r' 
            if instPartons[-1]=='q_rb':  instPartons[-1]='q_r' 
            if instPartons[-1]=='rb_q':  instPartons[-1]='q_r' 
            if instPartons[-1]=='qb_r':  instPartons[-1]='q_r' 
            if instPartons[-1]=='r_qb':  instPartons[-1]='q_r' 
            if instPartons[-1]=='qb_rb': instPartons[-1]='q_r' 
            if instPartons[-1]=='rb_qb': instPartons[-1]='q_r' 
            if instPartons[-1]=='A_A': partlumi='photons' 
            ii = ii+1

    withTMDlib = 'no'
    ranInstMesh = 'no'
    tmdLines = []
    withPDFxTMD = 'no'
    withPDFxTMDTMD = 'no'
    withPDFxTMDCPDFA = 'no'
    withPDFxTMDCPDFB = 'no'
    withPDFxTMDCoupling = PDFxTMDCoupling
    if len(PDFxTMDSet)>0:
        withPDFxTMD = 'yes'
        withPDFxTMDTMD = 'yes'
        pdfTypeA = 'PDFxTMD0'
        pdfTypeB = 'PDFxTMD0'
        if len(PDFxTMDSet)<2: PDFxTMDSet = PDFxTMDSet + ['2212','2212']
        if len(PDFxTMDSet)<3: PDFxTMDSet = PDFxTMDSet + [PDFxTMDSet[1]]
        if len(PDFxTMDSet)<4: PDFxTMDSet = PDFxTMDSet + ['0']
        PDFxTMDSetB = PDFxTMDSet
        if withPDFxTMDCoupling=='yes' and PDFxTMDCouplingSet=='':
            PDFxTMDCouplingSet = PDFxTMDSet[0]
    elif len(TMDlibSet)==0:
        if myFloat(Qsqueez)<0.0:
            pdfTypeA = 'grid'
            pdfTypeB = 'grid'
        else:
            pdfTypeA = 'squeez'
            pdfTypeB = 'squeez'
        for line in inputLines:
            if cr.lhs('tmdPath','=',line) or cr.lhs('tmdpath','=',line) or cr.lhs('tmdTableDir','=',line):
                tmpPath = cr.rhs('=',line).split()
                tmdLines.append('  call set_tmdpath(&\n')
                tmdLines = tmdLines + frmt_char(tmpPath[0],cont=True)
                tmdLines.append('  )\n')
            if cr.lhs('tmdPdf','=',line) or cr.lhs('tmdpdf','=',line):
                tmpFile = cr.rhs('=',line).split()
                tmdLines.append('  call add_tmdpdf('+particle(tmpFile[0])+',\''+tmpFile[1]+'\')\n')
                if not os.path.exists(jpth(tmpPath[0],tmpFile[1])):
                    cr.prnt('')
                    cr.prnt('ERROR: the file '+jpth(tmpPath[0],tmpFile[1]))
                    cr.prnt('       does not seem to exist.')
                    cr.prnt('       KaTie does not provide TMDs.')
                    cr.prnt('       It needs rectangular grids with one of the two formats:')
                    cr.prnt('         ln(x) ln(kT^2) TMDvalue')
                    cr.prnt('         ln(x) ln(kT^2) ln(mu^2) TMDvalue')
                    cr.prnt('       Alternatively, you can use TMDlib.')
                    cr.prnt('')
                    sys.exit()
    else:
        withTMDlib = 'yes'
        pdfTypeA = 'TMDlib'
        pdfTypeB = 'TMDlib'
        if len(TMDlibSet)<2: TMDlibSet = TMDlibSet + ['2212','2212']
        if len(TMDlibSet)<3: TMDlibSet = TMDlibSet + [TMDlibSet[1]]
    if withLHAPDF=='yes' and PDFxTMDType=='cpdf' and withPDFxTMDCoupling=='yes' and not lhaSetSet:
        withLHAPDF = 'no'
    if PDFxTMDType=='cpdf':
        withPDFxTMDTMD = 'no'
    if offshell[1]=='0':
        if withPDFxTMD=='yes' and PDFxTMDType=='cpdf':
            pdfTypeA = 'PDFxTMDCPDF'
            withPDFxTMDCPDFA = 'yes'
        else:
            pdfTypeA = 'LHAPDF'
    if offshell[0]=='0':
        if withPDFxTMD=='yes' and PDFxTMDType=='cpdf':
            pdfTypeB = 'PDFxTMDCPDF'
            withPDFxTMDCPDFB = 'yes'
        else:
            pdfTypeB = 'LHAPDF'

    exec_all = ['Nproc='+str(len(procList))+'\n','prefix='+targetDir+'\n']
    recompile_all = ['Nmain='+str(1+len(procList))+'\n' \
                    ,'array[  0]="'+' '.join(compile_main_arglist(jpth(targetDir,'main.f90')))+'"\n']

    mainProcLines = []
    for ii in range(0,len(procList)):
        pfd.mkdir(jpth(targetDir,procDir[ii]))

        flavor = re.sub(r' *factor.*','',procList[ii]).split()
        factor = re.sub(r'Nf','Nflavors',re.sub(r' *groups.*','',cr.rhs('=',procList[ii])))
        groups = re.sub(r' *pNonQCD.*','',cr.rhs('=',cr.rhs('=',procList[ii]))).split()
        pNonQCD = cr.rhs('=',cr.rhs('=',cr.rhs('=',procList[ii]))).split()
        
        mainProcLines.append('  proc('+procNum[ii]+')%group(1:'+str(len(groups))+') = '+print_array(groups)+'\n')
        mainProcLines.append('  proc('+procNum[ii]+')%pNonQCD = '+print_array(pNonQCD)+'\n')
        mainProcLines.append('  proc('+procNum[ii]+')%Nfinst = '+str(len(flavor)-3)+'\n')
        mainProcLines.append('  call proc('+procNum[ii]+')%flav%put('+print_process(flavor)+',get_anti)\n')
        mainProcLines.append('  proc('+procNum[ii]+')%factor = '+factor+'\n')
        mainProcLines.append('  proc('+procNum[ii]+')%label = &\n')
        mainProcLines = mainProcLines + frmt_char(jpth(targetDir,procDir[ii],'MC_'),cont=False)
        if partlumi=='combined':
          mainProcLines.append('  proc('+procNum[ii]+')%pdf => pdf_'+instPartons[ii]+'\n')
        else:
          mainProcLines.append('  proc('+procNum[ii]+')%parton(-2:-1) = '+print_inst(flavor)+'\n')

        optimProcLines = []
        optimProcLines.append('  call flavor%put('+print_process(flavor)+',get_anti)\n')
        optimProcLines.append('  ampFac = '+factor+'\n')

        srcLines = []
        cr.addfile(jpth(srcDir,'optimize.f90'),srcLines)
        if int(Ngroup)==1:
           fpp.blck('cuts',cuts[0].src,srcLines)
        else:
           fpp.blck('cuts',cuts[int(groups[0])].src,srcLines)
        fpp.blck('pars',pars.src,srcLines)
        fpp.blck('processes',optimProcLines,srcLines)
        fpp.blck('tmds',tmdLines,srcLines)
        fpp.blck('dumpPath',frmt_char(jpth(targetDir,procDir[ii],'MC_'),cont=False),srcLines)
        fpp.xpnd('Nfinst',str(len(flavor)-3),srcLines)
        fpp.xpnd('Nflavors',Nflavors,srcLines)
        fpp.xpnd('lhaSet',lhaSet,srcLines)
        fpp.xpnd('Qsqueez',Qsqueez,srcLines)
        fpp.xpnd('offshell',print_array(offshell[0:2]),srcLines)
        fpp.xpnd('pNonQCD',print_array(pNonQCD),srcLines)
        fpp.xpnd('thrs',thrs,srcLines)
        fpp.xpnd('NstepGrid',NstepGrid,srcLines)
        fpp.xpnd('Nstep',Nstep,srcLines)
        fpp.xpnd('instPartons',instPartons[ii],srcLines)
        fpp.xpnd('partonB',particle(flavor[0]),srcLines)
        fpp.xpnd('partonA',particle(flavor[1]),srcLines)
        fpp.case('helicity',optimHel,srcLines)
        fpp.case('partlumi',partlumi,srcLines)
        fpp.case('fluxFactor',fluxFactor,srcLines)
        fpp.case('pdfTypeA',pdfTypeA,srcLines)
        fpp.case('pdfTypeB',pdfTypeB,srcLines)
        fpp.case('withTMDlib',withTMDlib,srcLines)
        fpp.case('withPDFxTMD0',withPDFxTMDTMD,srcLines)
        fpp.case('withPDFxTMDB','no',srcLines)
        fpp.case('withPDFxTMDCPDFA',withPDFxTMDCPDFA,srcLines)
        fpp.case('withPDFxTMDCPDFB',withPDFxTMDCPDFB,srcLines)
        fpp.case('withPDFxTMDCoupling',withPDFxTMDCoupling,srcLines)
        if withTMDlib=='yes':
            fpp.xpnd('TMDlibSet',TMDlibSet[0],srcLines)
            fpp.xpnd('kfA',TMDlibSet[2],srcLines)
            fpp.xpnd('kfB',TMDlibSet[1],srcLines)
        if withPDFxTMD=='yes':
            fpp.xpnd('PDFxTMDSetA',PDFxTMDSet[0],srcLines)
            fpp.xpnd('PDFxTMDkfA',PDFxTMDSet[2],srcLines)
            fpp.xpnd('PDFxTMDkfB',PDFxTMDSet[1],srcLines)
            fpp.xpnd('PDFxTMDMemberA',PDFxTMDSet[3],srcLines)
            fpp.xpnd('PDFxTMDSetB',PDFxTMDSetB[0],srcLines)
            fpp.xpnd('PDFxTMDMemberB',PDFxTMDSetB[3],srcLines)
        if withPDFxTMDCoupling=='yes':
            fpp.xpnd('PDFxTMDCouplingSet',PDFxTMDCouplingSet,srcLines)
        fpp.case('withLHAPDF',withLHAPDF,srcLines)
    
        fpp.incl(buildDir,srcLines)
        fpp.xpnd('realknd2','real(kind(1d0))',srcLines)
        fpp.xpnd('realknd1','real(kind(1e0))',srcLines)
        fpp.xpnd('complex2','complex(kind(1d0))',srcLines)
        fpp.xpnd('integer1','integer(selected_int_kind(4))',srcLines)
        fpp.xpnd('integer0','integer(selected_int_kind(2))',srcLines)

        #mainSrc = jpth(procDir[ii],'main'+procNum[ii])
        mainSrc = jpth(procDir[ii],'main')
        cr.wfile(jpth(targetDir,procDir[ii],'main.f90'),srcLines,delPattern='')
        compile_main(jpth(targetDir,procDir[ii],'main.f90'),procNum[ii])
        exec_all.append(jpth('array['+str(ii).rjust(3)+']="$prefix',procDir[ii],'main'+procNum[ii]+'.out') \
                        +' seed='+seed+' Noptim='+Noptim \
                        +' > '+jpth('$prefix',procDir[ii],'output')+'"\n')
        recompile_all.append('array['+str(ii+1).rjust(3)+']=' \
                   +'"'+' '.join(compile_main_arglist(jpth(targetDir,procDir[ii],'main.f90'),procNum[ii]))+'"\n')

    srcLines = ed.lines_from_file(jpth(srcDir,'run_optimize.sh'))
    ed.replace_block(['#Begin process list','#End process list'],exec_all,srcLines)
    cr.wfile(jpth(targetDir,'optimize.sh'),srcLines)
    srcLines = ed.lines_from_file(jpth(srcDir,'run_recompile.sh'))
    ed.replace_block(['#Begin main list','#End main list'],recompile_all,srcLines)
    cr.wfile(jpth(targetDir,'recompile.sh'),srcLines)

    info4eventFile = []
    info4eventFile.append('  write(eventUnit,\'(A)\') "LIST OF PROCESSES: Nf=Nflavor: '+Nflavors+'"\n')
    ii = 0
    for line in procList:
        line = re.sub(r'groups.*','',line)
        line = '  write(eventUnit,\'(A)\') "LIST OF PROCESSES: process '+procNum[ii]+': '+line
        if partlumi=='combined': line = line+'partlumi = pdf_'+instPartons[ii]
        line = line+'"\n'
        info4eventFile.append(line)
        ii = ii+1
    if partlumi=='combined':
        info4eventFile.append('  write(eventUnit,\'(A)\') '
        +'"LIST OF PROCESSES: The partlumi are defined in src/partlumi.f90"\n')
    info4eventFile.append('  write(eventUnit,\'(A)\') "NUMBER OF GROUPS: '+Ngroup+'"\n')
    info4eventFile.append('  write(eventUnit,\'(A)\') "NUMBER OF FINAL-STATE PARTICLES: '+' '.join(Nfinst)+'"\n')
    info4eventFile = info4eventFile+cuts[0].info+pars.info

    srcLines = []
    cr.addfile(jpth(srcDir,'mainMC.f90'),srcLines)
    
    fpp.blck('cuts',cuts[0].src,srcLines)
    fpp.blck('pars',pars.src,srcLines)
    fpp.blck('processes',mainProcLines,srcLines)
    fpp.blck('tmds',tmdLines,srcLines)
    fpp.blck('info4eventFile',info4eventFile,srcLines)
    fpp.xpnd('Ngroup',Ngroup,srcLines)
    fpp.xpnd('Nfinst',print_array(Nfinst),srcLines)
    fpp.xpnd('NprocTot',str(len(procList)),srcLines)
    fpp.xpnd('Nflavors',Nflavors,srcLines)
    fpp.xpnd('lhaSet',lhaSet,srcLines)
    fpp.xpnd('Qsqueez',Qsqueez,srcLines)
    fpp.xpnd('offshell',print_array(offshell[0:2]),srcLines)
    fpp.xpnd('sigma_eff',sigma_eff,srcLines)
    fpp.case('helicity',helicity,srcLines)
    fpp.case('partlumi',partlumi,srcLines)
    fpp.case('fluxFactor',fluxFactor,srcLines)
    fpp.case('pdfTypeA',pdfTypeA,srcLines)
    fpp.case('pdfTypeB',pdfTypeB,srcLines)
    fpp.case('withTMDlib',withTMDlib,srcLines)
    fpp.case('withPDFxTMD0',withPDFxTMDTMD,srcLines)
    fpp.case('withPDFxTMDB','no',srcLines)
    fpp.case('withPDFxTMDCPDFA',withPDFxTMDCPDFA,srcLines)
    fpp.case('withPDFxTMDCPDFB',withPDFxTMDCPDFB,srcLines)
    fpp.case('withPDFxTMDCoupling',withPDFxTMDCoupling,srcLines)
    if withTMDlib=='yes': 
        fpp.xpnd('TMDlibSet',TMDlibSet[0],srcLines)
        fpp.xpnd('kfA',TMDlibSet[2],srcLines)
        fpp.xpnd('kfB',TMDlibSet[1],srcLines)
    if withPDFxTMD=='yes':
        fpp.xpnd('PDFxTMDSetA',PDFxTMDSet[0],srcLines)
        fpp.xpnd('PDFxTMDkfA',PDFxTMDSet[2],srcLines)
        fpp.xpnd('PDFxTMDkfB',PDFxTMDSet[1],srcLines)
        fpp.xpnd('PDFxTMDMemberA',PDFxTMDSet[3],srcLines)
        fpp.xpnd('PDFxTMDSetB',PDFxTMDSetB[0],srcLines)
        fpp.xpnd('PDFxTMDMemberB',PDFxTMDSetB[3],srcLines)
    if withPDFxTMDCoupling=='yes':
        fpp.xpnd('PDFxTMDCouplingSet',PDFxTMDCouplingSet,srcLines)
    fpp.case('withLHAPDF',withLHAPDF,srcLines)
    
    fpp.incl(buildDir,srcLines)
    fpp.xpnd('realknd2','real(kind(1d0))',srcLines)
    fpp.xpnd('realknd1','real(kind(1e0))',srcLines)
    fpp.xpnd('complex2','complex(kind(1d0))',srcLines)
    fpp.xpnd('integer1','integer(selected_int_kind(4))',srcLines)
    fpp.xpnd('integer0','integer(selected_int_kind(2))',srcLines)
    cr.wfile(jpth(targetDir,'main.f90'),srcLines,delPattern='')
    compile_main(jpth(targetDir,'main.f90'))
