#!/usr/bin/env python

from create_def import *

def add_tmdpdf(beam,parton,gridfile):
    beamArg = '-1'
    #if beam in ('gridB','B','-2'): beamArg = '-2'
    if beam==2: beamArg = '-2'
    line = '  call add_tmdpdf('
    if parton in ('g','u','u~','d','d~','s','s~','c','c~','b','b~'):
        line = line + beamArg+','
        line = line + particle(parton)+','
    else:
        line = line + '\''+parton+'\''+','
    line = line +'\''+gridfile+'\')\n'
    return line


def parse_procLine(procLine):
    flavor = ['none']
    factor = 'none'
    groups = ['none']
    pNonQCD = ['none']
    if re.search('=',procLine):
        if re.search(',',procLine):
            argList = procLine.split(',')
            flavor = argList[0].split()
            for arg in argList:
                if cr.lhs('factor','=',arg): factor = re.sub(r'Nf','Nflavors',cr.rhs('=',arg))
                if cr.lhs('groups','=',arg): groups = cr.rhs('=',arg).split()
                if cr.lhs('pNonQCD','=',arg): pNonQCD = cr.rhs('=',arg).split()
        else:
            flavor = re.sub(r' *factor.*','',procLine).split()
            factor = re.sub(r'Nf','Nflavors',re.sub(r' *groups.*','',cr.rhs('=',procLine)))
            groups = re.sub(r' *pNonQCD.*','',cr.rhs('=',cr.rhs('=',procLine))).split()
            pNonQCD = cr.rhs('=',cr.rhs('=',cr.rhs('=',procLine))).split()
    else:
        flavor = procLine.split()
    return flavor,factor,groups,pNonQCD


def pureDISflav(flavor):
    #psFlav = flavor.copy()
    psFlav = flavor[:]
    psFlav[0] = 'A'
    del psFlav[-1]
    return psFlav


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
    ii = ii+1

offshell = ['0','0']
Ngroup = '1'
Nfinst = []
Nflavors = '5'
lhaSet = 'MSTW2008lo68'
lhaSetSet = False
withLHAPDF = 'yes'
TMDlibSet = []
TMDlibSetB = []
PDFxTMDSet = []
PDFxTMDSetB = []
PDFxTMDType = 'tmd'
PDFxTMDCoupling = 'no'
PDFxTMDCouplingSet = ''
tmdTableDir = '/home/user0/kTfac/krzysztof/'
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
Ecm = '0'
EbeamB = '0'
EbeamA = '0'
defNonQCD = ['0','0','0']
pars = cuts_and_pars.parameters()
instOption = '1'
xAlimits = ['0d0','1d0']
itmdf = 'no'
leadingColor = 'no'
xBbin = ['','']
QsqBin = ['','']
pureDIS = False
DISF2 = 'no'
withMINCAS = 'no'

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
    if cr.lhs(  'TMDlib[Ss]et','=',line): TMDlibSet = cr.rhs('=',line).split()
    if cr.lhs(  'TMDlib[Ss]et +A','=',line): TMDlibSet = cr.rhs('=',line).split()
    if cr.lhs(  'TMDlib[Ss]et +B','=',line): TMDlibSetB = cr.rhs('=',line).split()
    if cr.lhs( 'PDFxTMD[Ss]et','=',line): PDFxTMDSet = cr.rhs('=',line).split()
    if cr.lhs( 'PDFxTMD[Ss]et +A','=',line): PDFxTMDSet = cr.rhs('=',line).split()
    if cr.lhs( 'PDFxTMD[Ss]et +B','=',line): PDFxTMDSetB = cr.rhs('=',line).split()
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
    if cr.lhs(      'Nflavors','=',line): Nflavors = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(      'offshell','=',line): offshell = cr.rhs('=',line).split()
    if cr.lhs(     'sigma_eff','=',line): sigma_eff = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(           'Ecm','=',line): Ecm = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(         'Esoft','=',line): Esoft = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(        'EbeamB','=',line): EbeamB = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(        'EbeamA','=',line): EbeamA = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs('Ebeam[Nn]eg[Rr]ap','=',line): EbeamB = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs('Ebeam[Pp]os[Rr]ap','=',line): EbeamA = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs('Ebeam[Ee]lectron','=',line): EbeamB = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs('Ebeam[Hh]adron','=',line): EbeamA = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(   'part[Ll]umi','=',line): partlumi = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs(    'p[Nn]onQCD','=',line): defNonQCD = cr.rhs('=',line).split()
    if cr.lhs('flux *[fF]actor','=',line): fluxFactor = cr.rhs('=',line).split(' ',)[0]
    if cr.lhs('print skeleton','=',line): printskel = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(   'inst option','=',line): instOption = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(         'itmd' ,'=',line): itmdf = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(         'itmdf','=',line): itmdf = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(  'leadingColor','=',line): leadingColor = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(      'xAlimits','=',line): xAlimits = cr.rhs('=',line).split()
    if cr.lhs(  'xBjorken bin','=',line): xBbin = cr.rhs('=',line).split()
    if cr.lhs(   'Qsquare bin','=',line): QsqBin = cr.rhs('=',line).split()
    if cr.lhs(         'DISF2','=',line): DISF2 = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(   'with MINCAS','=',line): withMINCAS = cr.rhs('=',line).split(' ',1)[0]
    if cr.lhs(      'helicity','=',line): 
        helicity = cr.rhs('=',line).split(' ',1)[0]
        optimHel = helicity
    if cr.lhs('optimization helicity','=',line):
        optimHel = cr.rhs('=',line).split(' ',1)[0]

    if cr.lhs(    'mass','=',line): pars.mass(    cr.rhs('=',line))
    if cr.lhs(  'switch','=',line): pars.switch(  cr.rhs('=',line))
    if cr.lhs('coupling','=',line): pars.coupling(cr.rhs('=',line))

pars.esoft(Esoft)
if EbeamB=='0' and EbeamA=='0':
    if Ecm=='0':
        cr.prnt('ERROR: Ecm or both beam energies must be non-zero')
        sys.exit()
    pars.energy(Sdouble(Ecm))
    pars.ePosRap('Ecm/2')
    pars.eNegRap('Ecm/2')
elif EbeamB=='0' or EbeamA=='0':
    cr.prnt('ERROR: both beam energies must be non-zero')
    sys.exit()
else:
    if Ecm!='0':
        cr.prnt('ERROR: Ecm and beam energies cannot be set simultaneously')
        sys.exit()
    pars.eNegRap(Sdouble(EbeamB))
    pars.ePosRap(Sdouble(EbeamA))
    pars.energy('2*sqrt(EposRap*EnegRap)')
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
    if cr.lhs('cut source','=',line):
         rhsline = cr.rhs('=',line)
         rhs = rhsline.split()
         if rhs[0]=='group':
             cuts[int(rhs[1])].add_extra_cut(re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
         else:
             cuts[0].add_extra_cut(rhsline)
    if cr.lhs('weight factor','=',line):
         rhsline = cr.rhs('=',line)
         rhs = rhsline.split()
         if rhs[0]=='group':
             cuts[int(rhs[1])].add_extra_weight(re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
         else:
             cuts[0].add_extra_weight(rhsline)
    if cr.lhs('scale','=',line): 
         rhsline = cr.rhs('=',line)
         rhs = rhsline.split()
         if rhs[0]=='group':
             cuts[int(rhs[1])].set_scale('renScale',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
             cuts[int(rhs[1])].set_scale('scaleA',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
             cuts[int(rhs[1])].set_scale('scaleB',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline))
         elif rhs[0]=='entry':
             cuts[0].set_scale('renScale',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline),group=rhs[1])
             cuts[0].set_scale('scaleA',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline),group=rhs[1])
             cuts[0].set_scale('scaleB',re.sub(' *'+rhs[0]+' *'+rhs[1]+' *','',rhsline),group=rhs[1])
         else:
             cuts[0].set_scale('renScale',rhsline)
             cuts[0].set_scale('scaleA',rhsline)
             cuts[0].set_scale('scaleB',rhsline)
    if cr.lhs('ren.scale','=',line) or cr.lhs('renormalization scale','=',line): 
         rhsline = cr.rhs('=',line)
         cuts[0].set_scale('renScale',rhsline)
    if cr.lhs('scaleA','=',line): 
         rhsline = cr.rhs('=',line)
         cuts[0].set_scale('scaleA',rhsline)
    if cr.lhs('scaleB','=',line): 
         rhsline = cr.rhs('=',line)
         cuts[0].set_scale('scaleB',rhsline)

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
        if instPartons[-1] in ('r_q','qb_rb','rb_qb'):
            instPartons[-1]='q_r' 
        if instPartons[-1] in ('r_qb','qb_r','rb_q'):
            instPartons[-1]='q_rb' 
        if instPartons[-1]=='A_A':
            partlumi='photons' 
        if word[0][:3]=='DIS':
            if not partlumi in ('individual','DIS'):
                cr.prnt('ERROR from '+KATIEpath+'/create.py:')
                cr.prnt('  all processes or none of them must be DIS')
                sys.exit()
            if word[0]=='DIS':
                sgn='-' 
                procList[-1] = re.sub('DIS','e-',procList[-1])
            else:
                sgn=word[0][3]
                procList[-1] = re.sub('DIS','e',procList[-1])
            partlumi = 'DIS'
            scatproc = 'DIS'+sgn
            isComma = re.search(',',procList[-1])
            if isComma:
                firstComma = isComma.start()
                procList[-1] = procList[-1][:firstComma]+'e'+sgn+' '+procList[-1][firstComma:]
            else:
                if re.search('=',procList[-1]):
                    cr.prnt('ERROR: format for process lines without commas is not allowed for DIS')
                    sys.exit()
                procList[-1] = procList[-1]+' e'+sgn
            pureDIS = xBbin[0]!='' or QsqBin[0]!=''
        if len(Nfinst)==0: Nfinst = [str(len(parse_procLine(procList[-1])[0])-3)]
        ii = ii+1

if   partlumi=='photons': scatproc = 'photon-photon'
elif partlumi=='DIS'    : pass
else:                     scatproc = 'hadron-hadron'

partonSymb = ['g','u','u~','d','d~','s','s~','c','c~','b','b~'\
             ,'qg1','qg2','qg3','gg1','gg2','gg3','gg4','gg5','gg6','gg7']
partonLines = [[[] for _ in range(20)] for _ in range(3)]
withTMDlib = 'no'
withTMDlibB = 'no'
withPDFxTMD = 'no'
withPDFxTMDB = 'no'
withPDFxTMDTMD = 'no'
withPDFxTMDCPDFA = 'no'
withPDFxTMDCPDFB = 'no'
withPDFxTMDCoupling = PDFxTMDCoupling
tmdLines = []
if len(PDFxTMDSet)>0 or len(PDFxTMDSetB)>0:
    withPDFxTMD = 'yes'
    withPDFxTMDTMD = 'yes'
    if len(PDFxTMDSet)<2: PDFxTMDSet = PDFxTMDSet + ['2212','2212']
    if len(PDFxTMDSet)<3: PDFxTMDSet = PDFxTMDSet + [PDFxTMDSet[1]]
    if len(PDFxTMDSet)<4: PDFxTMDSet = PDFxTMDSet + ['0']
    if withPDFxTMDCoupling=='yes' and PDFxTMDCouplingSet=='':
        PDFxTMDCouplingSet = PDFxTMDSet[0]
    if len(PDFxTMDSetB)==0:
        pdfTypeA = 'PDFxTMD0'
        pdfTypeB = 'PDFxTMD0'
        if PDFxTMDType=='cpdf':
            PDFxTMDSetB = PDFxTMDSet
    else:
        withPDFxTMDB = 'yes'
        pdfTypeA = 'PDFxTMDA'
        pdfTypeB = 'PDFxTMDB'
        if len(PDFxTMDSetB)<2: PDFxTMDSetB = PDFxTMDSetB + ['2212','2212']
        if len(PDFxTMDSetB)<3: PDFxTMDSetB = PDFxTMDSetB + [PDFxTMDSetB[1]]
        if len(PDFxTMDSetB)<4: PDFxTMDSetB = PDFxTMDSetB + ['0']
        if len(PDFxTMDSetB)>1: PDFxTMDSet[1] = PDFxTMDSetB[1]
        if len(PDFxTMDSetB)>2: PDFxTMDSet[2] = PDFxTMDSetB[2]
elif len(TMDlibSet)==0:
    pdfTypeA = ''
    pdfTypeB = ''
    for line in inputLines:
        if cr.lhs('tmd[Pp]ath','=',line) or cr.lhs('tmdTableDir','=',line):
            tmpPath = cr.rhs('=',line).split()
            tmdLines.append('  call set_tmdpath(&\n')
            tmdLines = tmdLines + frmt_char(tmpPath[0],cont=True)
            tmdLines.append('  )\n')
        if cr.lhs('tmd[Pp]df','=',line) or cr.lhs('tmd[Pp]df +.','=',line):
            tmpFile = cr.rhs('=',line).split()
            tmdExists(tmpPath,tmpFile)
            if   cr.lhs('tmd[Pp]df +A','=',line):
                if pdfTypeA!='' and pdfTypeA!='gridA':
                    cr.prnt('ERROR: inconsistent beam grids')
                    sys.exit()
                pdfTypeA = 'gridA'
                beamId = 1
            elif cr.lhs('tmd[Pp]df +B','=',line):
                if pdfTypeB!='' and pdfTypeB!='gridB':
                    cr.prnt('ERROR: inconsistent beam grids')
                    sys.exit()
                pdfTypeB = 'gridB'
                beamId = 2
            else:
                if (pdfTypeA!='' and pdfTypeA!='grid0') or (pdfTypeB!='' and pdfTypeB!='grid0'):
                    cr.prnt('ERROR: inconsistent beam grids')
                    sys.exit()
                pdfTypeA = 'grid0'
                pdfTypeB = 'grid0'
                beamId = 0
            ii = partonSymb.index(tmpFile[0])
            partonLines[beamId][ii].append('  call set_tmdpath(&\n')
            partonLines[beamId][ii] = partonLines[beamId][ii] + frmt_char(tmpPath[0],cont=True)
            partonLines[beamId][ii].append('  )\n')
            partonLines[beamId][ii].append( add_tmdpdf(beamId,tmpFile[0],tmpFile[1]) )
            tmdLines.append(                add_tmdpdf(beamId,tmpFile[0],tmpFile[1]) )
else:
    withTMDlib = 'yes'
    if len(TMDlibSetB)==0:
        pdfTypeA = 'TMDlib0'
        pdfTypeB = 'TMDlib0'
        if len(TMDlibSet)<2: TMDlibSet = TMDlibSet + ['2212','2212']
        if len(TMDlibSet)<3: TMDlibSet = TMDlibSet + [TMDlibSet[1]]
    else:
        withTMDlibB = 'yes'
        pdfTypeA = 'TMDlibA'
        pdfTypeB = 'TMDlibB'
        if len(TMDlibSet)<2: TMDlibSet = TMDlibSet + ['2212','2212']
        if len(TMDlibSet)<3: TMDlibSet = TMDlibSet + [TMDlibSet[1]]
    if len(TMDlibSetB)>1: TMDlibSet[1] = TMDlibSetB[1]
    TMDlibKey = 'char'
    if TMDlibSet[0].isdigit(): TMDlibKey = 'int'
if len(offshell)==1:
    if offshell[0] in ('yes','Yes','YES','Y','y'): offshell[0]='1'
    offshell.insert(0,'0')
if PDFxTMDType=='cpdf':
    withPDFxTMDTMD = 'no'
if withLHAPDF=='yes' and PDFxTMDType=='cpdf' and withPDFxTMDCoupling=='yes' and not lhaSetSet:
    withLHAPDF = 'no'
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
compileOptions = []
if withMINCAS=='yes': compileOptions.append('withMINCAS')
recompile_all = ['Nmain='+str(1+len(procList))+'\n' \
                ,'array[  0]="'+cmplcmd_file_as_str(jpth(targetDir,'main.f90'),options=compileOptions)+'"\n']

mainKinLines = []
for ii in range(0,len(Nfinst)):
    word = '['
    if offshell[1]=='0': word = word+'1'
    else:                word = word+'2'
    for jj in range(0,int(Nfinst[ii])): word = word+',0'
    if offshell[0]=='0': word = word+',1'
    else:                word = word+',2'
    word = word + ']'
    mainKinLines.append('  call add_kinematics( kinID('+str(ii+1)+')'\
                       +' ,'+str(int(Nfinst[ii])+2)+' ,'+word+' )\n')

instB = 0
instA = 1 
mainProcLines = []
for ii in range(0,len(procList)):
    pfd.mkdir(jpth(targetDir,procDir[ii]))

    flavor,factor,groups,pNonQCD = parse_procLine(procList[ii])
    if factor=='none': factor='1'
    if groups==['none']: groups = ['1']
    if pNonQCD==['none']: pNonQCD = defNonQCD

    if not pureDIS: mainProcLines.append('  proc('+procNum[ii]+')%group(1:'+str(len(groups))+') = '+print_array(groups)+'\n')
    mainProcLines.append('  proc('+procNum[ii]+')%pNonQCD = '+print_array(pNonQCD)+'\n')
    if not pureDIS: mainProcLines.append('  proc('+procNum[ii]+')%Nfinst = '+str(len(flavor)-3)+'\n')
    mainProcLines.append('  call proc('+procNum[ii]+')%flav%put('+print_process(flavor)+',get_anti)\n')
    if pureDIS: mainProcLines.append('  call proc('+procNum[ii]+')%psflav%put('+print_process(pureDISflav(flavor))+',get_anti)\n')
    mainProcLines.append('  proc('+procNum[ii]+')%factor = '+factor+'\n')
    mainProcLines.append('  proc('+procNum[ii]+')%label = &\n')
    mainProcLines = mainProcLines + frmt_char(jpth(targetDir,procDir[ii],'MC_'),cont=False)
    if partlumi=='combined':
      mainProcLines.append('  proc('+procNum[ii]+')%pdf => pdf_'+instPartons[ii]+'\n')
    else:
      mainProcLines.append('  proc('+procNum[ii]+')%parton(-2:-1) = '+print_inst(flavor)+'\n')

    optimKinLines = re.sub( 'kinID\(.*?\)', 'kinID' ,mainKinLines[int(groups[0])-1] )
    optimProcLines = []
    optimProcLines.append('  call flavor%put('+print_process(flavor)+',get_anti)\n')
    if pureDIS: optimProcLines.append('  call psflavor%put('+print_process(pureDISflav(flavor))+',get_anti)\n')
    optimProcLines.append('  ampFac = '+factor+'\n')

    if partlumi=='combined' or withTMDlib=='yes' or withPDFxTMD=='yes' or itmdf=='yes':
        optimTmdLines = tmdLines
    elif partlumi=='DIS':
        optimTmdLines = []
        if pdfTypeA=='grid0':
            optimTmdLines = optimTmdLines + partonLines[0][partonSymb.index(flavor[instA])]
        if pdfTypeA=='gridA':
            optimTmdLines = optimTmdLines + partonLines[1][partonSymb.index(flavor[instA])]
    else:
        optimTmdLines = []
        if pdfTypeA+pdfTypeB=='grid0grid0':
            optimTmdLines = optimTmdLines + partonLines[0][partonSymb.index(flavor[instA])]
            if flavor[instA]!=flavor[instB]:
                optimTmdLines = optimTmdLines + partonLines[0][partonSymb.index(flavor[instB])]
        if pdfTypeA+pdfTypeB=='gridAgridB':
            optimTmdLines = optimTmdLines + partonLines[1][partonSymb.index(flavor[instA])]
            optimTmdLines = optimTmdLines + partonLines[2][partonSymb.index(flavor[instB])]
        if pdfTypeA+pdfTypeB=='grid0LHAPDF':
            optimTmdLines = optimTmdLines + partonLines[0][partonSymb.index(flavor[instA])]
        if pdfTypeA+pdfTypeB=='gridALHAPDF':
            optimTmdLines = optimTmdLines + partonLines[1][partonSymb.index(flavor[instA])]
        if pdfTypeA+pdfTypeB=='LHAPDFgrid0':
            optimTmdLines = optimTmdLines + partonLines[0][partonSymb.index(flavor[instB])]
        if pdfTypeA+pdfTypeB=='LHAPDFgridB':
            optimTmdLines = optimTmdLines + partonLines[2][partonSymb.index(flavor[instB])]

    srcLines = []
    if pureDIS:
        cr.addfile(jpth(srcDir,'optimizeDIS.f90'),srcLines)
    else:
        cr.addfile(jpth(srcDir,'optimize.f90'),srcLines)
    if int(Ngroup)==1:
        fpp.blck('cuts',frmt_list(cuts[0].src+cuts[0].extracuts),srcLines)
        fpp.blck('weights',frmt_list(cuts[0].extraweights),srcLines)
    else:
        fpp.blck('cuts',frmt_list(cuts[int(groups[0])].src+cuts[int(groups[0])].extracuts),srcLines)
        fpp.blck('weights',frmt_list(cuts[int(groups[0])].extraweights),srcLines)
    fpp.blck('pars',frmt_list(pars.src),srcLines)
    fpp.blck('kinematics',optimKinLines,srcLines)
    fpp.blck('processes',optimProcLines,srcLines)
    fpp.blck('tmds',optimTmdLines,srcLines)
    fpp.blck('dumpPath',frmt_char(jpth(targetDir,procDir[ii],'MC_'),cont=False),srcLines)
    fpp.xpnd('Nfinst',str(len(flavor)-3),srcLines)
    fpp.xpnd('Nflavors',Nflavors,srcLines)
    fpp.xpnd('lhaSet',lhaSet,srcLines)
    fpp.xpnd('offshell',print_array(offshell[0:2]),srcLines)
    fpp.xpnd('pNonQCD',print_array(pNonQCD),srcLines)
    fpp.xpnd('thrs',thrs,srcLines)
    fpp.xpnd('NstepGrid',NstepGrid,srcLines)
    fpp.xpnd('Nstep',Nstep,srcLines)
    fpp.xpnd('instPartons',instPartons[ii],srcLines)
    fpp.xpnd('instOption',instOption,srcLines)
    fpp.xpnd('xAmin',xAlimits[0],srcLines)
    fpp.xpnd('xAmax',xAlimits[1],srcLines)
    fpp.xpnd('partonB',particle(flavor[instB]),srcLines)
    fpp.xpnd('partonA',particle(flavor[instA]),srcLines)
    fpp.xpnd('xBmin',Sdouble(xBbin[0]),srcLines)
    fpp.xpnd('xBmax',Sdouble(xBbin[1]),srcLines)
    fpp.xpnd('QsqMin',Sdouble(QsqBin[0]),srcLines)
    fpp.xpnd('QsqMax',Sdouble(QsqBin[1]),srcLines)
    fpp.case('helicity',optimHel,srcLines)
    fpp.case('partlumi',partlumi,srcLines)
    fpp.case('fluxFactor',fluxFactor,srcLines)
    fpp.case('itmdf',itmdf,srcLines)
    fpp.case('leadingColor',leadingColor,srcLines)
    fpp.case('DISF2',DISF2,srcLines)
    fpp.case('withMINCAS',withMINCAS,srcLines)
    fpp.case('pdfTypeA',pdfTypeA,srcLines)
    fpp.case('pdfTypeB',pdfTypeB,srcLines)
    fpp.case('pdfTypes',pdfTypeA+'_'+pdfTypeB,srcLines)
    fpp.case('withTMDlib0',withTMDlib,srcLines)
    fpp.case('withTMDlibB',withTMDlibB,srcLines)
    fpp.case('withPDFxTMD0',withPDFxTMDTMD,srcLines)
    fpp.case('withPDFxTMDB',withPDFxTMDB,srcLines)
    fpp.case('withPDFxTMDCPDFA',withPDFxTMDCPDFA,srcLines)
    fpp.case('withPDFxTMDCPDFB',withPDFxTMDCPDFB,srcLines)
    fpp.case('withPDFxTMDCoupling',withPDFxTMDCoupling,srcLines)
    if withTMDlib=='yes':
        fpp.case('TMDlibKey',TMDlibKey,srcLines)
        fpp.xpnd('TMDlibSetA',TMDlibSet[0],srcLines)
        fpp.xpnd('kfA',TMDlibSet[2],srcLines)
        fpp.xpnd('kfB',TMDlibSet[1],srcLines)
    if withTMDlibB=='yes':
        fpp.xpnd('TMDlibSetB',TMDlibSetB[0],srcLines)
    if withPDFxTMD=='yes':
        fpp.xpnd('PDFxTMDSetA',PDFxTMDSet[0],srcLines)
        fpp.xpnd('PDFxTMDkfA',PDFxTMDSet[2],srcLines)
        fpp.xpnd('PDFxTMDkfB',PDFxTMDSet[1],srcLines)
        fpp.xpnd('PDFxTMDMemberA',PDFxTMDSet[3],srcLines)
    if withPDFxTMDB=='yes' or withPDFxTMDCPDFB=='yes':
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
    compile_file(jpth(targetDir,procDir[ii],'main.f90')\
         ,target=jpth(targetDir,procDir[ii],'main'+procNum[ii]),options=compileOptions)
    exec_all.append(jpth('array['+str(ii).rjust(3)+']="$prefix',procDir[ii],'main'+procNum[ii]+'.out') \
                    +' seed='+seed+' Noptim='+Noptim \
                    +' > '+jpth('$prefix',procDir[ii],'output')+'"\n')
    recompile_all.append('array['+str(ii+1).rjust(3)+']=' \
               +'"'+cmplcmd_file_as_str(jpth(targetDir,procDir[ii],'main.f90')\
                                       ,target=jpth(targetDir,procDir[ii],'main'+procNum[ii])\
                                       ,options=compileOptions)+'"\n')

srcLines = ed.lines_from_file(jpth(srcDir,'run_optimize.sh'))
ed.replace_block(['#Begin process list','#End process list'],exec_all,srcLines)
cr.wfile(jpth(targetDir,'optimize.sh'),srcLines)
srcLines = ed.lines_from_file(jpth(srcDir,'run_recompile.sh'))
ed.replace_block(['#Begin main list','#End main list'],recompile_all,srcLines)
cr.wfile(jpth(targetDir,'recompile.sh'),srcLines)

info4eventFile = []
info4eventFile.append('  write(eventUnit,\'(A)\') "SCATTERING PROCESS: '+scatproc+'"\n')
info4eventFile.append('  write(eventUnit,\'(A)\') "LIST OF PROCESSES: Nf=Nflavor: '+Nflavors+'"\n')
ii = 0
for line in procList:
    line = re.sub(r'groups.*','',line)
    line = '  write(eventUnit,\'(A)\') "LIST OF PROCESSES: process '+procNum[ii]+': '+line
    if partlumi=='combined': line = line+' partlumi = pdf_'+instPartons[ii]
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
if pureDIS:
    cr.addfile(jpth(srcDir,'mainDIS.f90'),srcLines)
else:
    cr.addfile(jpth(srcDir,'mainMC.f90'),srcLines)

fpp.blck('cuts',frmt_list(cuts[0].src+cuts[0].extracuts),srcLines)
fpp.blck('weights',frmt_list(cuts[0].extraweights),srcLines)
fpp.blck('pars',frmt_list(pars.src),srcLines)
fpp.blck('kinematics',mainKinLines,srcLines)
fpp.blck('processes',mainProcLines,srcLines)
fpp.blck('tmds',tmdLines,srcLines)
fpp.blck('info4eventFile',info4eventFile,srcLines)
fpp.xpnd('Ngroup',Ngroup,srcLines)
fpp.xpnd('Nfinst',print_array(Nfinst),srcLines)
fpp.xpnd('NfinstPureDIS',str(Nfinst[0]),srcLines)
fpp.xpnd('NprocTot',str(len(procList)),srcLines)
fpp.xpnd('Nflavors',Nflavors,srcLines)
fpp.xpnd('lhaSet',lhaSet,srcLines)
fpp.xpnd('offshell',print_array(offshell[0:2]),srcLines)
fpp.xpnd('sigma_eff',sigma_eff,srcLines)
fpp.xpnd('instOption',instOption,srcLines)
fpp.xpnd('xAmin',xAlimits[0],srcLines)
fpp.xpnd('xAmax',xAlimits[1],srcLines)
fpp.xpnd('xBmin',Sdouble(xBbin[0]),srcLines)
fpp.xpnd('xBmax',Sdouble(xBbin[1]),srcLines)
fpp.xpnd('QsqMin',Sdouble(QsqBin[0]),srcLines)
fpp.xpnd('QsqMax',Sdouble(QsqBin[1]),srcLines)
fpp.case('helicity',helicity,srcLines)
fpp.case('partlumi',partlumi,srcLines)
fpp.case('fluxFactor',fluxFactor,srcLines)
fpp.case('itmdf',itmdf,srcLines)
fpp.case('leadingColor',leadingColor,srcLines)
fpp.case('DISF2',DISF2,srcLines)
fpp.case('withMINCAS',withMINCAS,srcLines)
fpp.case('pdfTypeA',pdfTypeA,srcLines)
fpp.case('pdfTypeB',pdfTypeB,srcLines)
fpp.case('pdfTypes',pdfTypeA+'_'+pdfTypeB,srcLines)
fpp.case('withTMDlib0',withTMDlib,srcLines)
fpp.case('withTMDlibB',withTMDlibB,srcLines)
fpp.case('withPDFxTMD0',withPDFxTMDTMD,srcLines)
fpp.case('withPDFxTMDB',withPDFxTMDB,srcLines)
fpp.case('withPDFxTMDCPDFA',withPDFxTMDCPDFA,srcLines)
fpp.case('withPDFxTMDCPDFB',withPDFxTMDCPDFB,srcLines)
fpp.case('withPDFxTMDCoupling',withPDFxTMDCoupling,srcLines)
if withTMDlib=='yes': 
    fpp.case('TMDlibKey',TMDlibKey,srcLines)
    fpp.xpnd('TMDlibSetA',TMDlibSet[0],srcLines)
    fpp.xpnd('kfA',TMDlibSet[2],srcLines)
    fpp.xpnd('kfB',TMDlibSet[1],srcLines)
if withTMDlibB=='yes': 
    fpp.xpnd('TMDlibSetB',TMDlibSetB[0],srcLines)
if withPDFxTMD=='yes':
    fpp.xpnd('PDFxTMDSetA',PDFxTMDSet[0],srcLines)
    fpp.xpnd('PDFxTMDkfA',PDFxTMDSet[2],srcLines)
    fpp.xpnd('PDFxTMDkfB',PDFxTMDSet[1],srcLines)
    fpp.xpnd('PDFxTMDMemberA',PDFxTMDSet[3],srcLines)
if withPDFxTMDB=='yes' or withPDFxTMDCPDFB=='yes':
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
compile_file(jpth(targetDir,'main.f90'),options=compileOptions)

srcLines = cmplcmd_file_as_strlist(jpth(targetDir,'create_eventfile.f90'),options=compileOptions)
for idx, item in enumerate(srcLines):
   if re.search(targetDir,item):
       item = re.sub(targetDir,'$here',item)
       srcLines[idx] = item
   if re.search(buildDir,item):
       item = re.sub(buildDir,'$katieBuild',item)
       srcLines[idx] = item
srcLines = ['#!/bin/bash\n','here='+targetDir+'\n','katieBuild='+buildDir+'\n'] \
         + srcLines
srcLines = srcLines + ['#\n','if [ "$#" -gt 0 ]; then\n' \
                            ,'  $here/create_eventfile.out "$@"\n' \
                            ,'fi\n']
cr.wfile(jpth(targetDir,'create_eventfile.sh'),srcLines)
pfd.cp(jpth(srcDir,'create_eventfile.f90'),targetDir)
compile_file(jpth(targetDir,'create_eventfile.f90'),options=compileOptions)
