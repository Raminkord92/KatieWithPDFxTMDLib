from __future__ import print_function
import sys,re

def pdf_particle(word):
    if   word=='g':  return '0'
    elif word=='u':  return '1' 
    elif word=='u~': return '-1'
    elif word=='d':  return '2'
    elif word=='d~': return '-2'
    elif word=='s':  return '3'
    elif word=='s~': return '-3'
    elif word=='c':  return '4'
    elif word=='c~': return '-4'
    elif word=='b':  return '5'
    elif word=='b~': return '-5'
    elif word=='t':  return '6'
    elif word=='t~': return '-6'

def particle(word):
    if   word=='ve':    return 'eleNeu'
    elif word=='ve~':   return '-eleNeu'
    elif word=='e':     return 'eleon'  
    elif word=='e-':    return 'eleon'  
    elif word=='e+':    return '-eleon'
    elif word=='u':     return 'uQuark' 
    elif word=='u~':    return '-uQuark'
    elif word=='d':     return 'dQuark' 
    elif word=='d~':    return '-dQuark'
    elif word=='vmu':   return 'muNeu'  
    elif word=='vmu~':  return '-muNeu'
    elif word=='mu':    return 'muon'   
    elif word=='mu-':   return 'muon'   
    elif word=='mu+':   return '-muon'
    elif word=='c':     return 'cQuark' 
    elif word=='c~':    return '-cQuark'
    elif word=='s':     return 'sQuark' 
    elif word=='s~':    return '-sQuark'
    elif word=='vtau':  return 'tauNeu' 
    elif word=='vtau~': return '-tauNeu'
    elif word=='tau':   return 'tauon'  
    elif word=='tau-':  return 'tauon'  
    elif word=='tau+':  return '-tauon'
    elif word=='t':     return 'tQuark' 
    elif word=='t~':    return '-tQuark'
    elif word=='b':     return 'bQuark' 
    elif word=='b~':    return '-bQuark'
    elif word=='W':     return 'Wboson' 
    elif word=='W+':    return 'Wboson' 
    elif word=='W-':    return '-Wboson'
    elif word=='A':     return 'photon'
    elif word=='Z':     return 'Zboson'
    elif word=='g':     return 'gluon'
    elif word=='H':     return 'Higgs'
    elif word=='q':     return 'uQuark' 
    elif word=='q~':    return '-uQuark'
    elif word=='r':     return 'dQuark' 
    elif word=='r~':    return '-dQuark'
    else: return word

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

def frmt_info(key,line):
    if len(line)<90-len(key):
        return ['  write(eventUnit,\'(A)\') "'+key+': '+line+'"\n']
    else:
        return ['  write(eventUnit,\'(A)\') "'+key+': '+line[:88-len(key)]+'"//&\n'] \
              +frmt_char(line[88-len(key):],linelen=120,cont=False)

def frmt(string,lineLen=84,indent=4):
    stringLen = len(string)
    lines = []
    offset = 0
    while indent+stringLen-offset>lineLen:
        temp = string[offset:offset+lineLen-indent-1]
        ii =-1
        while temp[ii] not in [' ','+','-','*','/','(',')','%',',']: ii=ii-1
        lines.append(' '*indent+temp[0:ii]+'&\n')
        offset = offset+lineLen-indent-1+ii
    lines.append(' '*indent+string[offset:stringLen]+'\n')
    return lines

def frmt_list(lines,lineLen=84,indent=4):
    linesOut = []
    for line in lines: linesOut = linesOut+frmt(line,lineLen,indent)
    return linesOut

def _limit(string):
    if string=='<': fout='.gt.'   # indeed, it is the opposite
    elif string=='>': fout='.lt.' # "bail out if not satisfied"
    return fout

def Sdouble(string):
    if   'D' in string: fout=string
    elif 'd' in string: fout=re.sub('d','D',string) 
    elif 'E' in string: fout=re.sub('E','D',string)
    elif 'e' in string: fout=re.sub('e','D',string)
    else: fout = string+'D0'
    return fout

def _observable(string):
    array = re.sub('{','',re.sub('}','' \
                  ,re.sub('{\|','',re.sub('\|}','',string)))).split('|')
    if len(array)==3: array = array+['']
    if len(array)==2: array = array+['','']
    if len(array)==1: array = array+['','','']
    return array

class cuts:

    def __init__(self):
        self.src = []
        self.extracuts = []
        self.extraweights = []
        self.info = []
        self.__evalnList = []
        self.__arrayList = []
        self.__ordrdList = []

    def __hashEvaln(self,func,arg):
        if arg+func in self.__evalnList:
            ii = self.__evalnList.index(arg+func)
        else:
            self.__evalnList.append(arg+func)
            ii = len(self.__evalnList)-1
            line = 'evalnTbl('+str(ii+1)+') = '
            args = '('
            arglist = arg.split(',')
            for jj in range(0,len(arglist)):
                word = arglist[jj].split('+')
                for jj in range(0,len(word)):
                    if re.search('breit',func): part = 'pBreit(:,'
                    else:                       part = 'psPoint%b('
                    if   word[jj]=='electron': args = args+part+'DISlepton)+'
                    elif word[jj]=='positron': args = args+part+'DISlepton)+'
                    elif word[jj]=='lepton':   args = args+part+'DISlepton)+'
                    elif word[jj]=='A': args = args+part+'-1)+'
                    elif word[jj]=='B': args = args+part+'-2)+'
                    else:               args = args+part+word[jj]+')+'
                args = re.sub('\+$',',',args)
            args = re.sub(',$','',args)
            args = args+')'
            if   func=='pT':        line=line+'psPoint%Trans'+args
            elif func=='ET':        line=line+'psPoint%ETrans'+args
            elif func=='energy':    line=line+'psPoint%Energy'+args
            elif func=='rapidity':  line=line+'psPoint%rapidity'+args 
            elif func=='pseudoRap': line=line+'psPoint%pseudoRap'+args
            elif func=='mass':      line=line+'psPoint%Mass'+args
            elif func=='deltaR':    line=line+'psPoint%deltaR'+args
            elif func=='deltaPhi':  line=line+'psPoint%deltaPhi'+args
            elif func=='phi':       line=line+'psPoint%phi'+args
            elif func=='theta':     line=line+'psPoint%theta'+args
            elif func=='plus':      line=line+'psPoint%plus'+args
            elif func=='minus':     line=line+'psPoint%minus'+args
            elif func=='angle':     line=line+'psPoint%angle'+args
            elif func=='Qsquare':   line=line+'Qsquare'
            elif func=='xBjorken':  line=line+'xBjorken'
            elif func=='inelast':   line=line+'yInelast'
            elif func=='pTbreit':   line=line+'pTrans'+args
            elif func=='deltaRbreit':line=line+'deltaR'+args
            else:
                 print('ERROR: cuts for '+func+' not defined')
                 sys.exit()
            self.src.append(line)
        return ii

    def __hashArray(self,func,arg):
        if arg+func in self.__arrayList:
            ii = self.__arrayList.index(arg+func)
        else:
            self.__arrayList.append(arg+func)
            ii = len(self.__arrayList)-1
            splitArg = arg.split(',')
            nn = len(splitArg)
            for jj in range(0,nn):
                kk = self.__hashEvaln(func,splitArg[jj])
                self.src.append('arrayTbl('+str(jj+1)+','+str(ii+1)+')' \
                                       +' = evalnTbl('+str(kk+1)+')')
        return ii

    def __hashOrdrd(self,func,arg):
        if arg+func in self.__ordrdList:
            ii = self.__ordrdList.index(arg+func)
        else:
            self.__ordrdList.append(arg+func)
            ii = len(self.__ordrdList)-1
            splitArg = arg.split(',')
            nn = len(splitArg)
            self.src.append('call order_big2small(ordrdTbl(1:'+str(nn)+','+str(ii+1)+'),' \
                                                 +'arrayTbl(1:'+str(nn)+','+str(ii+1)+'))')
        return ii

    def __fortranExpr(self,line):
        string = ''
        iPrev = 0
        iStart = [i for i, letter in enumerate(line) if letter == '{']
        iEnd = [i for i, letter in enumerate(line) if letter == '}']
        for jj in range(0,len(iStart)):
            argA,argB,argC,argD = _observable(line[iStart[jj]:iEnd[jj]+1])
            if argD=='':
                if argC=='':
                    aa = str(1+self.__hashEvaln(argA,argB))
                    string = string+line[iPrev:iStart[jj]]+'evalnTbl('+aa+')'
                else:
                    aa = str(1+self.__hashArray(argA,argC))
                    bb = str(1+self.__hashOrdrd(argA,argC))
                    string = string+line[iPrev:iStart[jj]]+'arrayTbl(ordrdTbl('+argB+','+bb+'),'+aa+')'
            else:
                aa = str(1+self.__hashArray(argC,argD)) # to make sure this array exists. ii is dummy here
                aa = str(1+self.__hashArray(argA,argD))
                bb = str(1+self.__hashOrdrd(argC,argD))
                string = string+line[iPrev:iStart[jj]]+'arrayTbl(ordrdTbl('+argB+','+bb+'),'+aa+')'
            iPrev = iEnd[jj]+1
        string = string+line[iPrev:]
        return re.sub('REJECT','goto 111',string)

    def __hashExpr(self,expr):
        if expr in self.__evalnList:
            ii = self.__evalnList.index(expr)
        else:
            self.__evalnList.append(expr)
            ii = len(self.__evalnList)-1
            line = 'evalnTbl('+str(ii+1)+') = '+self.__fortranExpr(expr)
            self.src.append(line)
        return ii

    def add_cut(self,line):
        funcarg,lim,val = line.split()
        argA,argB,argC,argD = _observable(funcarg)
        if argD=='':
            if argC=='':
                ii = self.__hashEvaln(argA,argB)
                self.src.append('if (evalnTbl('+str(ii+1)+')' \
                                       +_limit(lim)+Sdouble(val)+') goto 111')
            else: 
                ii = self.__hashArray(argA,argC)
                jj = self.__hashOrdrd(argA,argC)
                self.src.append('if (arrayTbl(ordrdTbl('+argB+','+str(jj+1)+'),'+str(ii+1)+')' \
                                       +_limit(lim)+Sdouble(val)+') goto 111')
        else:
            ii = self.__hashArray(argC,argD) # to make sure this array exists. ii is dummy here
            ii = self.__hashArray(argA,argD)
            jj = self.__hashOrdrd(argC,argD)
            self.src.append('if (arrayTbl(ordrdTbl('+argB+','+str(jj+1)+'),'+str(ii+1)+')' \
                                   +_limit(lim)+Sdouble(val)+') goto 111')

        self.info = self.info+frmt_info('PHASE SPACE CUT',line)

    def add_extra_cut(self,line):
        self.extracuts.append(self.__fortranExpr(line))
        self.info = self.info+frmt_info('PHASE SPACE CUT',line)

    def add_extra_weight(self,line):
        self.extraweights.append('weight = weight*'+self.__fortranExpr(line))
        self.info = self.info+frmt_info('WEIGHT FACTOR',line)

    def set_scale(self,label,line,group='0'):
        if group=='0': lout = label+' = evalnTbl('+str(1+self.__hashExpr(line))+')'
        else:          lout = label+'('+group+') = evalnTbl('+str(1+self.__hashExpr(line))+')'
        for item in self.src:
            if item.startswith(label+' ='): self.src.remove(item)
        for item in self.info:
            if item.startswith('  write(eventUnit,\'(A)\') "HARD SCALE: '+label+':'): self.info.remove(item)
        self.src.append(lout)
        self.info = self.info+frmt_info('HARD SCALE',line)

class parameters:
    src = []
    info = []
    __done = []

    def __init__(self): pass

    def energy(self,line):
        word = line.split()
        if 'Ecm' not in self.__done: self.__done.append('Emc')
        self.src.append('Ecm = '+word[0])

    def ePosRap(self,line):
        word = line.split()
        if 'EposRap' not in self.__done: self.__done.append('EposRap')
        self.src.append('EposRap = '+word[0])

    def eNegRap(self,line):
        word = line.split()
        if 'EnegRap' not in self.__done: self.__done.append('EnegRap')
        self.src.append('EnegRap = '+word[0])

    def esoft(self,line):
        word = line.split()
        if 'Esoft' not in self.__done: self.__done.append('Esoft')
        self.src.append('Esoft = '+Sdouble(word[0]))

    def mass(self,line):
        word = line.split()
        if len(word)==2: word.append('0')
        if word[0][-1] in ['0','+','-','~']: word[0]=word[0][:-1]
        if 'mass_'+word[0] not in self.__done:
            self.__done.append('mass_'+word[0])
        word[1] = Sdouble(word[1])
        word[2] = Sdouble(word[2])
        lout = 'call set_mass_and_width('+particle(word[0]) \
                                         +','+word[1]+','+word[2]+')'
        self.src.append(lout)
        self.info = self.info+frmt_info('MASS AND WIDTH',line)

    def switch(self,line):
        word = line.split()
        if word[0] not in self.__done: self.__done.append(word[0])
        val = '0'
        if word[1] in ['1','ON','On','on','Y','YES','Yes','yes' \
                      ,'T','TRUE','True','true']:
            val = '1'
        self.src.append('call set_'+word[0]+'('+val+')')
        self.info = self.info+frmt_info('INTERACTION SWITCH',line)

    def coupling(self,line):
        word = line.split()
        if word[0] not in self.__done: self.__done.append(word[0])
        if word[0]=='alphaEW':
          lout = 'call set_gQED(sqrt(4*r1PI*'+Sdouble(word[1])+'))'
          self.src.append(lout)
        elif word[0]=='Gfermi':
          if ('mass_W' not in self.__done) or ('mass_Z' not in self.__done):
              print('ERROR in cuts_and_pars.parameters.coupling:') 
              print('      W and Z mass must be set before Gfermi') 
              sys.exit()
          lout = 'call set_Gfermi('+Sdouble(word[1])+')'
          self.src.append(lout)
        elif word[0]=='Higgs-gluon':
          lout = 'call set_gHG('+Sdouble(word[1])+')'
          self.src.append(lout)
        elif word[0]=='Higgs-photon':
          lout = 'call set_gHA('+Sdouble(word[1])+')'
          self.src.append(lout)
        self.info = self.info+frmt_info('COUPLING CONSTANT',line)


