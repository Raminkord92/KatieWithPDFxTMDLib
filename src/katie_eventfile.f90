module katie_eventfile
  use katie_sudakov
  implicit none
  private
!
  public :: initialize,read_event,write_event
  public :: Nsize,exitLoop
  public :: eventWeight
  public :: pInstA,flavorA,colorA,anticA,helicityA
  public :: pInstB,flavorB,colorB,anticB,helicityB
  public :: pFinst,flavorF,colorF,anticF,helicityF
  public :: matrixElement,partlumi,alphaS,renScale
  public :: pdfA,xA,kTA,scaleA ,pdfB,xB,kTB,scaleB
  public :: Ecm,EposRap,EnegRap
  public :: Nproc,procNr,Nflavor,Nfinst
  public :: targetDir,label
  public :: partialAmpSqr
!
  integer,parameter :: charLen=256,NeventBlock=7,Nsize=13
!
  integer,protected,save :: flavorF(Nsize),colorF(Nsize),anticF(Nsize),helicityF(Nsize)
  integer,protected,save :: flavorA,colorA,anticA,helicityA,flavorB,colorB,anticB,helicityB
  real(kind(1d0)),protected,save :: pInstA(0:4),pInstB(0:4),pFinst(0:4,Nsize)
  real(kind(1d0)),protected,save :: matrixElement,partlumi,alphaS,renScale
  real(kind(1d0)),protected,save :: pdfA,xA,kTA,scaleA ,pdfB,xB,kTB,scaleB
  real(kind(1d0)),protected,save :: Ecm,EposRap,EnegRap
  integer,protected,save :: exitLoop=0
!
  integer,save :: rawUnit=21,eventUnit=22,iEvent=0
  character(charLen),save :: tmpLine,lastLine,task,lhaSet,targetDir,filename,label
  character(charLen),save :: scatproc
  character(charLen),allocatable,save :: rawFile(:),initBlock(:),lines(:)
  real(kind(1d0)),save :: Nevents,sumWevent,sumW0tot,sumW0pos,sumW1pos,sumW2pos
  real(kind(1d0)),save :: sumW0neg,sumW1neg,sumW2neg,xSecPos,sigSqPos,xSecNeg
  real(kind(1d0)),save :: sigSqNeg,xSect,sigma,sigPos,sigNeg,overallFac,eventWeight
  real(kind(1d0)),save :: x1,x2,x3,x4,x5,x6,x7,x8,x9,alphaEW,sudFac,sudNormFac
  real(kind(1d0)),allocatable,save :: sumW1proc(:),sumW2proc(:),maxWproc(:)
  real(kind(1d0)),allocatable,save :: partialAmpSqr(:)
  integer,save :: NinitBlock,iFile,Nfile,iLine,Nline,ii,jj,kk,iLast
  integer,save :: Ngroup,Nwords,Ntask,Nfinst,NfinstMax,NfinstTot
  integer,save :: iBgn(charLen),iEnd(charLen),posEq
  integer,protected,save :: Nproc,procNr,Nflavor,Npartial
  integer,allocatable,save :: NlineTot(:),Nev(:)
  integer,allocatable,save :: process(:,:),NfinstProc(:)
  logical,save :: nbYES,pbYES,customYES,lhefYES,sudakovYES,individualYES
  logical,save :: eventfileYES
  

contains

  subroutine initialize( rawFile_unit ,eventFile_unit )
  integer,intent(in),optional :: rawFile_unit,eventFile_unit
  if (present(rawFile_unit)) rawUnit = rawFile_unit
  if (present(eventFile_unit)) eventUnit = eventFile_unit
  if (COMMAND_ARGUMENT_COUNT().le.1) then
    write(*,*) ''
    write(*,*) 'USAGE:'
    write(*,*) ''
    write(*,*) 'To create an LHEF using all raw files:'
    write(*,*) '$ ./create_evenfile.out lhef raw*'
    write(*,*) ''
    write(*,*) 'Only using selected raw files:'
    write(*,*) '$ ./create_evenfile.out lhef raw215531.dat raw792327.dat'
    write(*,*) ''
    write(*,*) 'Evenfile in the custom file format:'
    write(*,*) '$ ./create_evenfile.out custom raw*'
    write(*,*) ''
    write(*,*) 'With the cross section in pico barn:'
    write(*,*) '$ ./create_evenfile.out custom,pb raw*'
    write(*,*) ''
    write(*,*) 'Setting the directory to which the file is written:'
    write(*,*) '$ ./create_evenfile.out custom,pb,dir=/tmp raw*'
    write(*,*) ''
    write(*,*) 'Setting the file name:'
    write(*,*) '$ ./create_evenfile.out custom,pb,name=events.dat raw*'
    write(*,*) ''
    write(*,*) 'Setting both:'
    write(*,*) '$ ./create_evenfile.out custom,pb,dir=/tmp,name=events.dat raw*'
    write(*,*) ''
    write(*,*) 'Do not create eventfile:'
    write(*,*) '$ ./create_evenfile.out custom,noEventFile raw*'
    write(*,*) ''
    stop
  endif
  call get_command_argument( 1 ,task )
  call split_string( Ntask ,iBgn ,iEnd ,trim(adjustl(task)) ,separ=',' )
  nbYES=.false.;pbYES=.false.;customYES=.false.;lhefYES=.false.;sudakovYES=.false.
  individualYES=.true.;eventfileYES=.true.
  targetDir = './'
  label = ''
  filename = 'eventfile.dat'
  Npartial = 0
  do ii=1,Ntask
    posEq = index(task(iBgn(ii):iEnd(ii)),'=')
    if     (task(iBgn(ii):iEnd(ii)).eq.'nb')         then;        nbYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'pb')         then;        pbYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'custom')     then;    customYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'merge')      then;    customYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'lhef')       then;      lhefYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'sudakov')    then;   sudakovYES=.true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'combined')   then;individualYES=.false.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'noEventFile')then; eventFileYES=.false.
    elseif (posEq.gt.0) then
       if     (task(iBgn(ii):iBgn(ii)+posEq-2).eq.'dir')   then ;targetDir = task(iBgn(ii)+posEq:iEnd(ii))
       elseif (task(iBgn(ii):iBgn(ii)+posEq-2).eq.'name')  then ; filename = task(iBgn(ii)+posEq:iEnd(ii))
       elseif (task(iBgn(ii):iBgn(ii)+posEq-2).eq.'label') then ;    label = task(iBgn(ii)+posEq:iEnd(ii))
       elseif (task(iBgn(ii):iBgn(ii)+posEq-2).eq.'Npartial') then ;read(task(iBgn(ii)+posEq:iEnd(ii)),*) Npartial
       else
         write(*,*) "ERROR in katie_eventfile: task '"//task(iBgn(ii):iEnd(ii))//"' not defined."
         write(*,*) "      You must give a task before the raw files."
         stop
       endif
    else
      write(*,*) "ERROR in katie_eventfile: task '"//task(iBgn(ii):iEnd(ii))//"' not defined."
      write(*,*) "      You must give a task before the raw files."
      stop
    endif
  enddo
  Nfile = command_argument_count()-1
  allocate(rawFile(Nfile),NlineTot(Nfile),Nev(Nfile))
  do iFile=1,Nfile
    call get_command_argument( 1+iFile ,rawFile(iFile) )
  enddo

  if (Npartial.gt.0) allocate(partialAmpSqr(Npartial))

! open eventfile early to check if it exists already
  if (eventfileYES) then
    open(unit=eventUnit,file=trim(targetDir)//'/'//trim(filename),status='new')
  endif

! Determine last line with accumulated sums for each file
  Nevents=0 ;sumWevent=0 ;sumW0tot=0
  sumW0pos=0 ;sumW1pos=0 ;sumW2pos=0
  sumW0neg=0 ;sumW1neg=0 ;sumW2neg=0
  do iFile=1,Nfile
    open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
    iLine = 0
    do ;iLine=iLine+1
      read(rawUnit,'(A)',end=111) tmpLine
      if (tmpLine(1:17).eq.'ACCUMULATED SUMS:') lastLine = tmpLine
    enddo
    111 continue
    close(rawUnit)
    NlineTot(iFile) = iLine
    read(lastLine(18:),*) x1,x2,x3,x4,x5,x6,x7,x8,x9
    Nev(iFile) = x1
    Nevents   = Nevents   + x1
    sumWevent = sumWevent + x2
    sumW0tot  = sumW0tot  + x3
    sumW0pos  = sumW0pos  + x4
    sumW1pos  = sumW1pos  + x5
    sumW2pos  = sumW2pos  + x6
    sumW0neg  = sumW0neg  + x7
    sumW1neg  = sumW1neg  + x8
    sumW2neg  = sumW2neg  + x9
  enddo  
! Accumulated quantities over all files
  xSecPos    = sumW1pos/sumW0tot
  sigSqPos   = (sumW2pos/sumW0tot-xSecPos**2)/(sumW0tot-1)
  xSecNeg    = sumW1neg/sumW0tot
  sigSqNeg   = (sumW2neg/sumW0tot-xSecNeg**2)/(sumW0tot-1)
  xSect      = xSecPos-xSecNeg
  sigma      = sqrt(sigSqPos+sigSqNeg)
  sigPos     = sqrt(sigSqPos)
  sigNeg     = sqrt(sigSqNeg)
  overallFac = xSect*Nevents/sumWevent
  write(*,'(a26,e15.8,a4,e15.8)') 'Total cross section in nb:',xSect,' +/-',sigma

! Determine NinitBlock
  open(unit=rawUnit,file=trim(rawFile(1)),status='old')
  iLine = 0
  do ;iLine = iLine+1
    read(rawUnit,'(A)') tmpLine
    if (tmpLine(1:13).eq.'EVENT WEIGHT:') exit
  enddo
  close(rawUnit)
  NinitBlock = iLine-1
  allocate(initBlock(NinitBlock),lines(NinitBlock))
  
! Read initBlock
  do iFile=1,Nfile
    open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
    do iLine=1,NinitBlock
      read(rawUnit,'(A)') tmpLine
      if (iFile.eq.1) then
        initBlock(iLine) = tmpLine
      else
        if (tmpLine.ne.initBlock(iLine)) then
          write(*,*) 'ERROR in read_raw: input files not compatible'
          stop
        endif
      endif
    enddo
    close(rawUnit)
  enddo

! Extract information from initBlock
  lines(1) = '1'
  call read_initBlock( lines,Nline ,'NUMBER OF GROUPS:' )
  read(lines(1),*) Ngroup
  if (Ngroup.ne.1) then
    write(*,*) 'ERROR: Ngroup must be equal to 1'
    stop
  endif

  call read_initBlock( lines,Nline ,'SCATTERING PROCESS:' )
  read(lines(1),*) scatproc; scatproc = adjustl(scatproc)

  call read_initBlock( lines,Nline ,'NUMBER OF FINAL-STATE PARTICLES:' )
  read(lines(1),*) Nfinst
  NfinstMax = Nfinst
  NfinstTot = Nfinst

  call read_initBlock( lines,Nline ,'CENTER-OF-MASS ENERGY:' )
  read(lines(1),*) Ecm

  call read_initBlock( lines,Nline ,'NEGATIVE RAPIDITY BEAM ENERGY:' )
  read(lines(1),*) EnegRap

  call read_initBlock( lines,Nline ,'POSITIVE RAPIDITY BEAM ENERGY:' )
  read(lines(1),*) EposRap

  lines(1) = '-1'
  call read_initBlock( lines,Nline ,'ELECTRO-WEAK COUPLING:' )
  read(lines(1),*) alphaEW

  call read_initBlock( lines,Nproc ,'LIST OF PROCESSES: process ' )
  allocate(process(-2:NfinstMax,Nproc),NfinstProc(Nproc))
  allocate(sumW1proc(Nproc),sumW2proc(Nproc),maxWproc(Nproc))
  do iLine=1,Nproc
    iLast = index(lines(iLine),',')-1
    if (iLast.le.0) iLast = index(lines(iLine),'factor')-1
    if (iLast.le.0) iLast = index(lines(iLine),'type')-1
    if (iLast.le.0) iLast = index(lines(iLine),'partlumi')-1
    if (iLast.le.0) iLast = len(lines(iLine))
    call split_string( Nwords ,iBgn ,iEnd ,lines(iLine)(1:iLast) ,separ=' ' )
    NfinstProc(iLine) = Nwords-4 ! This is Nfinst for this process
    process(-2,iLine) = translate( 'lhef' ,lines(iLine)(iBgn(2):iEnd(2)) )
    process(-1,iLine) = translate( 'lhef' ,lines(iLine)(iBgn(3):iEnd(3)) )
    do ii=1,NfinstProc(iLine)
      process(ii,iLine) = translate( 'lhef' ,lines(iLine)(iBgn(4+ii):iEnd(4+ii)) )
    enddo
  enddo

  call read_initBlock( lines,Nline ,'LIST OF PROCESSES: Nf=Nflavor:' )
  read(lines(1),*) Nflavor
  call read_initBlock( lines,Nline ,'LHAPDF SET: ' )
  read(lines(1),*) lhaSet

! Re-allocate lines
  if (NfinstTot+NeventBlock.gt.NinitBlock) then
    deallocate(lines)
    allocate(lines(NfinstTot+NeventBlock))
  endif

! Prepare Sudakov re-weighting
  if (sudakovYES) then
    call sudakov_init(trim(adjustl(lhaSet)),Nflavor)
    x1 = 0
    x2 = 0
    do iFile=1,Nfile
      open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
      do ii=1,NinitBlock
        read(rawUnit,*)
      enddo
      do ii=1,Nev(iFile)
        call read_event_block
        call get_sudFac
        if (sudFac.ne.1d0) then
          x1 = x1 + eventWeight
          x2 = x2 + eventWeight*sudFac
        endif
      enddo
      close(rawUnit)
    enddo
    sudNormFac = x1/x2
  endif
!
  if (customYES.and.pbYES) overallFac = overallFac*1000
  if (lhefYES.and..not.nbYES) overallFac = overallFac*1000
!
  if (eventfileYES) then
  if (lhefYES) then
    write(eventUnit,'(A)') '<LesHouchesEvents version="1.0">'
    write(eventUnit,'(A)') '<header>'
    write(eventUnit,'(A)') '<!-- individually designed XML tags, in fancy XML style -->'
    write(eventUnit,'(A)') '<KaTieInfo>'
    write(eventUnit,'(A)') '<!--'
    open(unit=rawUnit,file=trim(rawFile(1)),status='old')
    do ii=1,NinitBlock
      read(rawUnit,'(A)') tmpLine
      write(eventUnit,'(A)') '# '//trim(tmpLine)
    enddo
    close(rawUnit)
    write(eventUnit,'(A)') '-->'
    write(eventUnit,'(A)') '</KaTieInfo>'
    write(eventUnit,'(A)') '</header>'
    write(eventUnit,'(A)') '<init>'
      if     (trim(scatproc).eq.'DIS-') then
        write(eventUnit,'(a7,2e15.8,a12,i4)') '2212 11',EposRap,EnegRap,' 3 3 41 41 1',Nproc
      elseif (trim(scatproc).eq.'DIS+') then
        write(eventUnit,'(a8,2e15.8,a12,i4)') '2212 -11',EposRap,EnegRap,' 3 3 41 41 1',Nproc
      else
        write(eventUnit,'(a9,2e15.8,a12,i4)') '2212 2212',EposRap,EnegRap,' 3 3 41 41 1',Nproc
      endif
      sumW1proc = 0
      sumW2proc = 0
      maxWproc = 0
      do iFile=1,Nfile
        open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
        do ii=1,NinitBlock
          read(rawUnit,*)
        enddo
        do ii=1,Nev(iFile)
          call read_event_block
          if (sudakovYES) then
            call get_sudFac
            if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
          endif
          eventWeight = eventWeight*overallFac
          sumW1proc(procNr) = sumW1proc(procNr) + eventWeight
          sumW2proc(procNr) = sumW2proc(procNr) + eventWeight**2
          if (eventWeight.gt.maxWproc(procNr)) maxWproc(procNr) = eventWeight
        enddo
        close(rawUnit)
      enddo
      do ii=1,Nproc
        x1 = sumW1proc(ii)/Nevents
        x2 = sumW2proc(ii)/Nevents - x1*x1
        if (x2.gt.0) then
          x2 = sqrt(x2/(Nevents-1))
        else
          x2 = abs(x1)
        endif
        write(eventUnit,'(3e24.16,i4)') x1,x2,maxWproc(ii),ii 
      enddo
    write(eventUnit,'(A)') '</init>'
  endif
  endif
!
  iFile = 1
  iEvent = 1
  exitLoop = 0
  open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
  end subroutine


  subroutine read_event
  if (customYES) then
    if (iEvent.eq.1) then
      if (iFile.eq.1) then
        if (eventfileYES) then
          write(eventUnit,'(a31,2i4)') 'LHEF ENCODING: Nproc,NfinstMax:',Nproc,NfinstMax
          write(eventUnit,'(a69)') 'LHEF ENCODING: LEGEND:   iProc Nfinst  instB  instA finst1 finst2 ...'
          do ii=1,Nproc
            write(eventUnit,'(a23,99i7)') 'LHEF ENCODING: PROCESS:' &
                 ,ii,NfinstProc(ii),process(-2:-1,ii),process(1:NfinstProc(ii),ii)
          enddo
        endif
        do ii=1,NinitBlock
          read(rawUnit,'(A)') tmpLine
          if (eventfileYES) write(eventUnit,'(A)') trim(tmpLine)
        enddo
        if (eventfileYES) then
          write(eventUnit,'(a17,e16.8,a4,e15.8)') 'POSITIVE WEIGHTS:',xSecPos,' +/-',sigPos
          write(eventUnit,'(a17,e16.8,a4,e15.8)') 'NEGATIVE WEIGHTS:',xSecNeg,' +/-',sigNeg
          write(eventUnit,'(a20,e16.8,a4,e15.8)') 'TOTAL CROSS SECTION:',xSect,' +/-',sigma
          write(eventUnit,'(A)') 'GENERAL INFO:  cross section = sum(eventWeight)/Nevents'
          write(eventUnit,'(a17,i12)') 'NUMBER OF EVENTS:',sum(Nev(1:Nfile))
        endif
      else
        do ii=1,NinitBlock
          read(rawUnit,*)
        enddo
      endif
    endif
    call read_event_block
    if (sudakovYES) then
      call get_sudFac
      if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
    endif
    eventWeight = eventWeight*overallFac
!
  elseif (lhefYES) then
    if (iEvent.eq.1) then
      do ii=1,NinitBlock
        read(rawUnit,*)
      enddo
    endif
    call read_event_block
    if (sudakovYES) then
      call get_sudFac
      if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
    endif
    eventWeight = eventWeight*overallFac
    if (helicityA.eq.99) helicityA = 9
    if (helicityB.eq.99) helicityB = 9
    do jj=1,NfinstProc(procNr)
      if (helicityF(jj).eq.99) helicityF(jj) = 9
    enddo
!
  endif
  end subroutine


  subroutine write_event
  if (eventfileYES) then
  if (customYES) then
    write(eventUnit,'(a13,e24.16)') 'EVENT WEIGHT:',eventWeight
    write(eventUnit,'(i3)') procNr
    write(eventUnit,'(5e24.16,2i4,i3)') pInstB(0:4),colorB,anticB,helicityB
    write(eventUnit,'(5e24.16,2i4,i3)') pInstA(0:4),colorA,anticA,helicityA
    do kk=1,NfinstProc(procNr)
      write(eventUnit,'(5e24.16,2i4,i3)') pFinst(0:4,kk),colorF(kk),anticF(kk),helicityF(kk)
    enddo
    write(eventUnit,'(4e24.16)') matrixElement,partlumi,alphaS,renScale
    if (individualYES) then
      write(eventUnit,'(4e24.16)') pdfB,xB,kTB,scaleB
      write(eventUnit,'(4e24.16)') pdfA,xA,kTA,scaleA
    endif
!
  elseif (lhefYES) then
    write(eventUnit,'(A)') '<event>'
    write(eventUnit,'(i3,i4,4e24.16)') &
      NfinstProc(procNr)+2,procNr,eventWeight,renScale,alphaEW,alphaS
    write(eventUnit,'(i4,a7,2i4,5e24.16,a3,i3)') &
      flavorA,' -1 0 0',anticA,colorA,-pInstA(1:3),-pInstA(0),sqrt(abs(pInstA(4))),' 0.',helicityA
    write(eventUnit,'(i4,a7,2i4,5e24.16,a3,i3)') &
      flavorB,' -1 0 0',anticB,colorB,-pInstB(1:3),-pInstB(0),sqrt(abs(pInstB(4))),' 0.',helicityB
    do jj=1,NfinstProc(procNr)
      write(eventUnit,'(i4,a7,2i4,5e24.16,a3,i3)') &
        flavorF(jj),' +1 1 2',colorF(jj),anticF(jj) &
        ,pFinst(1:3,jj),pFinst(0,jj),sqrt(abs(pFinst(4,jj))),' 0.',helicityF(jj)
    enddo
    write(eventUnit,'(a11,e24.16)') '  #pdf1pdf2',partlumi
    if (individualYES) then
      write(eventUnit,'(a13,4e24.16)') '  #pdf_posRap',pdfA,xA,kTA,scaleA
      write(eventUnit,'(a13,4e24.16)') '  #pdf_negRap',pdfB,xB,kTB,scaleB
    endif
    write(eventUnit,'(A)') '</event>'
    if (iFile.eq.Nfile.and.iEvent.eq.Nev(Nfile)) then
      write(eventUnit,'(A)') '</LesHouchesEvents>'
    endif
!
  endif
  endif
!
  iEvent = iEvent+1
  if (iEvent.gt.Nev(iFile)) then
    iFile = iFile+1
    if (iFile.gt.Nfile) then
      close(rawUnit)
      if (eventfileYES) close(eventUnit)
      exitLoop = 1
    else
      close(rawUnit)
      iEvent = 1
      exitLoop = 0
      open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
    endif
  else
    exitLoop = 0
  endif
!
  end subroutine



  subroutine read_initBlock( lines,Nline ,pattern ) 
  character(*)      ,intent(in) :: pattern
  integer           ,intent(out) :: Nline
  character(charLen),intent(out) :: lines(NinitBlock)
  integer :: iLine,lenPattern
  Nline = 0
  lenPattern = len(pattern)
  do iLine=1,NinitBlock
    if (initBlock(iLine)(1:lenPattern).eq.pattern) then
      Nline = Nline+1
      lines(Nline)(1:) = initBlock(iLine)(lenPattern+1:)
    endif
  enddo
  end subroutine

  subroutine split_string( Nwords,iBgn,iEnd ,line ,separ )
  character(*),intent(in) :: line
  character(1),intent(in) :: separ
  integer,intent(out) :: Nwords,iBgn(charLen),iEnd(charLen)
  integer :: ii,lineLen
  lineLen = len(line)
  Nwords = 0
  ii = 0
  do
    do ;ii=ii+1 ;if(ii.gt.lineLen)exit 
      if(line(ii:ii).ne.separ)exit
    enddo
    if (ii.gt.lineLen) exit
    Nwords = Nwords+1
    iBgn(Nwords) = ii
    do ;ii=ii+1 ;if(ii.gt.lineLen)exit 
      if(line(ii:ii).eq.separ)exit
    enddo
    iEnd(Nwords) = ii-1
    if (ii.gt.lineLen) exit
  enddo
  end subroutine

  subroutine read_event_block
  read(rawUnit,'(A)') tmpLine
  read(tmpLine(14:),*) eventWeight
  read(rawUnit,*) procNr
  flavorB = process(-2,procNr)
  flavorA = process(-1,procNr)
  flavorF(1:NfinstProc(procNr)) = process(1:NfinstProc(procNr),procNr)
  read(rawUnit,*) pInstB(0:4),colorB,anticB,helicityB
  read(rawUnit,*) pInstA(0:4),colorA,anticA,helicityA
  do kk=1,NfinstProc(procNr)
    read(rawUnit,*) pFinst(0:4,kk),colorF(kk),anticF(kk),helicityF(kk)
  enddo
  read(rawUnit,*) matrixElement,partlumi,alphaS,renScale
  if (individualYES) then
    read(rawUnit,*) pdfB,xB,kTB,scaleB
    read(rawUnit,*) pdfA,xA,kTA,scaleA
  endif
  if (Npartial.gt.0) read(rawUnit,*) partialAmpSqr(1:Npartial)
  read(rawUnit,*)
  end subroutine

  subroutine get_sudFac
  real(kind(1d0)) :: kT
  integer :: jj
  sudFac = 1d0
  if (pInstA(4).ne.0d0) then
    kT = sqrt(abs(pInstA(4)))
    if (process(-1,procNr).eq.21) then
      sudFac = sudFac * sudakov_g(kT,scaleA)
    else
      sudFac = sudFac * sudakov_q(kT,scaleA)
    endif
  endif
  if (pInstB(4).ne.0d0) then
    kT = sqrt(abs(pInstB(4)))
    if (process(-2,procNr).eq.21) then
      sudFac = sudFac * sudakov_g(kT,scaleB)
    else
      sudFac = sudFac * sudakov_q(kT,scaleB)
    endif
  endif
  end subroutine


  function translate(option,string) result(encoding)
  character(*),intent(in) :: option,string
  integer :: encoding
  select case (option)
  case default 
    encoding = 0
  case ('lhef')
    select case (string)
    case default  ;encoding = 0
    case('ve'   ) ;encoding = 12
    case('ve~'  ) ;encoding =-12
    case('e-'   ) ;encoding = 11
    case('e+'   ) ;encoding =-11
    case('u'    ) ;encoding = 2
    case('u~'   ) ;encoding =-2
    case('d'    ) ;encoding = 1
    case('d~'   ) ;encoding =-1
    case('vmu'  ) ;encoding = 14
    case('vmu~' ) ;encoding =-14
    case('mu-'  ) ;encoding = 13
    case('mu+'  ) ;encoding =-13
    case('c'    ) ;encoding = 4
    case('c~'   ) ;encoding =-4
    case('s'    ) ;encoding = 3
    case('s~'   ) ;encoding =-3
    case('vtau' ) ;encoding = 16
    case('vtau~') ;encoding =-16
    case('tau'  ) ;encoding = 15
    case('tau~' ) ;encoding =-15
    case('t'    ) ;encoding = 6
    case('t~'   ) ;encoding =-6
    case('b'    ) ;encoding = 5
    case('b~'   ) ;encoding =-5
    case('W+'   ) ;encoding = 24
    case('W-'   ) ;encoding =-24
    case('A'    ) ;encoding = 22
    case('Z'    ) ;encoding = 23
    case('g'    ) ;encoding = 21
    case('H'    ) ;encoding = 37
    case('q' ) ;encoding = 2
    case('q~') ;encoding =-2
    case('r' ) ;encoding = 1
    case('r~') ;encoding =-1
    end select
  end select
  end function


 end module

