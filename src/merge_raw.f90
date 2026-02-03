program read_raw
  use katie_sudakov
  implicit none
  integer,parameter :: charLen=256,rawUnit=21,eventUnit=22
  character(charLen) :: tmpLine,lastLine,task,lhaSet,scatproc
  character(charLen),allocatable :: rawFile(:),initBlock(:),lines(:)
  real(kind(1d0)) :: Nevents,sumWevent,sumW0tot,sumW0pos,sumW1pos,sumW2pos
  real(kind(1d0)) :: sumW0neg,sumW1neg,sumW2neg,xSecPos,sigSqPos,xSecNeg
  real(kind(1d0)) :: sigSqNeg,xSect,sigma,sigPos,sigNeg,overallFac,eventWeight
  real(kind(1d0)) :: x1,x2,x3,x4,x5,x6,x7,x8,x9,alphaEW,sudFac,sudNormFac
  real(kind(1d0)) :: Ecm,EposRap,EnegRap
  real(kind(1d0)),allocatable :: sumW1proc(:),sumW2proc(:),maxWproc(:)
  real(kind(1d0)),allocatable :: pp(:,:,:)
  real(kind(1d0)) :: pdfA,xA,kTA,scaleA ,pdfB,xB,kTB,scaleB
  integer :: NinitBlock,iFile,Nfile,iLine,Nline,ii,iProc(13),Nproc,jj,kk,iLast
  integer :: Ngroup,Nwords,Ntask,Nfinst(13),NfinstMax,NfinstTot,Nflavor,cntr,nn
  integer :: posEq,IDWTUP
  integer :: iBgn(charLen),iEnd(charLen)
  integer,allocatable :: process(:,:),col1(:,:),col2(:,:),hel(:,:)
  integer,allocatable :: NlineTot(:),Nev(:)
  real(kind(1d0)),allocatable :: matrixElement(:),partlumi(:),alphaS(:),muScale(:)
  logical :: nbYES,pbYES,mergeYES,lhefYES,sudakovYES,individualYES
  
  call get_command_argument( 1 ,task )
  call split_string( Ntask ,iBgn ,iEnd ,trim(adjustl(task)) ,separ=',' )
  nbYES=.false.;pbYES=.false.;mergeYES=.false.;lhefYES=.false.;sudakovYES=.false.
  individualYES=.true.;IDWTUP=1
  do ii=1,Ntask
    posEq = index(task(iBgn(ii):iEnd(ii)),'=')
    if     (task(iBgn(ii):iEnd(ii)).eq.'nb')       then;         nbYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'pb')       then;         pbYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'merge')    then;      mergeYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'custom')   then;      mergeYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'lhef')     then;       lhefYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'sudakov')  then;    sudakovYES = .true.
    elseif (task(iBgn(ii):iEnd(ii)).eq.'combined') then; individualYES = .false.
    elseif (posEq.gt.0) then
       if (task(iBgn(ii):iBgn(ii)+posEq-2).eq.'IDWTUP') then ;read(task(iBgn(ii)+posEq:iEnd(ii)),*) IDWTUP
       else
         write(*,*) "ERROR in merge_raw: task '"//task(iBgn(ii):iEnd(ii))//"' not defined."
         write(*,*) "      You must give a task before the raw files."
         stop
       endif
    else
      write(*,*) "ERROR in merge_raw: task '"//task(iBgn(ii):iEnd(ii))//"' not defined."
      write(*,*) "      You must give a task before the raw files."
      stop
    endif
  enddo
  Nfile = command_argument_count()-1
  allocate(rawFile(Nfile),NlineTot(Nfile),Nev(Nfile))
  do iFile=1,Nfile
    call get_command_argument( 1+iFile ,rawFile(iFile) )
  enddo

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
  call read_initBlock( lines,Nline ,'SCATTERING PROCESS:' )
  read(lines(1),*) scatproc; scatproc = adjustl(scatproc)

  lines(1) = '1'
  call read_initBlock( lines,Nline ,'NUMBER OF GROUPS:' )
  read(lines(1),*) Ngroup

  call read_initBlock( lines,Nline ,'NUMBER OF FINAL-STATE PARTICLES:' )
  read(lines(1),*) Nfinst(1:Ngroup)
  NfinstMax = maxval(Nfinst(1:Ngroup))
  NfinstTot = sum(Nfinst(1:Ngroup))
  allocate(pp(0:4,2+NfinstMax,Ngroup))
  allocate(col1(2+NfinstMax,Ngroup),col2(2+NfinstMax,Ngroup),hel(2+NfinstMax,Ngroup))
  allocate(matrixElement(Ngroup),partlumi(Ngroup),alphaS(Ngroup),muScale(Ngroup))

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
  allocate(process(0:2+NfinstMax,Nproc))
  allocate(sumW1proc(Nproc),sumW2proc(Nproc),maxWproc(Nproc))
  do iLine=1,Nproc
    iLast = index(lines(iLine),',')-1
    if (iLast.le.0) iLast = index(lines(iLine),'factor')-1
    if (iLast.le.0) iLast = index(lines(iLine),'type')-1
    if (iLast.le.0) iLast = len(lines(iLine))
    call split_string( Nwords ,iBgn ,iEnd ,lines(iLine)(1:iLast) ,separ=' ' )
    process(0,iLine) = Nwords-4 ! This is Nfinst for this process
    process(1,iLine) = translate( 'lhef' ,lines(iLine)(iBgn(2):iEnd(2)) )
    process(2,iLine) = translate( 'lhef' ,lines(iLine)(iBgn(3):iEnd(3)) )
    do ii=1,process(0,iLine)
      process(2+ii,iLine) &
      = translate( 'lhef' ,lines(iLine)(iBgn(4+ii):iEnd(4+ii)) )
    enddo
  enddo

  call read_initBlock( lines,Nline ,'LIST OF PROCESSES: Nf=Nflavor:' )
  read(lines(1),*) Nflavor
  call read_initBlock( lines,Nline ,'LHAPDF SET: ' )
  read(lines(1),*) lhaSet

! Re-allocate lines
  if (NfinstTot+Ngroup*4+1.gt.NinitBlock) then
    deallocate(lines)
    allocate(lines(NfinstTot+Ngroup*4+1))
  endif

! Prepare Sudakov re-weighting
  if (sudakovYES) then
    if (Ngroup.gt.1) then
      write(*,*) 'ERROR: Sudakov re-weighting not implemented for MPI'
      stop
    endif
    call sudakov_init(trim(adjustl(lhaSet)),Nflavor)
    x1 = 0
    x2 = 0
    do iFile=1,Nfile
      open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
      do ii=1,NinitBlock
        read(rawUnit,*)
      enddo
      do ii=1,Nev(iFile)
        call read_event
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

! Write events to event file
  if (mergeYES) then
!!! My format
    if (pbYES) overallFac = overallFac*1000
    open(unit=eventUnit,file='eventfile.dat',status='new')
    do iFile=1,Nfile
      open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
      if (iFile.eq.1) then
        write(eventUnit,'(a31,2i4)') 'LHEF ENCODING: Nproc,NfinstMax:',Nproc,NfinstMax
        write(eventUnit,'(a69)') 'LHEF ENCODING: LEGEND:   iProc Nfinst  inst1  inst2 finst1 finst2 ...'
        do ii=1,Nproc
          write(eventUnit,'(a23,99i7)') 'LHEF ENCODING: PROCESS:',ii,process(0:2+process(0,ii),ii)
        enddo
        do ii=1,NinitBlock
          read(rawUnit,'(A)') tmpLine
          write(eventUnit,'(A)') trim(tmpLine)
        enddo
        write(eventUnit,'(a17,e16.8,a4,e15.8)') 'POSITIVE WEIGHTS:',xSecPos,' +/-',sigPos
        write(eventUnit,'(a17,e16.8,a4,e15.8)') 'NEGATIVE WEIGHTS:',xSecNeg,' +/-',sigNeg
        write(eventUnit,'(a20,e16.8,a4,e15.8)') 'TOTAL CROSS SECTION:',xSect,' +/-',sigma
        write(eventUnit,'(A)') 'GENERAL INFO:  cross section = sum(eventWeight)/Nevents'
        write(eventUnit,'(a17,i12)') 'NUMBER OF EVENTS:',sum(Nev(1:Nfile))
      else
        do ii=1,NinitBlock
          read(rawUnit,*)
        enddo
      endif
      do ii=1,Nev(iFile)
        call read_event
        if (sudakovYES) then
          call get_sudFac
          if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
        endif
        eventWeight = eventWeight*overallFac
        write(eventUnit,'(a13,e24.16)') 'EVENT WEIGHT:',eventWeight
        do jj=1,Ngroup
          nn = (jj-1)*10
          write(eventUnit,'(i3)') iProc(jj)
          do kk=1,process(0,iProc(jj))+2
            write(eventUnit,'(5e24.16,2i4,i3)') pp(0:4,kk,jj),col1(kk,jj)+nn,col2(kk,jj)+nn,hel(kk,jj)
          enddo
          write(eventUnit,'(4e24.16)') matrixElement(jj),partlumi(jj),alphaS(jj),muScale(jj)
        enddo
      enddo
      close(rawUnit)
    enddo
    close(eventUnit)

  elseif (lhefYES) then
!!! Les Houches format
    !if (Ngroup.gt.1) then
    !  write(*,*) 'ERROR: lhef format not available for Ngroup>1.'
    !  stop
    !endif
    overallFac = overallFac*1000
    if (nbYES) overallFac = overallFac/1000
    open(unit=eventUnit,file='eventfile.dat',status='new')
    write(eventUnit,'(A)') '<LesHouchesEvents version="1.0">'
    write(eventUnit,'(A)') '<header>'
    write(eventUnit,'(A)') '<!-- individually designed XML tags, in fancy XML style -->'
    write(eventUnit,'(A)') '<KaTieInfo>'
    open(unit=rawUnit,file=trim(rawFile(1)),status='old')
    do ii=1,NinitBlock
      read(rawUnit,'(A)') tmpLine
      write(eventUnit,'(A)') trim(tmpLine)
    enddo
    close(rawUnit)
    write(eventUnit,'(A)') '</KaTieInfo>'
    write(eventUnit,'(A)') '</header>'
    write(eventUnit,'(A)') '<init>'
      if (Ngroup.eq.1) then
        if     (trim(scatproc).eq.'DIS-') then
          write(eventUnit,'(a7,2e15.8,a10,i2,i4)') '2212 11',EposRap,EnegRap,' 3 3 41 41',IDWTUP,Nproc
        elseif (trim(scatproc).eq.'DIS+') then
          write(eventUnit,'(a8,2e15.8,a10,i2,i4)') '2212 -11',EposRap,EnegRap,' 3 3 41 41',IDWTUP,Nproc
        else
          write(eventUnit,'(a9,2e15.8,a10,i2,i4)') '2212 2212',EposRap,EnegRap,' 3 3 41 41',IDWTUP,Nproc
        endif
      else
        if     (trim(scatproc).eq.'DIS-') then
          write(eventUnit,'(a7,2e15.8,a10,i2,a2)') '2212 11',EposRap,EnegRap,' 3 3 41 41',IDWTUP,' 1'
        elseif (trim(scatproc).eq.'DIS+') then
          write(eventUnit,'(a8,2e15.8,a10,i2,a2)') '2212 -11',EposRap,EnegRap,' 3 3 41 41',IDWTUP,' 1'
        else
          write(eventUnit,'(a9,2e15.8,a10,i2,a2)') '2212 2212',EposRap,EnegRap,' 3 3 41 41',IDWTUP,' 1'
        endif
      endif
      sumW1proc = 0
      sumW2proc = 0
      maxWproc = 0
      if (Ngroup.eq.1) then
        do iFile=1,Nfile
          open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
          do ii=1,NinitBlock
            read(rawUnit,*)
          enddo
          do ii=1,Nev(iFile)
            call read_event
            if (sudakovYES) then
              call get_sudFac
              if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
            endif
            eventWeight = eventWeight*overallFac
            sumW1proc(iProc(1)) = sumW1proc(iProc(1)) + eventWeight
            sumW2proc(iProc(1)) = sumW2proc(iProc(1)) + eventWeight**2
            if (eventWeight.gt.maxWproc(iProc(1))) maxWproc(iProc(1)) = eventWeight
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
      else
        do iFile=1,Nfile
          open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
          do ii=1,NinitBlock
            read(rawUnit,*)
          enddo
          do ii=1,Nev(iFile)
            call read_event
            eventWeight = eventWeight*overallFac
            if (eventWeight.gt.maxWproc(1)) maxWproc(1) = eventWeight
          enddo
          close(rawUnit)
        enddo
        write(eventUnit,'(3e24.16,1x,a1)') xSect,sigma,maxWproc(1),'1' 
      endif
    write(eventUnit,'(A)') '</init>'
    do iFile=1,Nfile
      open(unit=rawUnit,file=trim(rawFile(iFile)),status='old')
      do ii=1,NinitBlock
        read(rawUnit,*)
      enddo
      do ii=1,Nev(iFile)
        call read_event
        if (sudakovYES) then
          call get_sudFac
          if (sudFac.ne.1d0) eventWeight = eventWeight*sudFac*sudNormFac
        endif
        eventWeight = eventWeight*overallFac
        write(eventUnit,'(A)') '<event>'
        if (Ngroup.eq.1) then
          write(eventUnit,'(i3,i4,4e24.16)') &
            Nfinst(1)+2,iProc(1),eventWeight,muScale(1),alphaEW,alphaS(1)
        else
          write(eventUnit,'(i3,1x,a1,4e24.16)') &
            sum(Nfinst(1:Ngroup))+2*Ngroup,'1',eventWeight,muScale(1),alphaEW,alphaS(1)
        endif
        cntr = 0
        do jj=1,Ngroup
          nn = (jj-1)*10
          do kk=1,process(0,iProc(jj))+2
            if (hel(kk,jj).eq.99) hel(kk,jj) = 9
          enddo
          cntr = cntr+1
          kk = 2
          write(eventUnit,'(i4,a9,2i4,5e24.16,a3,i3)') process(kk,iProc(jj)),' -1  0  0' &
          ,col2(kk,jj)+nn,col1(kk,jj)+nn,-pp(1:3,kk,jj),-pp(0,kk,jj),sqrt(abs(pp(4,kk,jj))),' 0.',hel(kk,jj)
          cntr = cntr+1
          kk = 1
          write(eventUnit,'(i4,a9,2i4,5e24.16,a3,i3)') process(kk,iProc(jj)),' -1  0  0' &
          ,col2(kk,jj)+nn,col1(kk,jj)+nn,-pp(1:3,kk,jj),-pp(0,kk,jj),sqrt(abs(pp(4,kk,jj))),' 0.',hel(kk,jj)
          do kk=3,process(0,iProc(jj))+2
            cntr = cntr+1
            write(eventUnit,'(i4,a3,2i3,2i4,5e24.16,a3,i3)') process(kk,iProc(jj)),'  1',cntr-kk+1,cntr-kk+2 &
            ,col1(kk,jj)+nn,col2(kk,jj)+nn,pp(1:3,kk,jj),pp(0,kk,jj),sqrt(abs(pp(4,kk,jj))),' 0.',hel(kk,jj)
          enddo
        enddo
        write(eventUnit,'(a12,99e24.16)') '  #pdf1pdf2',partlumi(1:Ngroup)
        write(eventUnit,'(A)') '</event>'
      enddo
      close(rawUnit)
    enddo
    write(eventUnit,'(A)') '</LesHouchesEvents>'
    close(eventUnit)

  endif

contains

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

  subroutine read_event
  read(rawUnit,'(A)') tmpLine
  read(tmpLine(14:),*) eventWeight
  do jj=1,Ngroup
    read(rawUnit,*) iProc(jj)
    do kk=1,process(0,iProc(jj))+2
      read(rawUnit,*) pp(0:4,kk,jj),col1(kk,jj),col2(kk,jj),hel(kk,jj)
    enddo
    read(rawUnit,*) matrixElement(jj),partlumi(jj),alphaS(jj),muScale(jj)
    if (individualYES) then
      read(rawUnit,*) pdfB,xB,kTB,scaleB
      read(rawUnit,*) pdfA,xA,kTA,scaleA
    endif
  enddo
  read(rawUnit,*)
  end subroutine

  subroutine get_sudFac
  real(kind(1d0)) :: kT
  integer :: jj
  sudFac = 1d0
  do jj=1,2
    if (pp(4,jj,1).ne.0d0) then
      kT = sqrt(abs(pp(4,jj,1)))
      if (process(1,iProc(1)).eq.21) then
        sudFac = sudFac * sudakov_g(kT,muScale(1))
      else
        sudFac = sudFac * sudakov_q(kT,muScale(1))
      endif
    endif
  enddo
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

end program
