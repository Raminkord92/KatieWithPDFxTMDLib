module read_events_mod
  use katie_histogramtools
  implicit none

  integer,parameter :: Nsize=13
  integer,save,protected :: Nfnst,Ngroups=1
  integer,save,protected :: NfinalState(Nsize),process(Nsize),procNr
  integer,save,protected,allocatable :: NfnstProc(:),lhefInst(:,:),lhefFnst(:,:)
  real(kind(1d0)),save,protected:: Ecm,EposRap,EnegRap,eventWeight
  real(kind(1d0)),save,protected:: xsection,errest
  real(kind(1d0)),save,protected:: pInst(0:3,2,Nsize),sInst(2,Nsize)
  real(kind(1d0)),save,protected:: pFnst(0:3,Nsize),sFnst(Nsize)
  real(kind(1d0)),save,protected:: matrixElement(Nsize),partonLumin(Nsize)
  real(kind(1d0)),save,protected:: alphaStrong(Nsize),muScale(Nsize)
  logical,save,protected :: exitLoop
  character(256),save,protected :: filename

  integer,save,private :: eventUnit=21
  character(256),save,private :: line

contains


  subroutine open_file(eventfileUnit)
  integer,intent(in),optional :: eventfileUnit
  integer :: ii,nn,jj,Nproc,NfinstMax,iProc
  if (present(eventfileUnit)) eventUnit = eventfileUnit
  call get_command_argument(1,filename)
  if (trim(filename).eq.'') then
    write(*,*)
    write(*,*) 'ERROR: you have to give an event file as argument'
    write(*,*)
    stop
  endif
  write(*,*) 'MESSAGE: reading '//trim(filename)
  open(eventUnit,file=trim(filename),status='old')
    ii = 0
    do ;ii = ii+1
      read(eventUnit,'(A)') line
      if (line(1:13).eq.'EVENT WEIGHT:') exit
    enddo
    nn = ii-1
  close(eventUnit)
  open(eventUnit,file=trim(filename),status='old')
    do ii=1,nn
      read(eventUnit,'(A)') line
      if     (line(1:22).eq.'CENTER-OF-MASS ENERGY:') then
        read(line(23:),*) Ecm
      elseif (line(1:30).eq.'NEGATIVE RAPIDITY BEAM ENERGY:') then 
        read(line(31:),*) EnegRap
      elseif (line(1:30).eq.'POSITIVE RAPIDITY BEAM ENERGY:') then 
        read(line(31:),*) EposRap
      elseif (line(1:17).eq.'NUMBER OF GROUPS:') then 
        read(line(18:),*) Ngroups
      elseif (line(1:32).eq.'NUMBER OF FINAL-STATE PARTICLES:') then
        read(line(33:),*) NfinalState(1:Ngroups)
        Nfnst = sum(NfinalState(1:Ngroups))
      elseif (line(1:20).eq.'TOTAL CROSS SECTION:') then
        jj = index(line,'+/-')
        read(line(21:jj-1),*) xsection
        read(line(jj+3:)  ,*) errest        
      elseif (line(1:14).eq.'LHEF ENCODING:') then
        if     (line(16:31).eq.'Nproc,NfinstMax:') then
          read(line(32:),*) Nproc,NfinstMax
          allocate(NfnstProc(         Nproc))
          allocate(lhefInst( 2       ,Nproc))
          allocate(lhefFnst(NfinstMax,Nproc))
          lhefFnst = 0
          iProc = 0
        elseif (line(16:23).eq.'PROCESS:') then
          iProc = iProc+1
          read(line(32:37),*) NfnstProc(iProc)
          read(line(39:51),*) lhefInst(1:2,iProc)
          read(line(53:  ),*) lhefFnst(1:NfnstProc(iProc),iProc)
        endif
      endif
    enddo
  end subroutine


  subroutine read_event
  integer :: kk,ii,jj
  exitLoop = .true.
  call readuntil('EVENT WEIGHT:')
  if (line(1:9).eq.'ENDOFFILE') return
  read(line,*) eventWeight
  kk = 0
  do ii=1,Ngroups
    read(eventUnit,*) process(ii)
    procNr = process(1)
    read(eventUnit,*) pInst(0:3,1,ii),sInst(1,ii)
    read(eventUnit,*) pInst(0:3,2,ii),sInst(2,ii)
    do jj=1,NfinalState(ii)
      kk = kk+1
      read(eventUnit,*) pFnst(0:3,kk),sFnst(kk)
    enddo
    read(eventUnit,*) &
      matrixElement(ii),partonLumin(ii),alphaStrong(ii),muScale(ii)
  enddo
  exitLoop = .false.
!
  contains
!
    subroutine readuntil( pattern )
    character(*),intent(in) :: pattern
    integer :: patlen
    patlen = len(pattern)
    do
      read(eventUnit,'(A)',end=999) line
      if (line(1:patlen).eq.pattern) exit
    enddo
    line(1:) = line(patlen+1:)
    return
    999 continue
    line = 'ENDOFFILE'
    end subroutine
!
  end subroutine


  subroutine close_file
  close(eventUnit)
  end subroutine
  

end module


