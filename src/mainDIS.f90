program main
  use avh_iounits
  use avh_iotools
  use avh_doloops
  use avh_prnt
  use avh_mathcnst
  use avh_random
  use avh_lagrangian
  use avh_coloredlgn
  use avh_kinematics
  use avh_kaleu_inst
  use avh_kaleu_base
  use katie_matrixelement
  use katie_DISinst
  use katie_partlumi
  use katie_pdfs  !{itmdf=no!}
  use katie_itmds !{itmdf=yes!}
  use katie_model_kaleu_qcd
  use katie_histogramtools, only: breit_type,deltaR,pTrans
  use avh_mctools
  use avh_othercnst
  use avh_lorentz, only: init_lorentz
  use avh_ranvar, only: init_ranvar
  use katie_events
  use katie_version
  implicit none

  integer,parameter :: NprocTot = !(NprocTot!)
  integer,parameter :: Nfinst = !(NfinstPureDIS!)

  type(psPoint_type) :: psp_q,psPoint
  integer :: Iev
  
  real(kind(1d0)) :: Ecm,Esoft,EposRap,EnegRap
  real(kind(1d0)) :: renScale,scaleA,scaleB

  real(kind(1d0)),save :: sumProcWght(0:NprocTot)

  abstract interface
    function pdfsum_interface( pdfB ,pdfA ) result(rslt)
    real(kind(1d0)),intent(in) :: pdfB(-6:6),pdfA(-6:6)
    real(kind(1d0)) :: rslt
    end function
  end interface

  type :: procVal_type
    integer :: pQCD,pNonQCD(3),parton(-2:-1)
    type(process_type) :: flav,psflav
    type(matrixelement_type) :: me
    real(kind(1d0)) :: factor,mXtrn(-2:9),sXtrn(-2:-1)
    procedure(pdfsum_interface),pointer,nopass :: pdf
    type(kaleu_base_type) :: kaleu
    type(DISinst_type) :: inst 
    character(256) :: label
  end type

  integer,parameter :: allInc(-2:10,10)=reshape( &
    [ 3,1 ,0 ,2,0,0,0,0,0,0,0, 0,0 , 4,1 ,0 ,2,3,0,0,0,0,0,0 ,0, 0 &
    , 5,1 ,0 ,2,3,4,0,0,0,0,0, 0,0 , 6,1 ,0 ,2,3,4,5,0,0,0,0 ,0, 0 &
    , 7,1 ,0 ,2,3,4,5,6,0,0,0, 0,0 , 8,1 ,0 ,2,3,4,5,6,7,0,0 ,0, 0 &
    , 9,1 ,0 ,2,3,4,5,6,7,8,0, 0,0 ,10,1 ,0 ,2,3,4,5,6,7,8,9 ,0, 0 &
    ,11,1 ,0 ,2,3,4,5,6,7,8,9,10,0 ,12,1 ,0 ,2,3,4,5,6,7,8,9,10,11 ] ,[13,10] )

  integer,parameter :: fivehund(0:12)=500*[0,1,1,1,1,1,1,1,1,1,1,1,1]

  integer,parameter :: ostask(-2:-1) = !(offshell!)
  integer,parameter :: Nflavors = !(Nflavors!)
  real(kind(1d0)) &
    :: wghtCnst,weight,rho,kTsq,crudeXsec,weightFac,fluxFac,wght_inst,alphaS &
      ,valPartLumi,matrixelement,ww,pInst(0:3,-2:-1),xA,qMOM(0:3),kMOM(0:3) &
      ,dsrPrecision,dsrNevents,rcc,rhl,evalnTbl(99),arrayTbl(11,9) &
      ,Qsquare,yInelast,xBjorken,Qvirtual(0:3),pBreit(0:3,13)
!{itmdf=yes
  real(kind(1d0)) :: tmdListA(10)
!}itmdf
  integer &
    :: Iop,ii,jj,iProc,seed,jGroup &
      ,eventUnit,colcon(2,Nfinst+2),icc(2),ordrdTbl(11,9) &
      ,helicity(Nfinst+2),writeHelicity(-2:Nfinst)
  integer :: kinID(1) = 0
  logical &
    :: discard,accept,stopMC
  character(maxLineLen) :: prefix,eventFile
  type(procVal_type) :: proc(NprocTot)
  type(kaleu_model_type) :: kaleuModel
  type(commandArgs_type) :: comarg
  integer,parameter :: DISlepton = Nfinst
  type(breit_type) :: breit

  call set_unit('banner',-1)

  call init_mathcnst
  call init_ranvar
  call init_lorentz

  call version

  call comarg%grab
  if (comarg%Nargs().lt.1) then
    write(*,*)
    write(*,*)'Execute as, for example:'
    write(*,*)
    write(*,*)'$ ./main.out seed=12345'
    write(*,*) 'or'
    write(*,*)'$ nohup ./main.out seed=12345 > output12345 &'
    write(*,*) 'or'
    write(*,*)'$ nohup ./main.out seed=12345 dir=R001/ > R001/output &'
    write(*,*)
    write(*,*)'The Monte Carlo will by default run until a precision of 0.001 (0.1%)'
    write(*,*)'is reached, and will try to create a raw event file with about 10^5'
    write(*,*)'weighted events. The rule is: the more events, the more their weights'
    write(*,*)'fluctuate. You can change these defaults with optional arguments, eg.'
    write(*,*)'  precision=0.01  Nevents=1e6'
    write(*,*)
    stop
  endif
  eventFile = comarg%find_char('seed=')
    if (eventFile.eq.'0') eventFile='54321'
  read(eventFile,*) seed
  eventFile = 'raw'//trim(eventFile)//'.dat'
  prefix = comarg%find_char('dir=')
    if (prefix.eq.'0') prefix='./'
  dsrPrecision = comarg%find_real('precision=')
    if (dsrPrecision.le.0) dsrPrecision=1d-3
  dsrNevents = comarg%find_real('Nevents=')
    if (dsrNevents.le.0) dsrNevents=1d5
  write(messu,*) 'MESSAGE: seed =',seed
  write(messu,*) 'MESSAGE: dir = ',trim(prefix)
  write(messu,*) 'MESSAGE: raw event file = ',trim(eventFile)
  write(messu,*) 'MESSAGE: desired precision = ',dsrPrecision
  write(messu,*) 'MESSAGE: desired Nevents = ',dsrNevents

  call pdfs_init_base(Nflavors)
!{withLHAPDF=yes
  call pdfs_init_lhapdf('!(lhaSet!)')
!}withLHAPDF
!{withTMDlib0=yes
  call pdfs_init_tmdlib('!(TMDlibSetA!)') !{TMDlibKey=char!}
  call pdfs_init_tmdlib(!(TMDlibSetA!))   !{TMDlibKey=int!}
!}withTMDlib0
!{withPDFxTMDCoupling=yes
  call pdfs_init_pdfxtmd_coupling('!(PDFxTMDCouplingSet!)')
!}withPDFxTMDCoupling
!{withPDFxTMD0=yes
  call pdfs_init_pdfxtmd('!(PDFxTMDSetA!)',!(PDFxTMDMemberA!))
!}withPDFxTMD0
!{withPDFxTMDCPDFA=yes
  call pdfs_init_pdfxtmd_cpdf('!(PDFxTMDSetA!)',!(PDFxTMDMemberA!))
!}withPDFxTMDCPDFA
![tmds
!]tmds
![pars
!]pars

![kinematics
  call add_kinematics( kinID(1) ,5 ,[1,0,0,0,1] )
!]kinematics

![processes
  proc(1)%pNonQCD = [0,0,0]
  call proc(1)%flav%put([eleon,uQuark],[uQuark,gluon,eleon],get_anti)
  call proc(1)%psflav%put([photon,uQuark],[uQuark,gluon],get_anti)
  proc(1)%factor = 1
  proc(1)%label = 'someLabel' 
  proc(1)%pdf => pdf_g_g
!]processes

  sumProcWght(0) = 0
  do ii=1,NprocTot
    associate( prc=>proc(ii) )
    prc%factor = prc%factor &
               / get_symFac(Nfinst,prc%flav%infi(1:Nfinst)) &
               / get_NcolDof(prc%flav%infi(-1))
    call prc%me%init( kinID(1) ,prc%flav%xtrn ,prc%pNonQCD &
                     ,helType='sum' &    !{helicity=sum!}
                     ,helType='no_sum' & !{helicity!=sum!}
!{itmdf=yes
                     ,itmdf=.true. &
                     ,leadingColor=.true. & !{leadingColor=yes!}
!}itmdf
                    )
    prc%pQCD = Nfinst-sum(prc%pNonQCD)+2*prc%pNonQCD(2)
    do jj=-2,Nfinst ;if(jj.eq.0)cycle
      prc%mXtrn(jj) = mass(abs(prc%flav%infi(jj)))
      if (jj.lt.0) prc%sXtrn(jj) = prc%mXtrn(jj)**2
    enddo
    open(99,file=trim(prc%label)//'result.dat',status='old')
    read(99,*) ww
    close(99)
    sumProcWght(ii) = sumProcWght(ii-1)+ww
    end associate
  enddo
  crudeXsec = sumProcWght(NprocTot)
  sumProcWght(1:NprocTot) = sumProcWght(1:NprocTot)/crudeXsec

  call rangen_init(seed)

  wghtCnst = 389379.66d0 * (r2PI)**(4-3*Nfinst) / (2*2)
  if (ostask(-1).eq.1) wghtCnst = wghtCnst/r1PI
  call psp_q%alloc(Nfinst-1)
  call psPoint%alloc(Nfinst)

  call model_kaleu( kaleuModel )
  call kaleu_set_unit( 'all' ,-1 )
  call kaleu_set_unit( 'progress' ,6 )

  do ii=1,NprocTot
    associate( prc=>proc(ii) )
    call prc%inst%init( ranfun ,ostask(-1) ,EposRap,EnegRap ,Nfinst-1 &
                       ,xBmin=!(xBmin!) ,xBmax=!(xBmax!) &
                       ,QsqMin=!(QsqMin!) ,QsqMax=!(QsqMax!) &
                       ,pathFile=prc%label ,readUnit=99 )
    call prc%kaleu%init( kaleuModel ,prc%psflav%apply(flavor_kaleu) ,Ecm,Esoft ,ranfun,ranfun &
                        ,fileName=prc%label ,readUnit=99 )
    end associate
  enddo

  eventUnit=98;open(eventUnit,file=trim(prefix)//trim(eventFile),status='new')
  !This is f2008: open(newunit=eventUnit,file=trim(prefix)//trim(eventFile),status='new')
  call init_events( rngenerator=ranfun ,writeUnit=eventUnit ,messageUnit=6 &
                   ,desiredPrecision=dsrPrecision ,desiredNevents=dsrNevents )
  write(eventUnit,'(a22,e23.16)') "CENTER-OF-MASS ENERGY:",Ecm
  write(eventUnit,'(a30,e23.16)') "NEGATIVE RAPIDITY BEAM ENERGY:",EnegRap
  write(eventUnit,'(a30,e23.16)') "POSITIVE RAPIDITY BEAM ENERGY:",EposRap
  write(eventUnit,'(A)') "LHAPDF SET: "//"!(lhaSet!)"
  write(eventUnit,'(A)') "TMDLIB SET: "//"!(TMDlibSetA!)" !{withTMDlib0=yes!}
  write(eventUnit,'(A)') "PDFxTMD SET: "//"!(PDFxTMDSetA!)" !{withPDFxTMD0=yes!}
![info4eventFile
!]info4eventFile
  write(eventUnit,'(a22,e23.16)') "ELECTRO-WEAK COUPLING:",gQED**2/(4*r1PI)

  Iev = 0
  do
    Iev = Iev+1
    weight = 0

    call gnrt_proc(iProc,rho)
    associate( prc => proc(iProc) )

    call prc%inst%gnrt( wght_inst ,xBjorken,Qsquare,yInelast ,xA,kTsq ,qMOM,kMOM )
    if (wght_inst.le.rZRO) goto 111

    fluxFac = 8*xA*EposRap*EnegRap

    pInst(0:3,-2)=-qMOM ;prc%sXtrn(-2)=-Qsquare
    pInst(0:3,-1)=-kMOM ;prc%sXtrn(-1)=-kTsq
    call prc%kaleu%gnrt( discard ,psp_q ,pInst ,prc%sXtrn ,mfins=prc%mXtrn(1:Nfinst-1) )
    if (discard) goto 111

    pInst(0:3,-2) = [-EnegRap,rZRO,rZRO,EnegRap]
    call psPoint%put( -2 ,pInst(:,-2) )
    call psPoint%put( -1 ,pInst(:,-1) )
    do jj=1,Nfinst-1
      call psPoint%copy( jj ,psp_q ,jj )
    enddo
    call psPoint%put( Nfinst ,-qMOM-pInst(:,-2) )
    call psPoint%complete

    call breit%init(qMOM)
    do ii=1,Nfinst
      pBreit(0:3,ii) = breit%act_ev(psPoint%p(psPoint%b(ii))%E &
                                   ,psPoint%p(psPoint%b(ii))%V(1:3))
    enddo

![cuts
!]cuts
    INCLUDE 'extra_cuts.h90'

!{itmdf=no
    valPartLumi = pdfunc_lhapdf( prc%parton(-1),xA,scaleA )              !{pdfTypeA=LHAPDF!}
    valPartLumi = pdfunc_grid(-1,prc%parton(-1),xA,scaleA,kTsq )         !{pdfTypeA=grid0!}
    valPartLumi = pdfunc_grid(-1,prc%parton(-1),xA,scaleA,kTsq )         !{pdfTypeA=gridA!}
    valPartLumi = pdfunc_tmdlib( prc%parton(-1),xA,scaleA,kTsq,!(kfA!) ) !{pdfTypeA=TMDlib0!}
    valPartLumi = pdfunc_pdfxtmd( prc%parton(-1),xA,scaleA,kTsq,!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
    valPartLumi = pdfunc_pdfxtmd_cpdf( prc%parton(-1),xA,scaleA,1 ) !{pdfTypeA=PDFxTMDCPDF!}
!}itmdf
!{itmdf=yes
    tmdListA(:) = tmdfunc( xA,scaleA,kTsq )
    valPartLumi = 1
!}itmdf
    if (valPartLumi.eq.rZRO) goto 111

!{helicity=sum
    call prc%me%helicitySummed( matrixelement ,psPoint &
                               ,tmdListA(:) & !{itmdf=yes!}
                               ,colcon=colcon(1:2,1:Nfinst+2) )
!}helicity
!{helicity=sampling
    call rangen(rhl)
    call prc%me%helicitySampled( matrixelement ,psPoint ,rhl &
                                ,tmdListA(:) & !{itmdf=yes!}
                                ,colcon=colcon(1:2,1:Nfinst+2) &
                                ,helicity=helicity(1:Nfinst+2) )
!}helicity

    alphaS = alphasFunc(renScale)

    weight = rONE
    weightFac = wghtCnst/fluxFac &
              * wght_inst * prc%kaleu%wght(psPoint) &
              * valPartLumi * alphaS**prc%pQCD &
              * matrixelement * prc%factor
    if (weightFac.lt.crudeXsec*1d-16) weight=0 ! Out-comment this line if you want 
                                               ! to keep negative-weight events,
                                               ! and recompile with recompile.sh.
    weight = weight*weightFac &
           / (sumProcWght(iProc)-sumProcWght(iProc-1))
    weight = weight * 8*r1PI/gQED**4 *xBjorken*Qsquare**2/(1+(1-yInelast)**2)/389379.66d0 !{DISF2=yes!}
![weights
!]weights
    INCLUDE 'extra_weights.h90'

    111 continue

    call update_events( discard ,accept ,stopMC ,weight )
    if (discard) cycle

    if (accept) then
      associate( psp=>psPoint )
      writeHelicity(-2:-1) = 99
      writeHelicity(1:Nfinst) = 99                                !{helicity!=sampling!}
      writeHelicity(1:Nfinst) = helicity(allInc(1:Nfinst,Nfinst)) !{helicity=sampling!}
      write(eventUnit,'(i3)') iProc
      do jj=-2,psp%Nfinst ;if (jj.eq.0) cycle
        icc(1:2) = colcon(1:2,allInc(jj,Nfinst))
        write(eventUnit,'(5e24.16,2i4,i3)') &
          psp%p(psp%b(jj))%E ,psp%p(psp%b(jj))%V ,psp%p(psp%b(jj))%S &
         ,fivehund(icc(1))+icc(1),fivehund(icc(2))+icc(2),writeHelicity(jj)
      enddo
      write(eventUnit,'(4e24.16)') matrixelement,valPartLumi*xA,alphaS,renScale
!{itmdf=no
      write(eventUnit,'(4e24.16)') rONE,rONE,rZRO,EnegRap
      write(eventUnit,'(4e24.16)') valPartLumi*xA,xA,sqrt(kTsq),scaleA
!}itmdf
!{itmdf=yes
      write(eventUnit,'(4e24.16)') rONE,rONE,rZRO,EnegRap
      write(eventUnit,'(4e24.16)') rONE,xA,sqrt(kTsq),scaleA
!}itmdf
      end associate
      call write_accumulated
    endif

    if (stopMC) exit

    end associate
  enddo

  call write_accumulated
  call write_result

  close(eventUnit)

contains

  subroutine gnrt_proc(ii,rho)
  integer,intent(out) :: ii
  real(kind(1d0)),intent(in) :: rho
  integer :: i0,i1
  i0 = 0
  i1 = NprocTot
  do while (i1-i0.gt.1)
    ii = (i0+i1)/2
    if (rho.le.sumProcWght(ii)) then
      i1 = ii
    else
      i0 = ii
    endif
  enddo
  ii = i1
  end subroutine

  


end program
