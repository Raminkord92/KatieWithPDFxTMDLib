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
  use katie_amplitudes, only: katamp_set_directions
  use katie_matrixelement
  use katie_ranInst
  use katie_partlumi
  use katie_pdfs  !{itmdf=no!}
  use katie_itmds !{itmdf=yes!}
  use katie_model_kaleu_qcd
  use katie_histogramtools, only: breit_type,deltaR,pTrans !{partlumi=DIS!}
  use avh_mctools
  use avh_othercnst
  use avh_lorentz, only: init_lorentz
  use avh_ranvar, only: init_ranvar
  use katie_events
  use katie_version
  implicit none

  integer,parameter :: Ngroup = !(Ngroup!)
  integer,parameter :: NprocTot = !(NprocTot!)
  integer,parameter :: Nfinst(Ngroup) = !(Nfinst!)
  real(kind(1d0)),parameter :: sigma_eff = !(sigma_eff!)

  type(psPoint_type) :: psPoint
  integer :: Iev
  
  real(kind(1d0)) :: Ecm,Esoft,EposRap,EnegRap
  real(kind(1d0)) :: renScale(Ngroup),scaleA(Ngroup),scaleB(Ngroup)

  integer,save :: Nproc(Ngroup)
  real(kind(1d0)),save :: sumProcWght(Ngroup,0:NprocTot)

  abstract interface
    function pdfsum_interface( pdfB ,pdfA ) result(rslt)
    real(kind(1d0)),intent(in) :: pdfB(-6:6),pdfA(-6:6)
    real(kind(1d0)) :: rslt
    end function
  end interface

  type :: procVal_type
    integer :: Nfinst,pQCD,pNonQCD(3),parton(-2:-1)
    integer :: group(Ngroup)
    type(process_type) :: flav
    type(matrixelement_type) :: me
    real(kind(1d0)) :: factor,mXtrn(-2:9),sXtrn(-2:9)
    procedure(pdfsum_interface),pointer,nopass :: pdf
    type(ranInst_type) :: inst
    type(kaleu_base_type) :: kaleu
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
    :: thrs,pp(0:4,1:2**10),wghtCnst(Ngroup),weight,pInstTot(0:3,-2:-1) &
      ,rho(4),kTsq(-2:-1,Ngroup),crudeXsec(Ngroup),weightFac &
      ,fluxFac(Ngroup),sHat(Ngroup),wght_inst(Ngroup),alphaS(Ngroup) &
      ,valPartLumi(Ngroup),matrixelement(Ngroup),ww &
      ,kTA(Ngroup,2),kTB(Ngroup,2),pInst(Ngroup,0:3,-2:-1),xA(Ngroup),xB(Ngroup) &
      ,dsrPrecision,dsrNevents,rcc,rhl,rpl(20),evalnTbl(99),arrayTbl(11,9)
!{itmdf=no
  real(kind(1d0)) :: pdfValA(-6:6),pdfValB(-6:6) !{partlumi=combined!}
  real(kind(1d0)) :: pdfValA(Ngroup),pdfValB(Ngroup) !{partlumi=individual,DIS!}
!}itmdf
!{itmdf=yes
  real(kind(1d0)) :: tmdListA(10,Ngroup)
!}itmdf
  integer &
    :: Iop,ii,jj,iProc(Ngroup),seed,offset,procID(Ngroup,NprocTot),jGroup &
      ,eventUnit,colcon(2,maxval(Nfinst)+2,Ngroup),icc(2),ordrdTbl(11,9) &
      ,helicity(maxval(Nfinst)+2,Ngroup),writeHelicity(-2:maxval(Nfinst))
  integer :: kinID(Ngroup)=0
  logical &
    :: discard,accept,stopMC
  character(maxLineLen) :: prefix,eventFile
  type(procVal_type) :: proc(NprocTot)
  type(kaleu_model_type) :: kaleuModel
  type(commandArgs_type) :: comarg
  type(psPoint_type),save :: psPnt(Ngroup)
!{partlumi=DIS
  integer,parameter :: DISlepton = Nfinst(1)
  real(kind(1d0)) :: Qsquare,yInelast,xBjorken,Qvirtual(0:3),pBreit(0:3,13)
  type(breit_type) :: breit
!}partlumi
!{withMINCAS=yes
  character(4175) :: zeug
  real(kind(1d0)) :: p1MINC(0:3),p2MINC(0:3) !initial momenta
  real(kind(1d0)) :: w1MINC,w2MINC !weights for the jets
  real(kind(1d0)) :: k1MINC(0:3),k2MINC(0:3) !final momenta
  interface
    subroutine mincasinit(inputname1,zeug1)
    character(5) :: inputname1
    character(4175) :: zeug1
    end subroutine
  end interface
  interface
    subroutine mincasleadingparticle(p0,px,py,pz,p10,k1x,k1y,p1z,w)
    real(kind(1d0)) :: p0,px,py,pz
    real(kind(1d0)) :: p10,k1x,k1y,p1z,w
    end subroutine
  end interface
!}withMINCAS


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
!{withTMDlibB=yes
  call pdfs_init_tmdlib('!(TMDlibSetB!)') !{TMDlibKey=char!}
  call pdfs_init_tmdlib(!(TMDlibSetB!))   !{TMDlibKey=int!}
!}withTMDlibB
!{withPDFxTMDCoupling=yes
  call pdfs_init_pdfxtmd_coupling('!(PDFxTMDCouplingSet!)')
!}withPDFxTMDCoupling
!{withPDFxTMD0=yes
  call pdfs_init_pdfxtmd('!(PDFxTMDSetA!)',!(PDFxTMDMemberA!))
!}withPDFxTMD0
!{withPDFxTMDB=yes
  call pdfs_init_pdfxtmd('!(PDFxTMDSetB!)',!(PDFxTMDMemberB!),2)
!}withPDFxTMDB
!{withPDFxTMDCPDFA=yes
  call pdfs_init_pdfxtmd_cpdf('!(PDFxTMDSetA!)',!(PDFxTMDMemberA!))
!}withPDFxTMDCPDFA
!{withPDFxTMDCPDFB=yes
  call pdfs_init_pdfxtmd_cpdf('!(PDFxTMDSetB!)',!(PDFxTMDMemberB!),2)
!}withPDFxTMDCPDFB
![tmds
!]tmds
![pars
!]pars

  do ii=1,NprocTot
    proc(ii)%group = 0
  enddo

  if (ostask(-1).eq.0.and.ostask(-2).ne.0) then
    call katamp_set_directions([-1d0,0d0,0d0,1d0],[-1d0,0d0,0d0,-1d0])
  endif

![kinematics
  call add_kinematics( kinID(1) ,5 ,[1,0,0,0,1] )
!]kinematics

![processes
  proc(1)%group(1) = [1]
  proc(1)%Nfinst = 3
  proc(1)%pNonQCD = [0,0,0]
  call proc(1)%flav%put([gluon,uQuark],[uQuark,gluon,gluon],get_anti)
  proc(1)%factor = 1
  proc(1)%label = 'someLabel' 
  proc(1)%pdf => pdf_g_g
!]processes

  sumProcWght(:,0) = 0
  do ii=1,NprocTot
    associate( prc=>proc(ii) )
    prc%factor = prc%factor &
               / get_symFac(prc%Nfinst,prc%flav%infi(1:prc%Nfinst)) &
               / get_NcolDof(prc%flav%infi(-1))/get_NcolDof(prc%flav%infi(-2))
    call prc%me%init( kinID(prc%group(1)) ,prc%flav%xtrn ,prc%pNonQCD &
                     ,helType='sum' &    !{helicity=sum!}
                     ,helType='no_sum' & !{helicity!=sum!}
!{itmdf=yes
                     ,itmdf=.true. &
                     ,leadingColor=.true. & !{leadingColor=yes!}
!}itmdf
                    )
    prc%pQCD = prc%Nfinst-sum(prc%pNonQCD)+2*prc%pNonQCD(2)
    do jj=-2,prc%Nfinst ;if(jj.eq.0)cycle
      prc%mXtrn(jj) = mass(abs(prc%flav%infi(jj)))
      prc%sXtrn(jj) = prc%mXtrn(jj)**2
    enddo
    open(99,file=trim(prc%label)//'result.dat',status='old')
    read(99,*) ww
    close(99)
    jj = 0 
    do ;jj=jj+1 ;if(jj.gt.Ngroup)exit
      jGroup = prc%group(jj)
      if (jGroup.eq.0) exit
      Nproc(jGroup) = Nproc(jGroup)+1
      procId(jGroup,Nproc(jGroup)) = ii
      sumProcWght(jGroup,Nproc(jGroup)) = sumProcWght(jGroup,Nproc(jGroup)-1) + ww
    enddo
    end associate
  enddo
  do ii=1,Ngroup
    crudeXsec(ii) = sumProcWght(ii,Nproc(ii))
    sumProcWght(ii,1:Nproc(ii)) = sumProcWght(ii,1:Nproc(ii))/crudeXsec(ii)
  enddo

  call rangen_init(seed)

  do jj=1,Ngroup
    wghtCnst(jj) = 389379.66d0 * (r2PI)**(4-3*Nfinst(jj)) / (2*2)
    if (ostask(-2).eq.1) wghtCnst(jj) = wghtCnst(jj)/r1PI
    if (ostask(-1).eq.1) wghtCnst(jj) = wghtCnst(jj)/r1PI
    call psPnt(jj)%alloc( Nfinst(jj) )
  enddo
  call psPoint%alloc( sum(Nfinst(1:Ngroup)) )

  call model_kaleu( kaleuModel )
  call kaleu_set_unit( 'all' ,-1 )
  call kaleu_set_unit( 'progress' ,6 )

  do ii=1,NprocTot
    associate( prc=>proc(ii) )
!{partlumi=individual,combined
    call prc%inst%init( ranfun ,ostask ,prc%Nfinst ,EposRap,EnegRap ,mass(abs(prc%flav%infi(1))) &
        ,option=!(instOption!) ,xmin=prc%Nfinst*Esoft/Ecm ,fileName=prc%label ,readUnit=99 )
!}partlumi
!{partlumi=DIS
    call prc%inst%init( ranfun ,ostask ,EposRap,EnegRap &
                       ,xmin=!(xAmin!) ,xmax=!(xAmax!) ,fileName=prc%label ,readUnit=99 )
!}partlumi
    call prc%kaleu%init( kaleuModel ,prc%flav%apply(flavor_kaleu) ,Ecm,Esoft ,ranfun,ranfun &
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
  write(eventUnit,'(A)') "TMDLIB SET (negrap): "//"!(TMDlibSetB!)" !{withTMDlibB=yes!}
  write(eventUnit,'(A)') "PDFxTMD SET: "//"!(PDFxTMDSetA!)" !{withPDFxTMD0=yes!}
  write(eventUnit,'(A)') "PDFxTMD SET (negrap): "//"!(PDFxTMDSetB!)" !{withPDFxTMDB=yes!}
![info4eventFile
!]info4eventFile
  write(eventUnit,'(a22,e23.16)') "ELECTRO-WEAK COUPLING:",gQED**2/(4*r1PI)
!{withMINCAS=yes
  call mincasinit("input",zeug)
  write(eventUnit,'(A)') zeug 
!}withMINCAS

  Iev = 0
  do
    Iev = Iev+1
    weight = 0

    call rangen(rho(1:Ngroup),Ngroup)
    do ii=1,Ngroup
      call gnrt_proc(ii,iProc(ii),rho(ii))
    enddo

    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
!{partlumi=combined,individual
      call prc%inst%gnrt_AB( & !{pdfTypes=gridA_gridB,grid0_grid0,TMDlib0_TMDlib0,TMDlibA_TMDlibB,PDFxTMD0_PDFxTMD0,PDFxTMDA_PDFxTMDB!}
      call prc%inst%gnrt_A0( & !{pdfTypes=gridA_LHAPDF,grid0_LHAPDF,TMDlib0_LHAPDF,PDFxTMD0_LHAPDF,PDFxTMDA_LHAPDF,PDFxTMD0_PDFxTMDCPDF,PDFxTMDA_PDFxTMDCPDF!}
      call prc%inst%gnrt_0B( & !{pdfTypes=LHAPDF_gridB,LHAPDF_grid0,LHAPDF_TMDlib0,LHAPDF_PDFxTMD0,LHAPDF_PDFxTMDB,PDFxTMDCPDF_gridB,PDFxTMDCPDF_grid0,PDFxTMDCPDF_TMDlib0,PDFxTMDCPDF_PDFxTMD0,PDFxTMDCPDF_PDFxTMDB!}
      call prc%inst%gnrt_00( & !{pdfTypes=LHAPDF_LHAPDF,PDFxTMDCPDF_LHAPDF,LHAPDF_PDFxTMDCPDF,PDFxTMDCPDF_PDFxTMDCPDF!}
                             wght_inst(ii) ,xA(ii),kTA(ii,:) ,xB(ii),kTB(ii,:) &
                           )
      if (wght_inst(ii).le.rZRO) goto 111
!}partlumi
!{partlumi=DIS
      call prc%inst%gnrt_A( & !{pdfTypeA=gridA,grid0,TMDlib0,PDFxTMD0,PDFxTMDA!}
      call prc%inst%gnrt_0( & !{pdfTypeA=LHAPDF!}
      call prc%inst%gnrt_0( & !{pdfTypeA=PDFxTMDCPDF!}
                           wght_inst(ii) ,xA(ii),kTA(ii,:) &
                          )
      if (wght_inst(ii).le.rZRO) goto 111
      xB(ii)=1 ;kTB(ii,:)=0
!}partlumi
!{partlumi=photons
      xA(ii)=1 ;kTA(ii,:)=0 ;xB(ii)=1 ;kTB(ii,:)=0 ;wght_inst(ii)=1
!}partlumi
      end associate
    enddo
    if (sum(xA(1:Ngroup)).gt.rONE) goto 111
    if (sum(xB(1:Ngroup)).gt.rONE) goto 111

    pInstTot = 0
    do ii=1,Ngroup
      pInst(ii,0:3,-1) = [-xA(ii)*EposRap,-kTA(ii,1),-kTA(ii,2),-xA(ii)*EposRap]
      pInst(ii,0:3,-2) = [-xB(ii)*EnegRap,-kTB(ii,1),-kTB(ii,2), xB(ii)*EnegRap]
      pInstTot(:,:) = pInstTot(:,:) + pInst(ii,:,:)
    enddo

    do ii=1,Ngroup
      sHat(ii) = loc_sqr( pInst(ii,:,-1) + pInst(ii,:,-2) )
      if (Nfinst(ii).gt.1.and.sHat(ii).le.0) goto 111
    enddo

    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
      kTsq(-1,ii) = kTA(ii,1)**2 + kTA(ii,2)**2
      kTsq(-2,ii) = kTB(ii,1)**2 + kTB(ii,2)**2
!{fluxFactor=textbook
      fluxFac(ii) = ( sHat(ii) + kTsq(-1,ii) + kTsq(-2,ii) )**2 - 4*kTsq(-1,ii)*kTsq(-2,ii)
      if (fluxFac(ii).le.0) goto 111
      fluxFac(ii) = 2*sqrt(fluxFac(ii))
!}fluxFactor
!{fluxFactor=collinear
      fluxFac(ii) = 8*pInst(ii,0,-2)*pInst(ii,0,-1)
!}fluxFactor
      end associate
    enddo

    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
      prc%sXtrn(-2:-1) =-kTsq(-2:-1,ii)
      call prc%kaleu%gnrt( discard ,psPnt(ii) ,pInst(ii,:,:) &
                          ,prc%sXtrn ,mfins=prc%mXtrn(1:prc%Nfinst) )
      if (discard) goto 111
      end associate
    enddo

    call psPoint%put( -2 ,pInstTot(:,-2) )
    call psPoint%put( -1 ,pInstTot(:,-1) )
    offset = 0
    do ii=1,Ngroup
      do jj=1,Nfinst(ii)
        call psPoint%copy( offset+jj ,psPnt(ii) ,jj )
      enddo
      offset = offset+Nfinst(ii)
    enddo
    call psPoint%complete

!{partlumi=DIS
    Qsquare = 2*EnegRap*psPoint%minus(psPoint%b(DISlepton))
    yInelast = 1 - psPoint%plus(psPoint%b(DISlepton))/(2*EnegRap)
    xBjorken = Qsquare/(4*EnegRap*EposRap*yInelast) 
    Qvirtual = [ EnegRap-psPoint%p(psPoint%b(DISlepton))%E &
               ,        -psPoint%p(psPoint%b(DISlepton))%V(1) &
               ,        -psPoint%p(psPoint%b(DISlepton))%V(2) &
               ,-EnegRap-psPoint%p(psPoint%b(DISlepton))%V(3) ]
    call breit%init(Qvirtual)
    do ii=1,Nfinst(1)
      pBreit(0:3,ii) = breit%act_ev(psPoint%p(psPoint%b(ii))%E &
                                   ,psPoint%p(psPoint%b(ii))%V(1:3))
    enddo
!}partlumi
!{withMINCAS=yes
  p1MINC(0)=psPoint%p(psPoint%b(1))%E;p1MINC(1:3)=psPoint%p(psPoint%b(1))%V(1:3)
  p2MINC(0)=psPoint%p(psPoint%b(2))%E;p2MINC(1:3)=psPoint%p(psPoint%b(2))%V(1:3)
  call mincasleadingparticle( p1MINC(0),p1MINC(1),p1MINC(2),p1MINC(3) &
                             ,k1MINC(0),k1MINC(1),k1MINC(2),k1MINC(3) ,w1MINC )
  call mincasleadingparticle( p2MINC(0),p2MINC(1),p2MINC(2),p2MINC(3) &
                             ,k2MINC(0),k2MINC(1),k2MINC(2),k2MINC(3) ,w2MINC )
  psPoint%p(psPoint%b(1))%E=k1MINC(0);psPoint%p(psPoint%b(1))%V(1:3)=k1MINC(1:3)
  psPoint%p(psPoint%b(2))%E=k2MINC(0);psPoint%p(psPoint%b(2))%V(1:3)=k2MINC(1:3)
  weight = weight*w1MINC*w2MINC
!}withMINCAS
![cuts
!]cuts
    INCLUDE 'extra_cuts.h90'

    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
!{itmdf=no
!{partlumi=individual
      pdfValA(ii)= pdfunc_lhapdf( prc%parton(-1),xA(ii),scaleA(ii) )                     !{pdfTypeA=LHAPDF!}
      pdfValA(ii)= pdfunc_grid(-1,prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=grid0!}
      pdfValA(ii)= pdfunc_grid(-1,prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=gridA!}
      pdfValA(ii)= pdfunc_tmdlib( prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(kfA!) ) !{pdfTypeA=TMDlib0!}
      pdfValA(ii)= pdfunc_tmdlibset(prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(kfA!),!(TMDlibSetA!) ) !{pdfTypeA=TMDlibA!}
      pdfValA(ii)= pdfunc_pdfxtmd( prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
      pdfValA(ii)= pdfunc_pdfxtmdset(prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!) ) !{pdfTypeA=PDFxTMDA!}
      pdfValA(ii)= pdfunc_pdfxtmd_cpdf(prc%parton(-1),xA(ii),scaleA(ii),1) !{pdfTypeA=PDFxTMDCPDF!}
      pdfValB(ii)= pdfunc_lhapdf( prc%parton(-2),xB(ii),scaleB(ii) )                     !{pdfTypeB=LHAPDF!}
      pdfValB(ii)= pdfunc_grid(-1,prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii) )         !{pdfTypeB=grid0!}
      pdfValB(ii)= pdfunc_grid(-2,prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii) )         !{pdfTypeB=gridB!}
      pdfValB(ii)= pdfunc_tmdlib( prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii),!(kfB!) ) !{pdfTypeB=TMDlib0!}
      pdfValB(ii)= pdfunc_tmdlibset(prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii),!(kfB!),!(TMDlibSetB!) ) !{pdfTypeB=TMDlibB!}
      pdfValB(ii)= pdfunc_pdfxtmd( prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii),!(PDFxTMDkfB!),2 ) !{pdfTypeB=PDFxTMD0!}
      pdfValB(ii)= pdfunc_pdfxtmdset(prc%parton(-2),xB(ii),scaleB(ii),kTsq(-2,ii),!(PDFxTMDkfB!),!(PDFxTMDSetB!),!(PDFxTMDMemberB!) ) !{pdfTypeB=PDFxTMDB!}
      pdfValB(ii)= pdfunc_pdfxtmd_cpdf(prc%parton(-2),xB(ii),scaleB(ii),2) !{pdfTypeB=PDFxTMDCPDF!}
      valPartLumi(ii) = pdfValB(ii)*pdfValA(ii)
!}partlumi
!{partlumi=combined
      pdfValA = pdfvec_lhapdf( xA(ii),scaleA(ii) )                     !{pdfTypeA=LHAPDF!}
      pdfValA = pdfvec_grid(-1,xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=grid0!}
      pdfValA = pdfvec_grid(-1,xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=gridA!}
      pdfValA = pdfvec_tmdlib( xA(ii),scaleA(ii),kTsq(-1,ii),!(kfA!) ) !{pdfTypeA=TMDlib0!}
      pdfValA = pdfvec_tmdlibset(xA(ii),scaleA(ii),kTsq(-1,ii),!(kfA!),!(TMDlibSetA!) ) !{pdfTypeA=TMDlibA!}
      pdfValA = pdfvec_pdfxtmd( xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
      pdfValA = pdfvec_pdfxtmdset(xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!) ) !{pdfTypeA=PDFxTMDA!}
      pdfValA = pdfvec_pdfxtmd_cpdf( xA(ii),scaleA(ii),1 ) !{pdfTypeA=PDFxTMDCPDF!}
      pdfValB = pdfvec_lhapdf( xB(ii),scaleB(ii) )                     !{pdfTypeB=LHAPDF!}
      pdfValB = pdfvec_grid(-1,xB(ii),scaleB(ii),kTsq(-2,ii) )         !{pdfTypeB=grid0!}
      pdfValB = pdfvec_grid(-2,xB(ii),scaleB(ii),kTsq(-2,ii) )         !{pdfTypeB=gridB!}
      pdfValB = pdfvec_tmdlib( xB(ii),scaleB(ii),kTsq(-2,ii),!(kfB!) ) !{pdfTypeB=TMDlib0!}
      pdfValB = pdfvec_tmdlibset(xB(ii),scaleB(ii),kTsq(-2,ii),!(kfB!),!(TMDlibSetB!) ) !{pdfTypeB=TMDlibB!}
      pdfValB = pdfvec_pdfxtmd( xB(ii),scaleB(ii),kTsq(-2,ii),!(PDFxTMDkfB!),2 ) !{pdfTypeB=PDFxTMD0!}
      pdfValB = pdfvec_pdfxtmdset(xB(ii),scaleB(ii),kTsq(-2,ii),!(PDFxTMDkfB!),!(PDFxTMDSetB!),!(PDFxTMDMemberB!) ) !{pdfTypeB=PDFxTMDB!}
      pdfValB = pdfvec_pdfxtmd_cpdf( xB(ii),scaleB(ii),2 ) !{pdfTypeB=PDFxTMDCPDF!}
      valPartLumi(ii) = prc%pdf(pdfValB,pdfValA)
!}partlumi
!{partlumi=DIS
      pdfValA(ii) = pdfunc_lhapdf( prc%parton(-1),xA(ii),scaleA(ii) )                     !{pdfTypeA=LHAPDF!}
      pdfValA(ii) = pdfunc_grid(-1,prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=grid0!}
      pdfValA(ii) = pdfunc_grid(-1,prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii) )         !{pdfTypeA=gridA!}
      pdfValA(ii) = pdfunc_tmdlib( prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(kfA!) ) !{pdfTypeA=TMDlib0!}
      pdfValA(ii) = pdfunc_pdfxtmd( prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
      pdfValA(ii) = pdfunc_pdfxtmdset(prc%parton(-1),xA(ii),scaleA(ii),kTsq(-1,ii),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!) ) !{pdfTypeA=PDFxTMDA!}
      pdfValA(ii) = pdfunc_pdfxtmd_cpdf(prc%parton(-1),xA(ii),scaleA(ii),1) !{pdfTypeA=PDFxTMDCPDF!}
      valPartLumi(ii) = pdfValA(ii)
!}partlumi
!}itmdf
!{itmdf=yes
      tmdListA(:,ii) = tmdfunc( xA(ii),scaleA(ii),kTsq(-1,ii) )
      valPartLumi(ii) = pdfunc_lhapdf( prc%parton(-2),xB(ii),scaleB(ii) )
!}itmdf
      if (valPartLumi(ii).eq.rZRO) goto 111
      end associate
    enddo

    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
!{helicity=sum
      call prc%me%helicitySummed( matrixelement(ii) ,psPnt(ii) &
                                 ,tmdListA(:,ii) & !{itmdf=yes!}
                                 ,colcon=colcon(1:2,1:prc%Nfinst+2,ii) )
!}helicity
!{helicity=sampling
      call rangen(rhl)
      call prc%me%helicitySampled( matrixelement(ii) ,psPnt(ii) ,rhl &
                                  ,tmdListA(:,ii) & !{itmdf=yes!}
                                  ,colcon=colcon(1:2,1:prc%Nfinst+2,ii) &
                                  ,helicity=helicity(1:prc%Nfinst+2,ii) )
!}helicity
      end associate
    enddo

    do ii=1,Ngroup
      alphaS(ii) = alphasFunc(renScale(ii))
    enddo

    weight = rONE/sigma_eff
    do ii=1,Ngroup
      associate( prc => proc(procId(ii,iProc(ii))) )
      weightFac = wghtCnst(ii)/fluxFac(ii) &
                * wght_inst(ii) * prc%kaleu%wght(psPnt(ii)) &
                * valPartLumi(ii) * alphaS(ii)**prc%pQCD &
                * matrixelement(ii) * prc%factor
      if (weightFac.lt.crudeXsec(ii)*1d-16) then ! Out-comment these four lines
        weight = 0                               ! if you want to keep negative
        exit                                     ! weight events, and recompile
      endif                                      ! with recompile.sh.    
      weight = weight*weightFac &
             / (sumProcWght(ii,iProc(ii))-sumProcWght(ii,iProc(ii)-1))
      end associate
    enddo
![weights
!]weights
    INCLUDE 'extra_weights.h90'

    111 continue

    call update_events( discard ,accept ,stopMC ,weight )
    if (discard) cycle

    if (accept) then
      do ii=1,Ngroup
        associate( psp=>psPnt(ii) ,prc=>proc(procId(ii,iProc(ii))))
        writeHelicity(-2:-1) = 99
        writeHelicity(1:prc%Nfinst) = 99                                           !{helicity!=sampling!}
        writeHelicity(1:prc%Nfinst) = helicity(allInc(1:prc%Nfinst,prc%Nfinst),ii) !{helicity=sampling!}
        write(eventUnit,'(i3)') procId(ii,iProc(ii))
        do jj=-2,psp%Nfinst ;if (jj.eq.0) cycle
          icc(1:2) = colcon(1:2,allInc(jj,prc%Nfinst),ii)
          write(eventUnit,'(5e24.16,2i4,i3)') &
            psp%p(psp%b(jj))%E ,psp%p(psp%b(jj))%V ,psp%p(psp%b(jj))%S &
           ,fivehund(icc(1))+icc(1),fivehund(icc(2))+icc(2),writeHelicity(jj)
        enddo
        write(eventUnit,'(4e24.16)') matrixelement(ii),valPartLumi(ii)*xA(ii)*xB(ii),alphaS(ii),renScale(ii)
!{itmdf=no
        write(eventUnit,'(4e24.16)') pdfValB(ii)*xB(ii),xB(ii),sqrt(kTsq(-2,ii)),scaleB(ii) !{partlumi=individual!}
        write(eventUnit,'(4e24.16)') pdfValA(ii)*xA(ii),xA(ii),sqrt(kTsq(-1,ii)),scaleA(ii) !{partlumi=individual!}
        write(eventUnit,'(4e24.16)') rONE              ,rONE  ,rZRO             ,EnegRap    !{partlumi=DIS!}
        write(eventUnit,'(4e24.16)') pdfValA(ii)*xA(ii),xA(ii),sqrt(kTsq(-1,ii)),scaleA(ii) !{partlumi=DIS!}
!}itmdf
!{itmdf=yes
        write(eventUnit,'(4e24.16)') valPartLumi(ii)*xB(ii),xB(ii),sqrt(kTsq(-2,ii)),scaleB(ii) !{partlumi=individual!}
        write(eventUnit,'(4e24.16)') rONE                  ,xA(ii),sqrt(kTsq(-1,ii)),scaleA(ii) !{partlumi=individual!}
        write(eventUnit,'(4e24.16)') rONE                  ,rONE  ,rZRO             ,EnegRap    !{partlumi=DIS!}
        write(eventUnit,'(4e24.16)') rONE                  ,xA(ii),sqrt(kTsq(-1,ii)),scaleA(ii) !{partlumi=DIS!}
!}itmdf
        end associate
      enddo
      call write_accumulated
    endif

    if (stopMC) exit

  enddo

  call write_accumulated
  call write_result

  close(eventUnit)

contains

  function loc_sqr(pp) result(rslt)
  intent(in) :: pp
  real(kind(1d0)) :: pp(0:3),rslt
  rslt = (pp(0)-pp(3))*(pp(0)+pp(3))-pp(1)*pp(1)-pp(2)*pp(2)
  end function

  subroutine gnrt_proc(jj,ii,rho)
  integer,intent(in) :: jj
  integer,intent(out) :: ii
  real(kind(1d0)),intent(in) :: rho
  integer :: i0,i1
  i0 = 0
  i1 = Nproc(jj)
  do while (i1-i0.gt.1)
    ii = (i0+i1)/2
    if (rho.le.sumProcWght(jj,ii)) then
      i1 = ii
    else
      i0 = ii
    endif
  enddo
  ii = i1
  end subroutine

end program
