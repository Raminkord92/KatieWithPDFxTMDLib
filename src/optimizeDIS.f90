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
  use avh_kaleu_stats
  use avh_kaleu_base
  use katie_matrixelement
  use katie_DISinst
  use katie_pdfs  !{itmdf=no!}
  use katie_itmds !{itmdf=yes!}
  use katie_model_kaleu_qcd
  use katie_histogramtools, only: breit_type,deltaR,pTrans
  use avh_mctools
  use avh_othercnst
  use avh_lorentz, only: init_lorentz
  use avh_ranvar, only: init_ranvar
  use katie_version
  implicit none

  type(psPoint_type) :: psPoint,psp_q
  integer :: Iev
  
  real(kind(1d0)) :: Ecm,Esoft,EposRap,EnegRap
  real(kind(1d0)) :: renScale,scaleA,scaleB

  integer,parameter :: ostask(-2:-1) = !(offshell!)
  integer,parameter :: Nfinst = !(Nfinst!)
  integer,parameter :: Nflavors = !(Nflavors!)
  integer,parameter :: pNonQCD(3) = !(pNonQCD!)
  real(kind(1d0)) &
    :: wghtCnst,weight,kTsq,rho,pInst(0:3,-2:-1),mXtrn(-2:17),sXtrn(-2:-1) &
      ,fluxFac,wght_inst,ampFac,valPartLumi,matrixelement,qMOM(0:3),kMOM(0:3) &
      ,evalnTbl(99),arrayTbl(11,9),xA
  real(kind(1d0)) :: tmdListA(10) !{itmdf=yes!}
  real(kind(1d0)) :: TMDgetQ2max !{withTMDlib0=yes!}
  integer &
    :: Iop,ii,jj,pQCD,seed,Noptim,Nbatch,NbatchGrid,ordrdTbl(11,9) &
      ,kinID=0
  logical &
    :: discard
  type(process_type) :: flavor,psflavor
  type(kaleu_model_type) :: kaleuModel
  type(kaleu_inst_type) :: kaleu_inst
  type(DISinst_type) :: DISinst
  type(kaleu_stats_type) :: stats
  type(kaleu_base_type) :: kaleu
  type(matrixelement_type) :: meObject
  type(commandArgs_type) :: comarg
  character(256),parameter :: dumpPath=&
![dumpPath
!]dumpPath
  integer,parameter :: DISlepton = Nfinst
  real(kind(1d0)) :: Qsquare,yInelast,xBjorken,pBreit(0:3,13)
  type(breit_type) :: breit

  call set_unit('banner',-1)
  call init_mathcnst
  call init_ranvar
  call init_lorentz

  call version

  call comarg%grab
  seed=comarg%find_int('seed=') ;if(seed.eq.0)seed=12345    
  Noptim=comarg%find_int('Noptim=') ;if(Noptim.eq.0)Noptim=100000

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
  call add_kinematics( kinID ,5 ,[1,0,0,0,1] )
!]kinematics

![processes
  flavor%put([eleon,uQuark],[uQuark,gluon,eleon],get_anti)
  psflavor%put([photon,uQuark],[uQuark,gluon],get_anti)
!]processes
  if (Nfinst.ne.flavor%Nfinst) then
    if (errru.ge.0) write(errru,*) 'ERROR in optimize: Nfinst incompatible'
    stop
  endif
  ampFac = ampFac/get_symFac(Nfinst,flavor%infi(1:Nfinst))
  pQCD = Nfinst-sum(pNonQCD)+2*pNonQCD(2)

  call rangen_init(seed)

  wghtCnst = 389379.66d0 &
           * (r2PI)**(4-3*Nfinst) &
           / (2*2) &
           / get_NcolDof(flavor%infi(-1))
  if (ostask(-1).eq.1) wghtCnst = wghtCnst/r1PI

  call meObject%init( kinID ,flavor%xtrn ,pNonQCD &
                     ,helType='sum' &    !{helicity=sum!}
                     ,helType='no_sum' & !{helicity!=sum!}
!{itmdf=yes
                     ,itmdf=.true. &
                     ,leadingColor=.true. & !{leadingColor=yes!}
!}itmdf
                    )
  call meObject%print_itmd( 6 ) !{itmdf=yes!}

  call psp_q%alloc(Nfinst-1)
  call psPoint%alloc( Nfinst )

  do ii=-2,Nfinst ;if(ii.eq.0)cycle
    mXtrn(ii) = mass(abs(flavor%infi(ii)))
    if (ii.lt.0) sXtrn(ii) = mXtrn(ii)**2
  enddo

  call model_kaleu( kaleuModel )
  call kaleu_set_unit( 'all' ,-1 )
  call kaleu_set_unit( 'progress' ,6 )
  call stats%init( ChThrs=!(thrs!) &
                  ,Nbatch=Noptim/max(1,!(Nstep!)) &
                  ,Nstep=!(Nstep!) &
                  ,NbatchGrid=Noptim/max(1,!(NstepGrid!)) &
                  ,NstepGrid=!(NstepGrid!) )
  call DISinst%init( ranfun ,ostask(-1) ,EposRap,EnegRap ,Nfinst=Nfinst-1 &
                    ,xBmin=!(xBmin!) ,xBmax=!(xBmax!) &
                    ,QsqMin=!(QsqMin!) ,QsqMax=!(QsqMax!) )
  call kaleu%init( kaleuModel ,psflavor%apply(flavor_kaleu) ,Ecm,Esoft &
                  ,ranfun,ranfun ,stats ,cancel=discard )
  if (discard) stop
  
  Iev = 0
  do while (stats%optimPhase)
    Iev = Iev+1

    call generate_event

    call stats%collect( weight )
    call DISinst%adapt( stats%absWeight )
    call kaleu%adapt( stats )
    call stats%update()

  enddo

  call DISinst%dump( trim(dumpPath) ,99 )
  !call ranInst%plot( 99 )
  call kaleu%dump( trim(dumpPath) ,99 )
  open(99,file=trim(dumpPath)//'result.dat' ,status='replace')
  write(99,*) stats%optAvePos
  close(99)

contains

  subroutine generate_event
    weight = 0
    call DISinst%gnrt( wght_inst ,xBjorken,Qsquare,yInelast ,xA,kTsq ,qMOM,kMOM )
    if (wght_inst.le.rZRO) goto 111

    fluxFac = 8*xA*EposRap*EnegRap
    pInst(0:3,-2)=-qMOM ;sXtrn(-2)=-Qsquare
    pInst(0:3,-1)=-kMOM ;sXtrn(-1)=-kTsq
    call kaleu%gnrt( discard ,psp_q ,pInst ,sXtrn ,mfins=mXtrn(1:Nfinst-1) )
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
    INCLUDE '../extra_cuts.h90'

!{itmdf=no
    valPartLumi = pdfunc_lhapdf( !(partonA!),xA,scaleA )              !{pdfTypeA=LHAPDF!}
    valPartLumi = pdfunc_grid(-1,!(partonA!),xA,scaleA,kTsq )         !{pdfTypeA=grid0!}
    valPartLumi = pdfunc_grid(-1,!(partonA!),xA,scaleA,kTsq )         !{pdfTypeA=gridA!}
    valPartLumi = pdfunc_tmdlib( !(partonA!),xA,scaleA,kTsq,!(kfA!) ) !{pdfTypeA=TMDlib0!}
    valPartLumi = pdfunc_pdfxtmd( !(partonA!),xA,scaleA,kTsq,!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
    valPartLumi = pdfunc_pdfxtmd_cpdf( !(partonA!),xA,scaleA,1 ) !{pdfTypeA=PDFxTMDCPDF!}
!}itmdf
!{itmdf=yes
    tmdListA(:) = tmdfunc( xA ,scaleA) ,kTsq )
    valPartLumi = 1
!}itmdf

    if (valPartLumi.eq.rZRO) goto 111
!{helicity=sum
    call meObject%helicitySummed( matrixelement ,psPoint &
                                 ,tmdListA & !{itmdf=yes!}
                                )
!}helicity
!{helicity=sampling
    call rangen(rho)
    call meObject%helicitySampled( matrixelement ,psPoint ,rho &
                                  ,tmdListA & !{itmdf=yes!}
                                 )
!}helicity

    weight = wghtCnst/fluxFac * kaleu%wght(psPoint) * wght_inst &
           * valPartLumi * matrixelement * ampFac * alphasFunc(renScale)**pQCD
    weight = weight * 8*r1PI/gQED**4 *xBjorken*Qsquare**2/(1+(1-yInelast)**2)/389379.66d0 !{DISF2=yes!}
![weights
!]weights
    INCLUDE '../extra_weights.h90'

    111 continue
  end subroutine


end program
