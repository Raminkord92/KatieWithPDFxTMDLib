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
  use katie_version
  implicit none

  type(psPoint_type) :: psPoint
  real(kind(1d0)) :: xA,xB
  integer :: Iev
  
  real(kind(1d0)) :: Ecm,Esoft,EposRap,EnegRap
  real(kind(1d0)) :: renScale(1),scaleA(1),scaleB(1)

  integer,parameter :: ostask(-2:-1) = !(offshell!)
  integer,parameter :: Nfinst = !(Nfinst!)
  integer,parameter :: Nflavors = !(Nflavors!)
  integer,parameter :: pNonQCD(3) = !(pNonQCD!)
  real(kind(1d0)) &
    :: thrs,pp(0:4,1:2**10),wghtCnst,weight,kTA(2),kTB(2) &
      ,rho(4),polarRho(20),pInst(0:3,-2:-1),mXtrn(-2:17),sXtrn(-2:17) &
      ,fluxFac,sHat,wght_inst,ampFac,matel &
      ,valPartLumi,matrixelement &
      ,evalnTbl(99),arrayTbl(11,9)
  real(kind(1d0)) :: pdfValA(-6:6),pdfValB(-6:6) !{itmdf=no!}
  real(kind(1d0)) :: tmdListA(10)                !{itmdf=yes!}
  real(kind(1d0)) :: TMDgetQ2max !{withTMDlib0=yes!}
  integer &
    :: Iop,ii,jj,pQCD,seed,Noptim,Nbatch,NbatchGrid,ordrdTbl(11,9) &
      ,kinID=0
  logical &
    :: discard
  type(process_type) :: flavor
  type(kaleu_model_type) :: kaleuModel
  type(kaleu_inst_type) :: kaleu_inst
  type(ranInst_type) :: ranInst
  type(kaleu_stats_type) :: stats
  type(kaleu_base_type) :: kaleu
  type(matrixelement_type) :: meObject
  type(commandArgs_type) :: comarg
  character(256),parameter :: dumpPath=&
![dumpPath
!]dumpPath
!{partlumi=DIS
  integer,parameter :: DISlepton = Nfinst
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

  if (ostask(-1).eq.0.and.ostask(-2).ne.0) then
    call katamp_set_directions([-1d0,0d0,0d0,1d0],[-1d0,0d0,0d0,-1d0])
  endif

![kinematics
  call add_kinematics( kinID ,5 ,[1,0,0,0,1] )
!]kinematics

![processes
  call flavor%put([gluon,uQuark],[uQuark,gluon,gluon],get_anti) ;ampFac=1
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
           / get_NcolDof(flavor%infi(-1))/get_NcolDof(flavor%infi(-2))
  if (ostask(-2).eq.1) wghtCnst = wghtCnst/r1PI
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

  call psPoint%alloc( Nfinst )

  do ii=-2,Nfinst ;if(ii.eq.0)cycle
    mXtrn(ii) = mass(abs(flavor%infi(ii)))
    sXtrn(ii) = mXtrn(ii)**2
  enddo

  call model_kaleu( kaleuModel )
  call kaleu_set_unit( 'all' ,-1 )
  call kaleu_set_unit( 'progress' ,6 )
  call stats%init( ChThrs=!(thrs!) &
                  ,Nbatch=Noptim/max(1,!(Nstep!)) &
                  ,Nstep=!(Nstep!) &
                  ,NbatchGrid=Noptim/max(1,!(NstepGrid!)) &
                  ,NstepGrid=!(NstepGrid!) )

!{partlumi=individual,combined
  call ranInst%init( ranfun ,ostask ,Nfinst ,EposRap,EnegRap ,mass(abs(flavor%infi(1))) &
                    ,option=!(instOption!) ,xmin=Nfinst*Esoft/Ecm )
!}partlumi
!{partlumi=DIS
  call ranInst%init( ranfun ,ostask ,EposRap,EnegRap ,xmin=!(xAmin!) ,xmax=!(xAmax!) )
!}partlumi
  call kaleu%init( kaleuModel ,flavor%apply(flavor_kaleu) ,Ecm,Esoft &
                  ,ranfun,ranfun ,stats ,cancel=discard )
  if (discard) stop
  call mincasinit("../input",zeug) !{withMINCAS=yes!}
  
  Iev = 0
  do while (stats%optimPhase)
    Iev = Iev+1

    call generate_event

    call stats%collect( weight )
    call ranInst%adapt( stats%absWeight )
    call kaleu%adapt( stats )
    call stats%update()

  enddo

  call ranInst%dump( trim(dumpPath) ,99 )
  !call ranInst%plot( 99 )
  call kaleu%dump( trim(dumpPath) ,99 )
  open(99,file=trim(dumpPath)//'result.dat' ,status='replace')
  write(99,*) stats%optAvePos
  close(99)


contains


  subroutine generate_event
    weight = 0
!{partlumi=combined,individual
    call ranInst%gnrt_AB( & !{pdfTypes=gridA_gridB,grid0_grid0,TMDlib0_TMDlib0,TMDlibA_TMDlibB,PDFxTMD0_PDFxTMD0,PDFxTMDA_PDFxTMDB!}
    call ranInst%gnrt_A0( & !{pdfTypes=gridA_LHAPDF,grid0_LHAPDF,TMDlib0_LHAPDF,PDFxTMD0_LHAPDF,PDFxTMDA_LHAPDF,PDFxTMD0_PDFxTMDCPDF,PDFxTMDA_PDFxTMDCPDF!}
    call ranInst%gnrt_0B( & !{pdfTypes=LHAPDF_gridB,LHAPDF_grid0,LHAPDF_TMDlib0,LHAPDF_PDFxTMD0,LHAPDF_PDFxTMDB,PDFxTMDCPDF_gridB,PDFxTMDCPDF_grid0,PDFxTMDCPDF_TMDlib0,PDFxTMDCPDF_PDFxTMD0,PDFxTMDCPDF_PDFxTMDB!}
    call ranInst%gnrt_00( & !{pdfTypes=LHAPDF_LHAPDF,PDFxTMDCPDF_LHAPDF,LHAPDF_PDFxTMDCPDF,PDFxTMDCPDF_PDFxTMDCPDF!}
                          wght_inst ,xA,kTA(:) ,xB,kTB(:) &
                        )
    if (wght_inst.le.rZRO) goto 111
!}partlumi
!{partlumi=DIS
    call ranInst%gnrt_A( & !{pdfTypeA=gridA,grid0,TMDlib0,PDFxTMD0,PDFxTMDA!}
    call ranInst%gnrt_0( & !{pdfTypeA=LHAPDF!}
    call ranInst%gnrt_0( & !{pdfTypeA=PDFxTMDCPDF!}
                        wght_inst ,xA,kTA(:) &
                       )
    if (wght_inst.le.rZRO) goto 111
    xB=1 ;kTB=0
!}partlumi
!{partlumi=photons
    xA=1 ;kTA=0 ;xB=1 ;kTB=0 ;wght_inst=1
!}partlumi

    pInst(0:3,-1) = [-xA*EposRap,-kTA(1),-kTA(2),-xA*EposRap]
    pInst(0:3,-2) = [-xB*EnegRap,-kTB(1),-kTB(2), xB*EnegRap]

    sHat = loc_sqr( pInst(:,-1) + pInst(:,-2) )
    if (Nfinst.gt.1.and.sHat.le.0) goto 111

    sXtrn(-1) =-( kTA(1)**2 + kTA(2)**2 )
    sXtrn(-2) =-( kTB(1)**2 + kTB(2)**2 )
!{fluxFactor=textbook
    fluxFac = ( sHat - sXtrn(-1) - sXtrn(-2) )**2 - 4*sXtrn(-1)*sXtrn(-2)
    if (fluxFac.le.0) goto 111
    fluxFac = 2*sqrt(fluxFac)
!}fluxFactor
!{fluxFactor=collinear
    fluxFac = 8*pInst(0,-2)*pInst(0,-1)
!}fluxFactor
    call kaleu%gnrt( discard ,psPoint ,pInst ,sXtrn ,mfins=mXtrn(1:Nfinst) )
    if (discard) goto 111

!{partlumi=DIS
    Qsquare = 2*EnegRap*psPoint%minus(psPoint%b(DISlepton))
    yInelast = 1 - psPoint%plus(psPoint%b(DISlepton))/(2*EnegRap)
    xBjorken = Qsquare/(4*EnegRap*EposRap*yInelast) 
    Qvirtual = [ EnegRap-psPoint%p(psPoint%b(DISlepton))%E &
               ,        -psPoint%p(psPoint%b(DISlepton))%V(1) &
               ,        -psPoint%p(psPoint%b(DISlepton))%V(2) &
               ,-EnegRap-psPoint%p(psPoint%b(DISlepton))%V(3) ]
    call breit%init(Qvirtual)
    do ii=1,Nfinst
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
    INCLUDE '../extra_cuts.h90'

    valPartLumi = 1
!{itmdf=no
!{partlumi=individual
    valPartLumi = valPartLumi &
      * pdfunc_lhapdf( !(partonA!),xA,scaleA(1) ) &                    !{pdfTypeA=LHAPDF!}
      * pdfunc_grid(-1,!(partonA!),xA,scaleA(1),-sXtrn(-1) ) &         !{pdfTypeA=grid0!}
      * pdfunc_grid(-1,!(partonA!),xA,scaleA(1),-sXtrn(-1) ) &         !{pdfTypeA=gridA!}
      * pdfunc_tmdlib( !(partonA!),xA,scaleA(1),-sXtrn(-1),!(kfA!) ) & !{pdfTypeA=TMDlib0!}
      * pdfunc_tmdlibset(!(partonA!),xA,scaleA(1),-sXtrn(-1),!(kfA!),!(TMDlibSetA!) ) & !{pdfTypeA=TMDlibA!}
      * pdfunc_pdfxtmd( !(partonA!),xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),1 ) & !{pdfTypeA=PDFxTMD0!}
      * pdfunc_pdfxtmdset(!(partonA!),xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!)) & !{pdfTypeA=PDFxTMDA!}
      * pdfunc_pdfxtmd_cpdf(!(partonA!),xA,scaleA(1),1) & !{pdfTypeA=PDFxTMDCPDF!}
      * pdfunc_lhapdf( !(partonB!),xB,scaleB(1) )                      !{pdfTypeB=LHAPDF!}
      * pdfunc_grid(-1,!(partonB!),xB,scaleB(1),-sXtrn(-2) )           !{pdfTypeB=grid0!}
      * pdfunc_grid(-2,!(partonB!),xB,scaleB(1),-sXtrn(-2) )           !{pdfTypeB=gridB!}
      * pdfunc_tmdlib( !(partonB!),xB,scaleB(1),-sXtrn(-2),!(kfB!) )   !{pdfTypeB=TMDlib0!}
      * pdfunc_tmdlibset(!(partonB!),xB,scaleB(1),-sXtrn(-2),!(kfB!),!(TMDlibSetB!) )   !{pdfTypeB=TMDlibB!}
      * pdfunc_pdfxtmd( !(partonB!),xB,scaleB(1),-sXtrn(-2),!(PDFxTMDkfB!),2 )   !{pdfTypeB=PDFxTMD0!}
      * pdfunc_pdfxtmdset(!(partonB!),xB,scaleB(1),-sXtrn(-2),!(PDFxTMDkfB!),!(PDFxTMDSetB!),!(PDFxTMDMemberB!))   !{pdfTypeB=PDFxTMDB!}
      * pdfunc_pdfxtmd_cpdf(!(partonB!),xB,scaleB(1),2)   !{pdfTypeB=PDFxTMDCPDF!}
!}partlumi
!{partlumi=combined
    pdfValA = pdfvec_lhapdf( xA,scaleA(1) )                    !{pdfTypeA=LHAPDF!}
    pdfValA = pdfvec_grid(-1,xA,scaleA(1),-sXtrn(-1) )         !{pdfTypeA=grid0!}
    pdfValA = pdfvec_grid(-1,xA,scaleA(1),-sXtrn(-1) )         !{pdfTypeA=gridA!}
    pdfValA = pdfvec_tmdlib( xA,scaleA(1),-sXtrn(-1),!(kfA!) ) !{pdfTypeA=TMDlib0!}
    pdfValA = pdfvec_tmdlibset(xA,scaleA(1),-sXtrn(-1),!(kfA!),!(TMDlibSetA!) ) !{pdfTypeA=TMDlibA!}
    pdfValA = pdfvec_pdfxtmd( xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
    pdfValA = pdfvec_pdfxtmdset(xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!) ) !{pdfTypeA=PDFxTMDA!}
    pdfValA = pdfvec_pdfxtmd_cpdf( xA,scaleA(1),1 ) !{pdfTypeA=PDFxTMDCPDF!}
    pdfValB = pdfvec_lhapdf( xB,scaleB(1) )                    !{pdfTypeB=LHAPDF!}
    pdfValB = pdfvec_grid(-1,xB,scaleB(1),-sXtrn(-2) )         !{pdfTypeB=grid0!}
    pdfValB = pdfvec_grid(-2,xB,scaleB(1),-sXtrn(-2) )         !{pdfTypeB=gridB!}
    pdfValB = pdfvec_tmdlib( xB,scaleB(1),-sXtrn(-2),!(kfB!) ) !{pdfTypeB=TMDlib0!}
    pdfValB = pdfvec_tmdlibset(xB,scaleB(1),-sXtrn(-2),!(kfB!),!(TMDlibSetB!) ) !{pdfTypeB=TMDlibB!}
    pdfValB = pdfvec_pdfxtmd( xB,scaleB(1),-sXtrn(-2),!(PDFxTMDkfB!),2 ) !{pdfTypeB=PDFxTMD0!}
    pdfValB = pdfvec_pdfxtmdset(xB,scaleB(1),-sXtrn(-2),!(PDFxTMDkfB!),!(PDFxTMDSetB!),!(PDFxTMDMemberB!) ) !{pdfTypeB=PDFxTMDB!}
    pdfValB = pdfvec_pdfxtmd_cpdf( xB,scaleB(1),2 ) !{pdfTypeB=PDFxTMDCPDF!}
    valPartLumi = pdf_!(instPartons!)(pdfValB,pdfValA)
!}partlumi
!{partlumi=DIS
    valPartLumi = pdfunc_lhapdf( !(partonA!),xA,scaleA(1) )                    !{pdfTypeA=LHAPDF!}
    valPartLumi = pdfunc_grid(-1,!(partonA!),xA,scaleA(1),-sXtrn(-1) )         !{pdfTypeA=grid0!}
    valPartLumi = pdfunc_grid(-1,!(partonA!),xA,scaleA(1),-sXtrn(-1) )         !{pdfTypeA=gridA!}
    valPartLumi = pdfunc_tmdlib( !(partonA!),xA,scaleA(1),-sXtrn(-1),!(kfA!) ) !{pdfTypeA=TMDlib0!}
    valPartLumi = pdfunc_pdfxtmd( !(partonA!),xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),1 ) !{pdfTypeA=PDFxTMD0!}
    valPartLumi = pdfunc_pdfxtmdset(!(partonA!),xA,scaleA(1),-sXtrn(-1),!(PDFxTMDkfA!),!(PDFxTMDSetA!),!(PDFxTMDMemberA!)) !{pdfTypeA=PDFxTMDA!}
    valPartLumi = pdfunc_pdfxtmd_cpdf(!(partonA!),xA,scaleA(1),1) !{pdfTypeA=PDFxTMDCPDF!}
!}partlumi
!}itmdf
!{itmdf=yes
      tmdListA(:) = tmdfunc( xA ,scaleA(1) ,-sXtrn(-1) )
      valPartLumi = pdfunc_lhapdf( !(partonB!) ,xB ,scaleB(1) )
!}itmdf

    if (valPartLumi.eq.rZRO) goto 111
!{helicity=sum
    call meObject%helicitySummed( matrixelement ,psPoint &
                                 ,tmdListA & !{itmdf=yes!}
                                )
!}helicity
!{helicity=sampling
    call rangen(rho(1))
    call meObject%helicitySampled( matrixelement ,psPoint ,rho(1) &
                                  ,tmdListA & !{itmdf=yes!}
                                 )
!}helicity

    weight = wghtCnst/fluxFac * kaleu%wght(psPoint) * wght_inst &
           * valPartLumi * matrixelement * ampFac * alphasFunc(renScale(1))**pQCD

![weights
!]weights
    INCLUDE '../extra_weights.h90'

    111 continue
  end subroutine


  function loc_sqr(pp) result(rslt)
  intent(in) :: pp
  real(kind(1d0)) :: pp(0:3),rslt
  rslt = (pp(0)-pp(3))*(pp(0)+pp(3))-pp(1)*pp(1)-pp(2)*pp(2)
  end function


end program
