module katie_model_kaleu_qcd
  use avh_kaleu_model
  use avh_lagrangian
  implicit none

contains

  function flavor_kaleu( ii ) result(rslt)
  integer,intent(in) :: ii
  integer :: rslt
  select case (ii)
  case default ;rslt = abs(ii)
  end select
  end function

  subroutine model_kaleu( model )
!********************************************************************
!********************************************************************
  type(kaleu_model_type),intent(out) :: model
!
  call addparticle( eleNeu )
  call addparticle( eleon  )
  call addparticle( muNeu  )
  call addparticle( muon   )
  call addparticle( tauNeu )
  call addparticle( tauon  )
  call addparticle( uQuark )
  call addparticle( dQuark )
  call addparticle( cQuark )
  call addparticle( sQuark )
  call addparticle( tQuark )
  call addparticle( bQuark )
  call addparticle( Wboson )
  call addparticle( photon )
  call addparticle( Zboson )
  call addparticle( gluon  )
  call addparticle( Higgs  )
!
  if (withQCD) then
    call model%addVertex(  gluon,gluon,gluon )
    call model%addVertex( uQuark,uQuark,gluon )
    call model%addVertex( dQuark,dQuark,gluon )
    call model%addVertex( cQuark,cQuark,gluon )
    call model%addVertex( sQuark,sQuark,gluon )
    call model%addVertex( tQuark,tQuark,gluon )
    call model%addVertex( bQuark,bQuark,gluon )
  endif
!
  if (withQED) then
    call model%addVertex(  gluon, gluon,photon )
    call model%addVertex( uQuark,uQuark,photon )
    call model%addVertex( dQuark,dQuark,photon )
    call model%addVertex( cQuark,cQuark,photon )
    call model%addVertex( sQuark,sQuark,photon )
    call model%addVertex( tQuark,tQuark,photon )
    call model%addVertex( bQuark,bQuark,photon )
    call model%addVertex( Wboson,Wboson,photon )
    call model%addVertex(  eleon, eleon,photon )
    call model%addVertex(   muon,  muon,photon )
    call model%addVertex(  tauon, tauon,photon )
  endif
!
  if (withWeak) then
    call model%addVertex(  gluon, gluon,Zboson )
    call model%addVertex( uQuark,uQuark,Zboson )
    call model%addVertex( dQuark,dQuark,Zboson )
    call model%addVertex( cQuark,cQuark,Zboson )
    call model%addVertex( sQuark,sQuark,Zboson )
    call model%addVertex( tQuark,tQuark,Zboson )
    call model%addVertex( bQuark,bQuark,Zboson )
!  
    call model%addVertex( uQuark,dQuark,Wboson )
    call model%addVertex( cQuark,sQuark,Wboson )
    call model%addVertex( tQuark,bQuark,Wboson )
!  
    call model%addVertex( Wboson,Wboson,Zboson )
!  
    call model%addVertex(  eleon, eleon,Zboson )
    call model%addVertex( eleNeu,eleNeu,Zboson )
    call model%addVertex(  eleon,eleNeu,Wboson )
!  
    call model%addVertex(   muon,  muon,Zboson )
    call model%addVertex(  muNeu, muNeu,Zboson )
    call model%addVertex(   muon, muNeu,Wboson )
!  
    call model%addVertex(  tauon, tauon,Zboson )
    call model%addVertex( tauNeu,tauNeu,Zboson )
    call model%addVertex(  tauon,tauNeu,Wboson )
  endif
!  
  if (withHiggs) then
    call model%addVertex(  Higgs, Higgs,Higgs )
    call model%addVertex( Zboson,Zboson,Higgs )
    call model%addVertex( Wboson,Wboson,Higgs )
    call model%addVertex( tQuark,tQuark,Higgs )
    call model%addVertex( bQuark,bQuark,Higgs )
    call model%addVertex(  tauon, tauon,Higgs )
!    call model%addVertex( cQuark,cQuark,Higgs )
!    call model%addVertex( sQuark,sQuark,Higgs )
!    call model%addVertex(   muon,  muon,Higgs )
!    call model%addVertex( uQuark,uQuark,Higgs )
!    call model%addVertex( dQuark,dQuark,Higgs )
!    call model%addVertex(  eleon, eleon,Higgs )
  endif
!
  if (withHG) then
    call model%addVertex( Higgs,gluon,gluon )
  endif
!
  if (withHA) then
    call model%addVertex( Higgs,photon,photon )
  endif
!
  contains
!
    subroutine addparticle( ii )
    integer,intent(in) :: ii
    call model%addParticle( ii ,prnt_particle(ii) ,mass(ii),width(ii) )
    end subroutine
!
  end subroutine

end module
