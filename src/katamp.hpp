// The following should match the values of the eponymous parameters
// in the  module avh_lagrangian  of  avhlib.
#define eleNeu 1
#define eleon  2
#define uQuark 3
#define dQuark 4
#define muNeu  5
#define muon   6
#define cQuark 7
#define sQuark 8
#define tauNeu 9
#define tauon  10
#define tQuark 11
#define bQuark 12
#define Wboson 27
#define photon 28
#define Zboson 29
#define gluon  30
#define Higgs  31
// Anti-particles and the W- boson get a minus sign.

// The following definitions work (probably only) with gcc.
#define katamp_set_mass_and_width   __katie_amplitudes_MOD_set_mass_and_width
#define katamp_set_withQCD          __katie_amplitudes_MOD_set_withqcd
#define katamp_set_withQED          __katie_amplitudes_MOD_set_withqed
#define katamp_set_withWeak         __katie_amplitudes_MOD_set_withweak
#define katamp_set_withHiggs        __katie_amplitudes_MOD_set_withhiggs
#define katamp_set_withHG           __katie_amplitudes_MOD_set_withhg
#define katamp_set_withHA           __katie_amplitudes_MOD_set_withha
#define katamp_set_directions  __katie_amplitudes_MOD_katamp_set_directions
#define katamp_add_kinematics  __katie_amplitudes_MOD_katamp_add_kinematics
#define katamp_add_process     __katie_amplitudes_MOD_katamp_add_process
#define katamp_prepare_itmd    __katie_amplitudes_MOD_katamp_prepare_itmd
#define katamp_print_itmd      __katie_amplitudes_MOD_katamp_print_itmd
#define katamp_put_momenta     __katie_amplitudes_MOD_katamp_put_momenta
#define katamp_evaluate        __katie_amplitudes_MOD_katamp_evaluate
#define katamp_sqr             __katie_amplitudes_MOD_katamp_sqr
#define katamp_itmd_sqr        __katie_amplitudes_MOD_katamp_itmd_sqr
// An executable can be created by linking object files, which were 
// created with g++ and gfortran, as follows:
// $ g++ example_cpp_src.o example_fortran_src.o -lgfortran

extern "C"
{

  // particle label as defined above, it must be positive
  void katamp_set_mass_and_width(const unsigned& label,
               const double& mass, const double& width);

  // include QCD interactions (quark-gluon and gluon-gluon)
  void katamp_set_withQCD(const unsigned& flag);

  // include photon interactions (fermion-photon)
  void katamp_set_withQED(const unsigned& flag);

  // include W and Z interactions (couplings to fermions, the photon, and self)
  void katamp_set_withWeak(const unsigned& flag);

  // include Higgs interactions (couplings to fermions, gauge bosons, and self)
  void katamp_set_withHiggs(const unsigned& flag);

  // include an effective Higgs-gluon coupling
  void katamp_set_withHG(const unsigned& flag);

  // include an effective Higgs-photon coupling
  void katamp_set_withHA(const unsigned& flag);

  // set non-default directions for space-like initial states 
  void katamp_set_directions( const double pA[], const double pB[] );

  // The output kinID is an identifier for the kinematics defined by
  // the input (see manual).
  void katamp_add_kinematics( unsigned& kinID,
               const unsigned& Nxtrn, const int inStates[] );

  // The output prcID is an identifier for the process defined by
  // the input (see manual).
  void katamp_add_process( const unsigned& kinID, unsigned& prcIC,
               const int process[],
               const unsigned pNonQCD[] );

  // All input (see manual).
  void katamp_prepare_itmd( const unsigned& kinID, const unsigned& prcIC,
               const unsigned& option );

  // All input (see manual).
  void katamp_print_itmd( const unsigned& kinID, const unsigned& prcIC,
               const unsigned& writeUnit );

  // All input (see manual).
  void katamp_put_momenta( const unsigned& kinID,
               const double momenta[][4],
               const double momSqr[] );

  // The output iEval is an identifier for the amplitude evaluation defined by
  // the input (see manual).
  void katamp_evaluate( const unsigned& kinID, const unsigned& prcIC,
               unsigned& iEval,
               const int helicities[] );

  // rslt output, rest input (see manual).
  void katamp_sqr( const unsigned& kinID, const unsigned& prcIC,
               const unsigned& iEval,
               double& rslt );

  // rslt output, rest input (see manual).
  void katamp_itmd_sqr( const unsigned& kinID, const unsigned& prcIC,
               const unsigned& iEval,
               const double tmds[],
               double& rslt );
}
