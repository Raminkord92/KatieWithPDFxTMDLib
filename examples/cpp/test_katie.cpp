//
// Let $KaTie be the KaTie-directory. Executing
//   $ $KaTie/run.sh katamp
// creates the sourefile katamp.f90 in the directory $KaTie/build .
// The header file katamp.hpp can be found in $KaTie/src .
// To create an executable a.out, execute
//   $ gfortran -c katamp.f90
//   $ g++ -c test_katie.cpp
//   $ g++ katamp.o test_katie.o -lgfortran
//
#include <iostream>
#include "katamp.hpp"
int main() {

unsigned int kinID, prcID, iEval[4];
double rslt[4];

unsigned int Nxtrn = 6;
int inStates[6] =  {2,0,0,0,0,2};
int process[6] = {gluon,uQuark,gluon,dQuark,-dQuark,-uQuark};
unsigned int pNonQCD[3] = {0,0,0};

int helicities1[6] = { 1,-1, 1,-1, 1, 1};
int helicities2[6] = {-1,-1, 1,-1, 1, 1};
int helicities3[6] = { 1,-1, 1,-1, 1,-1};
int helicities4[6] = {-1,-1, 1,-1, 1,-1};

double momenta[6][4] ={
{-0.3070889469862555e+3, 0.8702189678979154e+1, 0.1074153121985460e+1,-0.3070889469862555e+3},
{ 0.2630067717703687e+3, 0.7797115039865352e+2, 0.1891380599398776e+2,-0.2504702170825968e+3},
{ 0.3858906671044585e+2,-0.1124753735780160e+1,-0.3207140035914875e+2, 0.2143073208316230e+2},
{ 0.3966675878166514e+2,-0.2751169947612465e+2,-0.2441987877466492e+2,-0.1484006956886930e+2},
{ 0.2626285233665797e+3,-0.5138314615973688e+2, 0.4162921248704136e+2, 0.2541663279117554e+3},
{-0.2968021736428038e+3,-0.6653740705990984e+1,-0.5125892469200897e+1, 0.2968021736428038e+3}
};

katamp_add_kinematics( kinID ,Nxtrn ,inStates );
katamp_add_process( kinID ,prcID ,process ,pNonQCD );
katamp_put_momenta( kinID ,momenta);
katamp_evaluate( kinID ,prcID ,iEval[0] ,helicities1 );
katamp_evaluate( kinID ,prcID ,iEval[1] ,helicities2 );
katamp_evaluate( kinID ,prcID ,iEval[2] ,helicities3 );
katamp_evaluate( kinID ,prcID ,iEval[3] ,helicities4 );
katamp_sqr( kinID ,prcID ,iEval[0] ,rslt[0] );
katamp_sqr( kinID ,prcID ,iEval[1] ,rslt[1] );
katamp_sqr( kinID ,prcID ,iEval[2] ,rslt[2] );
katamp_sqr( kinID ,prcID ,iEval[3] ,rslt[3] );

double sum = 0.0;
for (int entry=0 ;entry<4 ;entry++) sum += rslt[entry];

std::cout << "Ratio result/pre-calculated: "
          << sum/0.4718171604704535e+3
          << std::endl;

return 0;
}
