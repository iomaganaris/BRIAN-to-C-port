//
// Created by Magkanaris Ioannis on 11.10.20.
//

#include "constants.hpp"

#define MxM

const double defaultclock_dt = 1*1e-3;	//ms

//const double taum_S = 10 * 1e-3; //ms
const double Ee = 0 * 1e-3; //mV
const double tuae = 2 * 1e-3; //ms
const double Fon = 50; //Hz
const double Foff = 3; //Hz

#ifdef NxM
const double s = 100 * 1e-10;//100*1e-10;
#endif
#ifdef MxM
const double s = 500000;	//for testing
#endif
const double Amax = 2.0;
const double Amin = 0;
const double Ainit = 0.1;
const double Umax = 1.0;
const double Umin = 0;
const double Uinit = 0.1;

const double dFBn = 0;
const double dFBp = 0;
const double dFFp = 0;

//#Short-term plasticity params
const double tau_u = 50 * 1e-3;	//ms
const double tau_r = 200 * 1e-3;	//ms

//#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
//const double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
const double AFBn = 0.1771;
const double tau_FBn = 0.0327 * 1e3 * 1e-3;	//ms
const double AFBp = 0.1548;
const double tau_FBp = 0.2302 * 1e3 * 1e-3;	//ms
const double AFFp = 0.0618;
const double tau_FFp = 0.0666 * 1e3 * 1e-3;	//ms
//#etaU = 0.35
const double etaU = 0.15;
const double etaA = 0.15;
//#etaA = 0.35

//# Adex Parameters
const double C = 281*1e-12;	//pF
const double gL = 30*1e-9; //nS
const double taum = 281*1e-12 / 30*1e-9;	// C/gL	// const double initilization of taum(?)
const double EL = -70.6*1e-3;	//mV
const double DeltaT = 2*1e-3;	//mV
const double vti = -50.4*1e-3;	//mV
//#vtrest = vti + 5 * DeltaT
const double vtrest = -45*1e-3;	//mV
const double VTmax = 18*1e-3;	//mV
const double tauvt = 50*1e-3;	//ms

const double tauw = 144*1e-3;	//ms
const double c = 4*1e-9;	//ns
const double b = 0.0805*1e-9;	//nA
const double Vr = -70.6*1e-3;	//mV