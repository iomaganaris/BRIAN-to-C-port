//
// Created by Magkanaris Ioannis on 11.10.20.
//

#pragma once

#define MxM

extern const double defaultclock_dt;

//extern const double taum_S = 10 * 1e-3; //ms
extern const double Ee;
extern const double tuae;
extern const double Fon;
extern const double Foff;

#ifdef NxM
extern const double s;
#endif
#ifdef MxM
extern const double s;
#endif
extern const double Amax;
extern const double Amin;
extern const double Ainit;
extern const double Umax;
extern const double Umin;
extern const double Uinit;

extern const double dFBn;
extern const double dFBp;
extern const double dFFp;

//#Short-term plasticity params
extern const double tau_u;
extern const double tau_r;

//#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
//extern const double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
extern const double AFBn;
extern const double tau_FBn;
extern const double AFBp;
extern const double tau_FBp;
extern const double AFFp;
extern const double tau_FFp;
//#etaU = 0.35
extern const double etaU;
extern const double etaA;
//#etaA = 0.35

//# Adex Parameters
extern const double C;
extern const double gL;
extern const double taum;
extern const double EL;
extern const double DeltaT;
extern const double vti;
//#vtrest = vti + 5 * DeltaT
extern const double vtrest;
extern const double VTmax;
extern const double tauvt;

extern const double tauw;
extern const double c;
extern const double b;
extern const double Vr;
