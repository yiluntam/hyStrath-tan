#include <map>

// --- Default libraries
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "fvOptions.H"

// --- Smoothing for Local-Time-Stepping
#include "fvcSmooth.H"

// --- High-temperature chemistry library
#include "rho2HTCModel.H"

// --- Relaxation processes
#include "relaxationTimeModel.H"
#include "relaxationTimeModelVV.H"
#include "relaxationTimeModelHE.H"
#include "relaxationTimeModeleV.H"

// --- Multi-species transport models
#include "mixingRule.H"
#include "multiSpeciesTransportModel.H"

// --- Rarefied gas library
#include "rarefactionParameter.H"

// --- Numerical flux functions
#include "numerics/directionInterpolate.H"
#include "numerics/KNP.H"
#include "numerics/KT.H"
#include "numerics/AUSM_Flux.H"
#include "numerics/MUSCL_Kernels.H"
#include "numerics/WENOZ_Kernels.H"
#include "numerics/highOrderKernels.H"

// --- MHD library
#include "mhdModel.H"

// --- Others
#include "wallFvPatch.H"
