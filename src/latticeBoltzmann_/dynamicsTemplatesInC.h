#ifndef DYNAMICSTEMPLATESINC_H_INCLUDED
#define DYNAMICSTEMPLATESINC_H_INCLUDED

#include "core_/globalDefsInC.h"
#include "latticeBoltzmann_/basicRoutineInC.h"
#include "latticeBoltzmann_/descriptorInC.h"

double bgk_ma2_collision_InC(struct UniformSequence cellPopulation, double rhoBar, struct UniformSequence_const j, double omega, struct DescriptorInC* descriptor_);

double bgk_ma2_equilibrium_InC (int iPop, double rhoBar, double invRho, struct UniformSequence_const j, double jSqr, struct DescriptorInC* descriptor_);

#endif // DYNAMICSTEMPLATESINC_H_INCLUDED
