


#ifndef BASICROUTINEINC_H_INCLUDED
#define BASICROUTINEINC_H_INCLUDED

#include "core_/globalDefsInC.h"
#include "latticeBoltzmann_/descriptorInC.h"

/*struct UniformSequence;
struct UniformSequence_const;
struct DescriptorInC;*/


//------------------------------- Geometric operations -------------------------
double normSqr_InC(struct UniformSequence_const);

//------------------------------- Physical fields ----------------------------

void get_rhoBar_j_InC(struct UniformSequence_const, double*, struct UniformSequence, struct DescriptorInC*);

#endif // BASICROUTINEINC_H_INCLUDED
