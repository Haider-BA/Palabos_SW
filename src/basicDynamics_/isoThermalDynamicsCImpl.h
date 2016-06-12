#ifndef ISOTHERMALDYNAMICSCIMPL_H_INCLUDED
#define ISOTHERMALDYNAMICSCIMPL_H_INCLUDED


#include "core_/globalDefsInC.h"
#include "core_/blockStatisticsInC.h"
#include "core_/dynamicsInC.h"

void collide_BGKDynamics (struct DynamicsInC* _this, struct UniformSequence cellPopulation, struct BlockStatisticsInC* statistics);

void collideExternal_BGKDynamics(struct DynamicsInC* _this, struct UniformSequence cellPopulation, double rhoBar, struct UniformSequence_const j, struct BlockStatisticsInC* statistics, double thetaBar);

double computeEquilibrium_BGKDynamics (struct DynamicsInC* _this, int iPop, double rhoBar, struct UniformSequence_const j, double jSqr, double thetaBar);


#endif // ISOTHERMALDYNAMICSCIMPL_H_INCLUDED
