
#include "basicDynamics_/isoThermalDynamicsCImpl.h"
#include "latticeBoltzmann_/basicRoutineInC.h"
#include "latticeBoltzmann_/dynamicsTemplatesInC.h"

void collide_BGKDynamics (struct DynamicsInC* _this, struct UniformSequence cellPopulation, struct BlockStatisticsInC* statistics)
{
    ///
    double rhoBar;
    double jSeq[(*(*_this).descriptor_).d];
    struct UniformSequence_const j={(*(*_this).descriptor_).d, sizeof(double), jSeq};
    struct UniformSequence_const cellPops = {cellPopulation.num, cellPopulation.length, cellPopulation.sequence};
    struct UniformSequence j2={(*(*_this).descriptor_).d, sizeof(double), jSeq};
    get_rhoBar_j_InC(cellPops, &rhoBar, j2, (*_this).descriptor_);
    double uSqr =bgk_ma2_collision_InC(cellPopulation, rhoBar, j, (*_this).omega, (*_this).descriptor_);
    ///Internal statistics is not implemented yet
    /*if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }*/
}

void collideExternal_BGKDynamics(struct DynamicsInC* _this, struct UniformSequence cellPopulation, double rhoBar, struct UniformSequence_const j, struct BlockStatisticsInC* statistics, double thetaBar)
{
    ///
    bgk_ma2_collision_InC(cellPopulation, rhoBar, j, (*_this).omega, (*_this).descriptor_);
}

double computeEquilibrium_BGKDynamics (struct DynamicsInC* _this, int iPop, double rhoBar, struct UniformSequence_const j, double jSqr, double thetaBar)
{
    ///
    double invRho = (*(*_this).descriptor_).invRho(rhoBar);
    return bgk_ma2_equilibrium_InC(iPop, rhoBar, invRho, j, jSqr, (*_this).descriptor_);
}




