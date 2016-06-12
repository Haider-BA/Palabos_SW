
#include "latticeBoltzmann_/basicRoutineInC.h"
#include "core_/dynamicsInC.h"

//----------------------- Dynamics -------------------------------
/*///This function did nothing
void computeRhoBarJ_Dynamics (UniformSequence_const cellPopulation, double& rhoBar, UniformSequence j)
{}

///This function did nothing
void collide_Dynamics (UniformSequence cellPopulation, BlockStatisticsInC& statistics);
{}

///This function did nothing
void computeEquilibrium_Dynamics (int iPop, double rhoBar, UniformSequence_const j, double jSqr, double thetaBar=0)
{}*/

void collideExternal_Dynanics(struct DynamicsInC* _this, struct UniformSequence cellPopulation, double rhoBar,
                              struct UniformSequence_const j, struct BlockStatisticsInC* statistics, double thetaBar)
{
    double oldRhoBar;
    double oldJSeq[j.num];
    struct UniformSequence oldJ={j.num, sizeof(double), oldJSeq};
    struct UniformSequence_const cellPop2={cellPopulation.num, cellPopulation.length, cellPopulation.sequence};
    (*_this).computeRhoBarJ(_this, cellPop2, &oldRhoBar, oldJ);
    struct UniformSequence_const oldJ2={j.num, sizeof(double), oldJSeq};
    double oldJsqr = normSqr_InC(oldJ2);
    double jSqr = normSqr_InC(j);
    double* cellPop= (double*) (cellPopulation.sequence);
    for (int iPop=0; iPop<cellPopulation.num; ++iPop) {
        double oldEq = (*_this).computeEquilibrium(_this, iPop, oldRhoBar, oldJ2, oldJsqr, thetaBar);
        double newEq = (*_this).computeEquilibrium(_this, iPop, rhoBar, j, jSqr, thetaBar);
        *(cellPop+iPop)= newEq - oldEq;
    }
    (*_this).collide(_this, cellPopulation, statistics);
}

void computeEquilibria_Dynamics (struct DynamicsInC* _this, struct UniformSequence fEq, double rhoBar,
                                 struct UniformSequence_const j, double jSqr, double thetaBar)
{

    double* cellPop= (double*) (fEq.sequence);
    for (int iPop=0; iPop<fEq.num; ++iPop) {
        *(cellPop+iPop)= (*_this).computeEquilibrium(_this, iPop, rhoBar, j, jSqr, thetaBar);
    }
}

//--------------------------- BasicBulkDynamics --------------------------------

void computeRhoBarJ_BasicBulkDynamics (struct DynamicsInC* _this, struct UniformSequence_const cellPopulation,
                                       double* rhoBar, struct UniformSequence j)
{
    get_rhoBar_j_InC(cellPopulation, rhoBar, j, (*_this).descriptor_);
}
