

#include "latticeBoltzmann_/basicRoutineInC.h"


//------------------------------- Geometric operations -------------------------
double normSqr_InC(struct UniformSequence_const vec_)
{
    const double* elems= (const double*) (vec_.sequence);
    double norm=0;
    for (int i=0; i<vec_.num; i++)
    {
        norm+=(*(elems+i))*(*(elems+i));
    }
    return norm;
}

//------------------------------- Physical fields ----------------------------

void get_rhoBar_j_InC(struct UniformSequence_const cellPopulation, double* rhoBar, struct UniformSequence j, struct DescriptorInC* descriptor_)
{
    const double* cellPop= (const double*) (cellPopulation.sequence);

    *rhoBar=*cellPop;
    for (int iPop=1; iPop < cellPopulation.num; ++iPop) {
        *rhoBar += *(cellPop+iPop);
    }

    double* j_= (double*) (j.sequence);
    for (int iD=0; iD < (*descriptor_).d; ++iD) {
        *(j_+iD) = (*cellPop)*(*descriptor_).c[0][iD];
    }
    for (int iPop=1; iPop < (*descriptor_).q; ++iPop) {
        for (int iD=0; iD < (*descriptor_).d; ++iD) {
           *(j_+iD) += (*(cellPop+iPop))*(*descriptor_).c[iPop][iD];
        }
    }
}
