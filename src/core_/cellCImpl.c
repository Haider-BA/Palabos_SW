

#include "core_/cellCImpl.h"

void cellRevert(struct UniformSequence cellPopulations) {

    double* cellPops = (double*) (cellPopulations.sequence);
    for (int iPop=1; iPop<=cellPopulations.num/2; ++iPop) {
        double tmp=*(cellPops+iPop);
        *(cellPops+iPop)=*(cellPops+iPop+cellPopulations.num/2);
        *(cellPops+iPop+cellPopulations.num/2)=tmp;
       // swap(*(cellPops+iPop),*(cellPops+iPop+cellPopulations.num/2));
    }
}
