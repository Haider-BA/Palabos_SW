#ifndef DYNAMICSPROCESSOR3DINC_H_INCLUDED
#define DYNAMICSPROCESSOR3DINC_H_INCLUDED


#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/dynamics.h"
#include "core_/dynamicsInC.h"

using namespace plb;

/* *************** Class ExternalRhoJcollideAndStream3DInC ******************* */

template<typename T, template<typename U> class Descriptor>
class ExternalRhoJcollideAndStream3DInC : public BoxProcessingFunctional3D
{
public:
    ExternalRhoJcollideAndStream3DInC(DynamicsInC&);
    // Block 0: lattice; Block 1: rhoBar; Block 2: j.
    virtual void processGenericBlocks( Box3D domain,
                                       std::vector<AtomicBlock3D*> atomicBlocks );
    virtual ExternalRhoJcollideAndStream3DInC<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    virtual void collide (
            BlockLattice3D<T,Descriptor>& lattice, Box3D const& domain,
            ScalarField3D<T> const& rhoBarField, Dot3D const& offset1,
            TensorField3D<T,3> const& jField, Dot3D const& offset2, BlockStatisticsInC& stat );
    virtual void bulkCollideAndStream (
            BlockLattice3D<T,Descriptor>& lattice, Box3D const& domain,
            ScalarField3D<T> const& rhoBarField, Dot3D const& offset1,
            TensorField3D<T,3> const& jField, Dot3D const& offset2, BlockStatisticsInC& stat );
    virtual void boundaryStream (
            BlockLattice3D<T,Descriptor>& lattice,
            Box3D const& bound, Box3D const& domain );
private:
    DynamicsInC& backGroundDym;
};



#endif // DYNAMICSPROCESSOR3DINC_H_INCLUDED
