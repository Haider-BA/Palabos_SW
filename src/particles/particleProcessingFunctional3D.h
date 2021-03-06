/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARTICLE_PROCESSING_FUNCTIONAL_3D_H
#define PARTICLE_PROCESSING_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "algorithm/functions.h"
#include <map>

namespace plb {

/// Count the number of particles, no matter which kind, found inside the domain.
template<typename T, template<typename U> class Descriptor>
class CountParticlesFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountParticlesFunctional3D<T,Descriptor>* clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint numParticlesId;
};

/// Count the number of particles, no matter which kind, found inside the domain.
template<typename T, template<typename U> class Descriptor>
class CountParticlesSelectiveFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountParticlesSelectiveFunctional3D(util::SelectInt* tags_);
    ~CountParticlesSelectiveFunctional3D();
    CountParticlesSelectiveFunctional3D(CountParticlesSelectiveFunctional3D<T,Descriptor> const& rhs);
    CountParticlesSelectiveFunctional3D<T,Descriptor>& operator=(CountParticlesSelectiveFunctional3D<T,Descriptor> const& rhs);
    void swap(CountParticlesSelectiveFunctional3D<T,Descriptor>& rhs);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountParticlesSelectiveFunctional3D<T,Descriptor>* clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint numParticlesId;
    util::SelectInt* tags;
};

/// Compute the average over all particle velocities.
template<typename T, template<typename U> class Descriptor>
class AverageParticleVelocityFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    AverageParticleVelocityFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual AverageParticleVelocityFunctional3D<T,Descriptor>* clone() const;
    Array<T,3> getAverageParticleVelocity() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Array<plint,3> averageVelocityId;
};

/// Copy particles of a certain tag from one field to another.
template<typename T, template<typename U> class Descriptor>
class CopySelectParticles3D : public BoxProcessingFunctional3D
{
public:
    CopySelectParticles3D(util::SelectInt* tags_);
    ~CopySelectParticles3D();
    CopySelectParticles3D(CopySelectParticles3D<T,Descriptor> const& rhs);
    CopySelectParticles3D<T,Descriptor>& operator=(CopySelectParticles3D<T,Descriptor> const& rhs);
    void swap(CopySelectParticles3D<T,Descriptor>& rhs);
    /// Arguments: [0] From Particle-field, [1] To Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CopySelectParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    util::SelectInt* tags;
};

/// Inject particles into the domain. The particles must be defined in a non-
///   parallel way, and duplicated over all processors.
template<typename T, template<typename U> class Descriptor>
class InjectParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    /// The particles are not consumed in this class. A clone of the particles is
    ///   automatically made as they are added into the domain.
    InjectParticlesFunctional3D(std::vector<Particle3D<T,Descriptor>*>& particles_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual InjectParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<Particle3D<T,Descriptor>*>& particles;
};

/// Generate a random number of particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell. All particles are identical clones (except for their position).
template<typename T, template<typename U> class Descriptor>
class InjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    InjectRandomParticlesFunctional3D(Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_);
    InjectRandomParticlesFunctional3D(InjectRandomParticlesFunctional3D<T,Descriptor> const& rhs);
    InjectRandomParticlesFunctional3D<T,Descriptor>&
        operator=(InjectRandomParticlesFunctional3D<T,Descriptor> const& rhs);
    void swap(InjectRandomParticlesFunctional3D<T,Descriptor>& rhs);
    ~InjectRandomParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual InjectRandomParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Particle3D<T,Descriptor>* particleTemplate;
    T probabilityPerCell;
};

/// Generate a random number of point-particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell.
template<typename T, template<typename U> class Descriptor, class DomainFunctional>
class AnalyticalInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    AnalyticalInjectRandomParticlesFunctional3D(Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_, DomainFunctional functional_);
    AnalyticalInjectRandomParticlesFunctional3D(AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs);
    AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>&
        operator=(AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs);
    void swap(AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>& rhs);
    ~AnalyticalInjectRandomParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Particle3D<T,Descriptor>* particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
};

/// Generate a random number of point-particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell. Additionally to analytically-inject, this functional uses a bit-
///   mask to decide where to inject.
template<typename T, template<typename U> class Descriptor, class DomainFunctional>
class MaskedInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    MaskedInjectRandomParticlesFunctional3D(Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_, DomainFunctional functional_, int flag_);
    MaskedInjectRandomParticlesFunctional3D(MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs);
    MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>&
        operator=(MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs);
    void swap(MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>& rhs);
    ~MaskedInjectRandomParticlesFunctional3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Particle3D<T,Descriptor>* particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
    int flag;
};

/// Generate equally spaced particles inside each cell. Every cell generates
///   nx particles in its x-direction, ny in its y-direction and nz in its
///   z-direction (nx * ny * nz in total in each cell).
///   All particles are identical clones (except for their position).
template<typename T, template<typename U> class Descriptor>
class InjectEquallySpacedParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    InjectEquallySpacedParticlesFunctional3D(Particle3D<T,Descriptor>* particleTemplate_, plint nx_, plint ny_, plint nz_);
    InjectEquallySpacedParticlesFunctional3D(InjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs);
    InjectEquallySpacedParticlesFunctional3D<T,Descriptor>&
        operator=(InjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs);
    void swap(InjectEquallySpacedParticlesFunctional3D<T,Descriptor>& rhs);
    ~InjectEquallySpacedParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual InjectEquallySpacedParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Particle3D<T,Descriptor>* particleTemplate;
    plint nx, ny, nz;
};

/// Generate equally spaced particles inside each cell. Every cell generates
///   nx particles in its x-direction, ny in its y-direction and nz in its
///   z-direction (nx * ny * nz in total in each cell). This functional uses
///   a bit-mask to decide where to inject. All particles are identical clones
///   (except for their position).
template<typename T, template<typename U> class Descriptor>
class MaskedInjectEquallySpacedParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    MaskedInjectEquallySpacedParticlesFunctional3D(Particle3D<T,Descriptor>* particleTemplate_,
            plint nx_, plint ny_, plint nz_, int flag_);
    MaskedInjectEquallySpacedParticlesFunctional3D(
            MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs);
    MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>&
        operator=(MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs);
    void swap(MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>& rhs);
    ~MaskedInjectEquallySpacedParticlesFunctional3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Particle3D<T,Descriptor>* particleTemplate;
    plint nx, ny, nz;
    int flag;
};

/// Remove all particles from a given domain.
template<typename T, template<typename U> class Descriptor>
class AbsorbParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    /// Argument: Particle-field.
    virtual AbsorbParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

/// Remove all particles from a given domain.
template<typename T, template<typename U> class Descriptor>
class AbsorbParticlesFunctionalSelective3D : public BoxProcessingFunctional3D
{
public:
    AbsorbParticlesFunctionalSelective3D(util::SelectInt* tags_);
    ~AbsorbParticlesFunctionalSelective3D();
    AbsorbParticlesFunctionalSelective3D(AbsorbParticlesFunctionalSelective3D<T,Descriptor> const& rhs);
    AbsorbParticlesFunctionalSelective3D<T,Descriptor>& operator=(AbsorbParticlesFunctionalSelective3D<T,Descriptor> const& rhs);
    void swap(AbsorbParticlesFunctionalSelective3D<T,Descriptor>& rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    /// Argument: Particle-field.
    virtual AbsorbParticlesFunctionalSelective3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    util::SelectInt* tags;
};

/// Find particles injected inside wall nodes and remove them.
template<typename T, template<typename U> class Descriptor>
class RemoveParticlesFromWall3D : public BoxProcessingFunctional3D
{
public:
    RemoveParticlesFromWall3D(int wallFlag_);
    /// Arguments: [0] Particle-field [1] Flag-matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual RemoveParticlesFromWall3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    int wallFlag; // Value that represents the wall in the flag matrix.
};

/// Find particles close to a wall and change their positions so they are pushed back to the flow field.
template<typename T, template<typename U> class Descriptor>
class PushParticlesAwayFromWall3D : public BoxProcessingFunctional3D
{
public:
    PushParticlesAwayFromWall3D(T cutOffValue_, T movingDistance_, int wallFlag_, int fluidFlag_);
    /// Arguments: [0] Particle-field [1] Flag-matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual PushParticlesAwayFromWall3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T cutOffValue; // When the speed of the particle drops below sqrt(cutOffValue), then this particle is a candidate for pushing.
    T movingDistance; // This is the distance the particles will be moved.
    int wallFlag; // Value that represents the wall nodes in the flag matrix.
    int fluidFlag; // Value that represents the fluid nodes in the flag matrix.
};

/// Execute the particle-fluid interaction step (during which the particles
///   don't move and the fluid doesn't change).
template<typename T, template<typename U> class Descriptor>
class FluidToParticleCoupling3D : public BoxProcessingFunctional3D
{
public:
    /// Particle speed = scaling*fluid speed.
    FluidToParticleCoupling3D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Fluid.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual FluidToParticleCoupling3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T scaling;
};

template<typename T, template<typename U> class Descriptor>
class VelocityToParticleCoupling3D : public BoxProcessingFunctional3D
{
public:
    /// Particle speed = scaling*fluid speed.
    VelocityToParticleCoupling3D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Velocity.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual VelocityToParticleCoupling3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T scaling;
};

template<typename T, template<typename U> class Descriptor>
class RhoBarJtoParticleCoupling3D : public BoxProcessingFunctional3D
{
public:
    /// Particle speed = scaling*fluid speed.
    RhoBarJtoParticleCoupling3D(bool velIsJ_, T scaling_);
    /// Arguments: [0] Particle-field; [1] rhoBarJ.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual RhoBarJtoParticleCoupling3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    bool velIsJ;
    T scaling;
};

/// Execute the iteration step during which particles advance.
template<typename T, template<typename U> class Descriptor>
class AdvanceParticlesFunctional3D : public BoxProcessingFunctional3D
{
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesFunctional3D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual AdvanceParticlesFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T cutOffValue;
};

/// Execute the iteration step during which particles advance, on the whole domain
/** The data processor's domain indication is being ignored. This works also with periodicity. **/
template<typename T, template<typename U> class Descriptor>
class AdvanceParticlesEveryWhereFunctional3D : public BoxProcessingFunctional3D
{
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesEveryWhereFunctional3D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T cutOffValue;
};


/* ******** VerletUpdateVelocity3D *********************************** */

/// Update the velocity to complete an iteration of the Verlet algorithm. Works
/// with Verlet particles only.
template<typename T, template<typename U> class Descriptor>
class VerletUpdateVelocity3D : public BoxProcessingFunctional3D
{
public:
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual VerletUpdateVelocity3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    bool projectForce;
    Array<T,3> planeNormal;
};


/* ******** VerletUpdateVelocitySelective3D *********************************** */

/// Update the velocity to complete an iteration of the Verlet algorithm. Works
/// with Verlet particles only. Acts only on particles with the specified tag.
template<typename T, template<typename U> class Descriptor>
class VerletUpdateVelocitySelective3D : public BoxProcessingFunctional3D
{
public:
    VerletUpdateVelocitySelective3D(util::SelectInt* tags_);
    ~VerletUpdateVelocitySelective3D();
    VerletUpdateVelocitySelective3D(VerletUpdateVelocitySelective3D<T,Descriptor> const& rhs);
    VerletUpdateVelocitySelective3D<T,Descriptor>& operator=(VerletUpdateVelocitySelective3D<T,Descriptor> const& rhs);
    void swap(VerletUpdateVelocitySelective3D<T,Descriptor>& rhs);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual VerletUpdateVelocitySelective3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    util::SelectInt* tags;
};

template< typename T,
          template<typename U> class Descriptor >
void addWallParticles (
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles, TriangleBoundary3D<T>& boundary );

template< typename T,
          template<typename U> class Descriptor, class ParticleFieldT >
void addWallParticlesGeneric (
    MultiParticleField3D<ParticleFieldT>& particles, TriangleBoundary3D<T>& boundary );

/// Count the number of particles at each cell node and add the result to the scalar field.
template<typename T, template<typename U> class Descriptor>
class CountAndAccumulateParticles3D : public BoxProcessingFunctional3D
{
public:
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountAndAccumulateParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

/// Count the number of particles with a given tag at each cell node and add the result to the scalar field.
template<typename T, template<typename U> class Descriptor>
class CountAndAccumulateTaggedParticles3D : public BoxProcessingFunctional3D
{
public:
    CountAndAccumulateTaggedParticles3D(plint tag_);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountAndAccumulateTaggedParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint tag;
};

/// Count the number of particles (with a given tag) at each refined cell node, and add the result
/// to the scalar field which is refined (defined on a refined grid with respect to the particle grid).
/// The particles which belong to each "sub-volume" of the refined scalar grid contained in the
/// "big volume" of the particle grid, must be identified, counted and accumulated.
template<typename T, template<typename U> class Descriptor>
class CountAndAccumulateTaggedParticlesRefined3D : public BoxProcessingFunctional3D
{
public:
    CountAndAccumulateTaggedParticlesRefined3D(plint tag_, plint dxScale_);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D coarseDomain, std::vector<AtomicBlock3D*> fields);
    virtual CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    /// These helper functions are also implemented in the ParticleField3D class,
    /// but we need to re-implement them here, since we need them for the refined
    /// scalar field.
    static plint nearestCell(T pos);
    static void computeGridPosition(Array<T,3> const& position, Dot3D const& location,
            plint& iX, plint& iY, plint& iZ);
private:
    plint tag;
    plint dxScale;
};

/// Count the number of particles with given tags at each cell node and place the result to the scalar field.
template<typename T, template<typename U> class Descriptor>
class CountTaggedParticles3D : public BoxProcessingFunctional3D
{
public:
    CountTaggedParticles3D(util::SelectInt* tags_);
    ~CountTaggedParticles3D();
    CountTaggedParticles3D(CountTaggedParticles3D<T,Descriptor> const& rhs);
    CountTaggedParticles3D<T,Descriptor>& operator=(CountTaggedParticles3D<T,Descriptor> const& rhs);
    void swap(CountTaggedParticles3D<T,Descriptor>& rhs);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual CountTaggedParticles3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    util::SelectInt* tags;
};

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint countParticles (
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain );

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint countParticles (
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, util::SelectInt* tags );

template<typename T, template<typename U> class Descriptor>
void injectParticles(std::vector<Particle3D<T,Descriptor>*>& injectedParticles,
                     MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles, Box3D domain);

/* Iterations of a passive-scalar fluid-particle system:
 * =====================================================
 *
 * Note: The difficulty comes from the fact that particle-fields may have a larger
 *   envelope than the fluid. When advancing particles on bulk and envelope, the
 *   velocity data from the fluid is therefore not necessarily locally available.
 *   The velocity is therefore first stored in the particle (in the bulk), and then
 *   communicated to the envelopes.
 *
 * --- Particles are at time t, fluid is at time t, defined on bulk and envelope. ---
 *  1. Fluid collideAndStream().
 *  2. Particle advance (bulk+envelope). ==> Particles at time t on bulk (needs no communication).
 *  3. Fluid communication ==> Fluid at time t+1.
 *  4. Particle interact (bulk domain) with velocity at time t+1.
 *  5. Particle communication ==> Particle at time t+1.
 */

}  // namespace plb

#endif  // PARTICLE_PROCESSING_FUNCTIONAL_3D_H
