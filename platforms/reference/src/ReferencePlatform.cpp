/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferencePlatform.h"
#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/Vec3.h"

using namespace OpenMM;
using namespace std;

ReferencePlatform::ReferencePlatform() {
    ReferenceKernelFactory* factory = new ReferenceKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(ApplyConstraintsKernel::Name(), factory);
    registerKernelFactory(VirtualSitesKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCMAPTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomHbondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCentroidBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCompoundBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCVForceKernel::Name(), factory);
    registerKernelFactory(CalcRMSDForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(CalcGayBerneForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateNoseHooverStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinMiddleStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
}

double ReferencePlatform::getSpeed() const {
    return 1;
}

bool ReferencePlatform::supportsDoublePrecision() const {
    return true;
}

void ReferencePlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    context.setPlatformData(new PlatformData(context.getSystem()));
}

void ReferencePlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

ReferencePlatform::PlatformData::PlatformData(const System& system) : time(0.0), stepCount(0), numParticles(system.getNumParticles()) {
    positions = new vector<Vec3>(numParticles);
    velocities = new vector<Vec3>(numParticles);
    forces = new vector<Vec3>(numParticles);

    forces_drl_bon = new vector<Vec3>(numParticles);
    forces_drl_ang = new vector<Vec3>(numParticles);
    forces_drl_tor = new vector<Vec3>(numParticles);
    forces_drl_n14 = new vector<Vec3>(numParticles);

    energies_drl_bon = new vector<vector<double>>();
    energies_drl_bon->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_bon->size(); ++Ix) {
        (*energies_drl_bon)[Ix].resize(numParticles);
    }

    energies_drl_ang = new vector<vector<double>>();
    energies_drl_ang->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_ang->size(); ++Ix) {
        (*energies_drl_ang)[Ix].resize(numParticles);
    }

    energies_drl_tor = new vector<vector<double>>();
    energies_drl_tor->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_tor->size(); ++Ix) {
        (*energies_drl_tor)[Ix].resize(numParticles);        
    }

    energies_drl_n14 = new vector<vector<double>>();
    energies_drl_n14->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_n14->size(); ++Ix) {
        (*energies_drl_n14)[Ix].resize(numParticles);
    }

    energies_drl_vdw = new vector<vector<double>>();
    energies_drl_vdw->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_vdw->size(); ++Ix) {
        (*energies_drl_vdw)[Ix].resize(numParticles);
    }

    energies_drl_cou = new vector<vector<double>>();
    energies_drl_cou->resize(numParticles);
    for (size_t Ix = 0; Ix < energies_drl_cou->size(); ++Ix) {
        (*energies_drl_cou)[Ix].resize(numParticles);
    }

    periodicBoxSize = new Vec3();
    periodicBoxVectors = new Vec3[3];
    constraints = new ReferenceConstraints(system);
    energyParameterDerivatives = new map<string, double>();
}

ReferencePlatform::PlatformData::~PlatformData() {
    delete positions;
    delete velocities;
    delete forces;

    delete forces_drl_bon;
    delete forces_drl_ang;
    delete forces_drl_tor;
    delete forces_drl_n14;

    delete energies_drl_bon;
    delete energies_drl_ang;
    delete energies_drl_tor;
    delete energies_drl_n14;

    delete energies_drl_vdw;
    delete energies_drl_cou;

    delete periodicBoxSize;
    delete[] periodicBoxVectors;
    delete constraints;
    delete energyParameterDerivatives;
}
