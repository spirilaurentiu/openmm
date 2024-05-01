/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2020 Stanford University and the Authors.      *
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

#include "openmm/serialization/StateProxy.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include <map>

using namespace std;
using namespace OpenMM;

StateProxy::StateProxy() : SerializationProxy("State") {

}

void StateProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    node.setStringProperty("openmmVersion", Platform::getOpenMMVersion());
    const State& s = *reinterpret_cast<const State*>(object);
    node.setDoubleProperty("time", s.getTime());
    Vec3 a,b,c;
    s.getPeriodicBoxVectors(a,b,c);
    SerializationNode& boxVectorsNode = node.createChildNode("PeriodicBoxVectors");
    boxVectorsNode.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    boxVectorsNode.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    boxVectorsNode.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    if ((s.getDataTypes()&State::Parameters) != 0) {
        s.getParameters();
        SerializationNode& parametersNode = node.createChildNode("Parameters");
        for (auto& param : s.getParameters())
            parametersNode.setDoubleProperty(param.first, param.second);
    }
    if ((s.getDataTypes()&State::Energy) != 0) {
        s.getPotentialEnergy();
        SerializationNode& energiesNode = node.createChildNode("Energies");
        energiesNode.setDoubleProperty("PotentialEnergy", s.getPotentialEnergy());
        energiesNode.setDoubleProperty("KineticEnergy", s.getKineticEnergy());
    }
    if ((s.getDataTypes()&State::Positions) != 0) {
        s.getPositions();
        SerializationNode& positionsNode = node.createChildNode("Positions");
        vector<Vec3> statePositions = s.getPositions();
        for (int i=0; i<statePositions.size();i++) {
           positionsNode.createChildNode("Position").setDoubleProperty("x", statePositions[i][0]).setDoubleProperty("y", statePositions[i][1]).setDoubleProperty("z", statePositions[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Velocities) != 0) {
        s.getVelocities();
        SerializationNode& velocitiesNode = node.createChildNode("Velocities");
        vector<Vec3> stateVelocities = s.getVelocities();
        for (int i=0; i<stateVelocities.size();i++) {
           velocitiesNode.createChildNode("Velocity").setDoubleProperty("x", stateVelocities[i][0]).setDoubleProperty("y", stateVelocities[i][1]).setDoubleProperty("z", stateVelocities[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Forces) != 0) {
        s.getForces();
        SerializationNode& forcesNode = node.createChildNode("Forces");
        vector<Vec3> stateForces = s.getForces();
        for (int i=0; i<stateForces.size();i++) {
            forcesNode.createChildNode("Force").setDoubleProperty("x", stateForces[i][0]).setDoubleProperty("y", stateForces[i][1]).setDoubleProperty("z", stateForces[i][2]);
        }
    }
    if ((s.getDataTypes()&State::IntegratorParameters) != 0) {
        node.getChildren().push_back(s.getIntegratorParameters());
    }

    // OPENMM DRILL drl

    if ((s.getDataTypes()&State::Forces_drl_bon) != 0) {
        s.getForces_drl_bon();
        SerializationNode& forcesNode_drl_bon = node.createChildNode("Forces_drl_bon");
        vector<Vec3> stateForces_drl_bon = s.getForces_drl_bon();
        for (int i=0; i<stateForces_drl_bon.size();i++) {
            forcesNode_drl_bon.createChildNode("Force").setDoubleProperty("x", stateForces_drl_bon[i][0]).setDoubleProperty("y", stateForces_drl_bon[i][1]).setDoubleProperty("z", stateForces_drl_bon[i][2]);
        }

        s.getEnergies_drl_bon(); // drl serialize 2D vector BEGIN
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_bon");
        vector<vector<double>> stateEnergies_drl_bon = s.getEnergies_drl_bon();
        for (int i=0; i<stateEnergies_drl_bon.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_bon[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_bon[i][j]);
            }
        } // drl serialize 2D vector END
    }

    if ((s.getDataTypes()&State::Forces_drl_ang) != 0) {
        s.getForces_drl_ang();
        SerializationNode& forcesNode_drl_ang = node.createChildNode("Forces_drl_ang");
        vector<Vec3> stateForces_drl_ang = s.getForces_drl_ang();
        for (int i=0; i<stateForces_drl_ang.size();i++) {
            forcesNode_drl_ang.createChildNode("Force").setDoubleProperty("x", stateForces_drl_ang[i][0]).setDoubleProperty("y", stateForces_drl_ang[i][1]).setDoubleProperty("z", stateForces_drl_ang[i][2]);
        }

        s.getEnergies_drl_ang(); // drl serialize 2D vector BEGIN
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_ang");
        vector<vector<double>> stateEnergies_drl_ang = s.getEnergies_drl_ang();
        for (int i=0; i<stateEnergies_drl_ang.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_ang[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_ang[i][j]);
            }
        } // drl serialize 2D vector END
    }

    if ((s.getDataTypes()&State::Forces_drl_tor) != 0) {
        s.getForces_drl_tor();
        SerializationNode& forcesNode_drl_tor = node.createChildNode("Forces_drl_tor");
        vector<Vec3> stateForces_drl_tor = s.getForces_drl_tor();
        for (int i=0; i<stateForces_drl_tor.size();i++) {
            forcesNode_drl_tor.createChildNode("Force").setDoubleProperty("x", stateForces_drl_tor[i][0]).setDoubleProperty("y", stateForces_drl_tor[i][1]).setDoubleProperty("z", stateForces_drl_tor[i][2]);
        }

        s.getEnergies_drl_tor(); // drl serialize 2D vector BEGIN
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_tor");
        vector<vector<double>> stateEnergies_drl_tor = s.getEnergies_drl_tor();
        for (int i=0; i<stateEnergies_drl_tor.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_tor[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_tor[i][j]);
            }
        } // drl serialize 2D vector END
    }

    if ((s.getDataTypes()&State::Forces_drl_n14) != 0) {
        s.getForces_drl_n14();
        SerializationNode& forcesNode_drl_n14 = node.createChildNode("Forces_drl_n14");
        vector<Vec3> stateForces_drl_n14 = s.getForces_drl_n14();
        for (int i=0; i<stateForces_drl_n14.size();i++) {
            forcesNode_drl_n14.createChildNode("Force").setDoubleProperty("x", stateForces_drl_n14[i][0]).setDoubleProperty("y", stateForces_drl_n14[i][1]).setDoubleProperty("z", stateForces_drl_n14[i][2]);
        }

        // drl serialize 2D vector BEGIN
        s.getEnergies_drl_n14();
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_n14");
        vector<vector<double>> stateEnergies_drl_n14 = s.getEnergies_drl_n14();
        for (int i=0; i<stateEnergies_drl_n14.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_n14[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_n14[i][j]);
            }
        } // drl serialize 2D vector END
    }

    if ((s.getDataTypes()&State::Forces_drl_vdw) != 0) {
        // s.getForces_drl_vdw();
        // SerializationNode& forcesNode_drl_vdw = node.createChildNode("Forces_drl_vdw");
        // vector<Vec3> stateForces_drl_vdw = s.getForces_drl_vdw();
        // for (int i=0; i<stateForces_drl_vdw.size();i++) {
        //     forcesNode_drl_vdw.createChildNode("Force").setDoubleProperty("x", stateForces_drl_vdw[i][0]).setDoubleProperty("y", stateForces_drl_vdw[i][1]).setDoubleProperty("z", stateForces_drl_vdw[i][2]);
        // }

        // drl serialize 2D vector BEGIN
        s.getEnergies_drl_vdw();
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_vdw");
        vector<vector<double>> stateEnergies_drl_vdw = s.getEnergies_drl_vdw();
        for (int i=0; i<stateEnergies_drl_vdw.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_vdw[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_vdw[i][j]);
            }
        } // drl serialize 2D vector END
    }

    if ((s.getDataTypes()&State::Forces_drl_cou) != 0) {
        // s.getForces_drl_cou();
        // SerializationNode& forcesNode_drl_cou = node.createChildNode("Forces_drl_cou");
        // vector<Vec3> stateForces_drl_cou = s.getForces_drl_cou();
        // for (int i=0; i<stateForces_drl_cou.size();i++) {
        //     forcesNode_drl_cou.createChildNode("Force").setDoubleProperty("x", stateForces_drl_cou[i][0]).setDoubleProperty("y", stateForces_drl_cou[i][1]).setDoubleProperty("z", stateForces_drl_cou[i][2]);
        // }

        // drl serialize 2D vector BEGIN
        s.getEnergies_drl_cou();
        SerializationNode& energiesNode_drl_bon = node.createChildNode("Energies_drl_cou");
        vector<vector<double>> stateEnergies_drl_cou = s.getEnergies_drl_cou();
        for (int i=0; i<stateEnergies_drl_cou.size();i++) {
            SerializationNode& innerNode = energiesNode_drl_bon.createChildNode("EnergyVector");
            for (size_t j = 0; j < stateEnergies_drl_cou[i].size(); j++) {
                innerNode.createChildNode("Element").setDoubleProperty("value", stateEnergies_drl_cou[i][j]);
            }
        } // drl serialize 2D vector END
    }

    // drl END

}

void* StateProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    double outTime = node.getDoubleProperty("time");
    const SerializationNode& boxVectorsNode = node.getChildNode("PeriodicBoxVectors");
    const SerializationNode& AVec = boxVectorsNode.getChildNode("A");
    Vec3 outAVec(AVec.getDoubleProperty("x"),AVec.getDoubleProperty("y"),AVec.getDoubleProperty("z"));
    const SerializationNode& BVec = boxVectorsNode.getChildNode("B");
    Vec3 outBVec(BVec.getDoubleProperty("x"),BVec.getDoubleProperty("y"),BVec.getDoubleProperty("z"));
    const SerializationNode& CVec = boxVectorsNode.getChildNode("C");
    Vec3 outCVec(CVec.getDoubleProperty("x"),CVec.getDoubleProperty("y"),CVec.getDoubleProperty("z"));
    int types = 0;
    vector<int> arraySizes;
    State::StateBuilder builder(outTime);
    for (auto& child : node.getChildren()) {
        if (child.getName() == "Parameters") {
            map<string, double> outStateParams;
            for (auto& param : child.getProperties())
                outStateParams[param.first] = child.getDoubleProperty(param.first);
            builder.setParameters(outStateParams);
        }
        else if (child.getName() == "Energies") {
            double potentialEnergy = child.getDoubleProperty("PotentialEnergy");
            double kineticEnergy = child.getDoubleProperty("KineticEnergy");
            builder.setEnergy(kineticEnergy, potentialEnergy);
        }
        else if (child.getName() == "Positions") {
            vector<Vec3> outPositions;
            for (auto& particle : child.getChildren())
                outPositions.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setPositions(outPositions);
            arraySizes.push_back(outPositions.size());
        }
        else if (child.getName() == "Velocities") {
            vector<Vec3> outVelocities;
            for (auto& particle : child.getChildren())
                outVelocities.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setVelocities(outVelocities);
            arraySizes.push_back(outVelocities.size());
        }
        else if (child.getName() == "Forces") {
            vector<Vec3> outForces;
            for (auto& particle : child.getChildren())
                outForces.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces(outForces);
            arraySizes.push_back(outForces.size());
        }
        else if (child.getName() == "IntegratorParameters") {
            builder.updateIntegratorParameters() = child;
        }
        
        else if (child.getName() == "Forces_drl_bon") {
            vector<Vec3> outForces_drl_bon;
            for (auto& particle : child.getChildren())
                outForces_drl_bon.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces_drl_bon(outForces_drl_bon);
            arraySizes.push_back(outForces_drl_bon.size());
        }
        else if (child.getName() == "Forces_drl_ang") {
            vector<Vec3> outForces_drl_ang;
            for (auto& particle : child.getChildren())
                outForces_drl_ang.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces_drl_ang(outForces_drl_ang);
            arraySizes.push_back(outForces_drl_ang.size());
        }
        else if (child.getName() == "Forces_drl_tor") {
            vector<Vec3> outForces_drl_tor;
            for (auto& particle : child.getChildren())
                outForces_drl_tor.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces_drl_tor(outForces_drl_tor);
            arraySizes.push_back(outForces_drl_tor.size());
        }
        else if (child.getName() == "Forces_drl_n14") {
            vector<Vec3> outForces_drl_n14;
            for (auto& particle : child.getChildren())
                outForces_drl_n14.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces_drl_n14(outForces_drl_n14);
            arraySizes.push_back(outForces_drl_n14.size());
        }
        // else if (child.getName() == "Forces_drl_vdw") {
        //     vector<Vec3> outForces_drl_vdw;
        //     for (auto& particle : child.getChildren())
        //         outForces_drl_vdw.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
        //     builder.setForces_drl_vdw(outForces_drl_vdw);
        //     arraySizes.push_back(outForces_drl_vdw.size());
        // }
        //         else if (child.getName() == "Forces_drl_cou") {
        //     vector<Vec3> outForces_drl_cou;
        //     for (auto& particle : child.getChildren())
        //         outForces_drl_cou.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
        //     builder.setForces_drl_cou(outForces_drl_cou);
        //     arraySizes.push_back(outForces_drl_cou.size());
        // }

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_bon") {
            vector<vector<double>> outEnergies_drl_bon;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_bon.push_back(innerVector);
            }           
            builder.setEnergies_drl_bon(outEnergies_drl_bon);
            arraySizes.push_back(outEnergies_drl_bon.size());
        }  
        // drl serialize 2D vector END

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_ang") {
            vector<vector<double>> outEnergies_drl_ang;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_ang.push_back(innerVector);
            }           
            builder.setEnergies_drl_ang(outEnergies_drl_ang);
            arraySizes.push_back(outEnergies_drl_ang.size());
        }  
        // drl serialize 2D vector END

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_tor") {
            vector<vector<double>> outEnergies_drl_tor;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_tor.push_back(innerVector);
            }           
            builder.setEnergies_drl_tor(outEnergies_drl_tor);
            arraySizes.push_back(outEnergies_drl_tor.size());
        }  
        // drl serialize 2D vector END

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_n14") {
            vector<vector<double>> outEnergies_drl_n14;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_n14.push_back(innerVector);
            }           
            builder.setEnergies_drl_n14(outEnergies_drl_n14);
            arraySizes.push_back(outEnergies_drl_n14.size());
        } 

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_vdw") {
            vector<vector<double>> outEnergies_drl_vdw;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_vdw.push_back(innerVector);
            }           
            builder.setEnergies_drl_vdw(outEnergies_drl_vdw);
            arraySizes.push_back(outEnergies_drl_vdw.size());
        } 

        // drl serialize 2D vector BEGIN
        else if (child.getName() == "Energies_drl_cou") {
            vector<vector<double>> outEnergies_drl_cou;
            for (const auto& innerNode : child.getChildren()) {
                vector<double> innerVector;                
                for (const auto& element : innerNode.getChildren()) {
                    double value = element.getDoubleProperty("value");
                    innerVector.push_back(value);
                }
                outEnergies_drl_cou.push_back(innerVector);
            }           
            builder.setEnergies_drl_cou(outEnergies_drl_cou);
            arraySizes.push_back(outEnergies_drl_cou.size());
        } 

        // drl serialize 2D vector END


    }
    for (int i = 1; i < arraySizes.size(); i++) {
        if (arraySizes[i] != arraySizes[i-1]) {
            throw(OpenMMException("State Deserialization Particle Size Mismatch, check number of particles in Forces, Velocities, Positions!"));
        }
    }
    builder.setPeriodicBoxVectors(outAVec, outBVec, outCVec);
    State *s = new State();
    *s = builder.getState();
    return s;
}