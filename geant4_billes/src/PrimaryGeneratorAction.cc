#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════
// VERSION : ÉMISSION DANS UN CÔNE DE 60°
// ═══════════════════════════════════════════════════════════════════════════
// L'émission se fait dans un cône de 60° (demi-angle au sommet)
// La normalisation temporelle est t = N / A_4π
// Le facteur de correction est f_corr = f_cone = (1 - cos(60°))/2 = 0.25
// ═══════════════════════════════════════════════════════════════════════════

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fLastEventGammaCount(0),
  fConeAngle(60.0 * deg),   // *** CÔNE DE 60° ***
  fSourcePosition(0.*mm, 0.*mm, 2.*cm)
{
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    
    // ═══════════════════════════════════════════════════════════════
    // DÉFINITION DU TYPE DE PARTICULE (GAMMA)
    // ═══════════════════════════════════════════════════════════════
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

    // ═══════════════════════════════════════════════════════════════
    // SPECTRE GAMMA DE L'EUROPIUM-152
    // ═══════════════════════════════════════════════════════════════
    fGammaEnergies = {
        40.12, 39.52, 121.78, 244.70, 344.28, 411.12,
        443.96, 778.90, 867.38, 964.08, 1112.07, 1408.01
    };

    fGammaIntensities = {
        37.7, 20.8, 28.5, 7.6, 26.5, 2.2,
        2.8, 12.9, 4.2, 14.6, 13.6, 21.0
    };

    // ═══════════════════════════════════════════════════════════════
    // CONVERSION EN PROBABILITÉS
    // ═══════════════════════════════════════════════════════════════
    fGammaProbabilities.clear();
    G4double totalIntensity = 0.0;
    for (size_t i = 0; i < fGammaIntensities.size(); ++i) {
        G4double prob = fGammaIntensities[i] / 100.0;
        fGammaProbabilities.push_back(prob);
        totalIntensity += fGammaIntensities[i];
    }

    // ═══════════════════════════════════════════════════════════════
    // CALCUL DU FACTEUR DE CORRECTION
    // ═══════════════════════════════════════════════════════════════
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;

    // ═══════════════════════════════════════════════════════════════
    // AFFICHAGE DU SPECTRE
    // ═══════════════════════════════════════════════════════════════
    G4cout << "\n╔═══════════════════════════════════════════════════════════════════╗\n";
    G4cout << "║        SPECTRE GAMMA EUROPIUM-152 (Eu-152)                        ║\n";
    G4cout << "║        *** ÉMISSION DANS UN CÔNE DE 60° ***                       ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  Énergie (keV)  │  Intensité (%)  │ Prob. d'émission              ║\n";
    G4cout << "╠═════════════════╪═════════════════╪═══════════════════════════════╣\n";

    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {
        char buffer[120];
        sprintf(buffer, "║    %7.2f      │      %5.1f      │      %6.4f                   ║",
                fGammaEnergies[i], fGammaIntensities[i], fGammaProbabilities[i]);
        G4cout << buffer << G4endl;
    }

    G4cout << "╠═══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Intensité totale : " << totalIntensity << " %                                     ║" << G4endl;
    G4cout << "║  Nombre moyen de gammas/désintégration : " << totalIntensity/100.0 << "                   ║" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  MÉTHODE : Tirage indépendant pour chaque raie.                   ║" << G4endl;
    G4cout << "║            Émission dans un CÔNE de " << fConeAngle/deg << "°                        ║" << G4endl;
    G4cout << "║  Fraction angle solide (Ω_cône/4π) : " << f_cone*100 << "%                      ║" << G4endl;
    G4cout << "║  Facteur de correction f_corr = " << f_cone << "                          ║" << G4endl;
    G4cout << "╚═══════════════════════════════════════════════════════════════════╝\n" << G4endl;

    fParticleGun->SetParticlePosition(fSourcePosition);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GenerateDirectionInCone(G4double coneAngle,
                                                     G4double& theta,
                                                     G4double& phi,
                                                     G4ThreeVector& direction)
{
    // ═══════════════════════════════════════════════════════════════════════
    // GÉNÉRATION DANS UN CÔNE D'ANGLE coneAngle
    // ═══════════════════════════════════════════════════════════════════════
    // Distribution uniforme dans le cône : cos(θ) uniforme sur [cos(coneAngle), 1]
    // Le cône est orienté vers +z (vers le détecteur)
    // ═══════════════════════════════════════════════════════════════════════

    G4double cosConeAngle = std::cos(coneAngle);
    G4double cosTheta = cosConeAngle + (1.0 - cosConeAngle) * G4UniformRand();
    theta = std::acos(cosTheta);
    phi = G4UniformRand() * 2.0 * M_PI;

    G4double sinTheta = std::sin(theta);
    G4double u_x = sinTheta * std::cos(phi);
    G4double u_y = sinTheta * std::sin(phi);
    G4double u_z = cosTheta;

    direction = G4ThreeVector(u_x, u_y, u_z);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // ═══════════════════════════════════════════════════════════════════════
    // SIMULATION RÉALISTE D'UNE DÉSINTÉGRATION Eu-152
    // ═══════════════════════════════════════════════════════════════════════

    fLastEventGammaCount = 0;

    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {

        G4double random = G4UniformRand();

        if (random < fGammaProbabilities[i]) {

            G4double energy = fGammaEnergies[i] * keV;
            fParticleGun->SetParticleEnergy(energy);

            // Direction dans le cône de 60°
            G4double theta, phi;
            G4ThreeVector direction;
            GenerateDirectionInCone(fConeAngle, theta, phi, direction);
            fParticleGun->SetParticleMomentumDirection(direction);

            fParticleGun->GeneratePrimaryVertex(anEvent);
            fLastEventGammaCount++;
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // AFFICHAGE DE DIAGNOSTIC
    // ═══════════════════════════════════════════════════════════════════════
    G4int eventID = anEvent->GetEventID();
    if (eventID < 5 || eventID % 100000 == 0) {
        G4cout << "PrimaryGenerator | Event " << eventID
               << " : " << fLastEventGammaCount << " gamma(s) generated (cone 60°)";

        if (fLastEventGammaCount > 0 && eventID < 5) {
            G4cout << " [";
            for (G4int v = 0; v < anEvent->GetNumberOfPrimaryVertex(); ++v) {
                G4PrimaryVertex* vertex = anEvent->GetPrimaryVertex(v);
                if (vertex && vertex->GetPrimary()) {
                    if (v > 0) G4cout << ", ";
                    G4cout << vertex->GetPrimary()->GetKineticEnergy()/keV << " keV";
                }
            }
            G4cout << "]";
        }
        G4cout << G4endl;
    }
}
