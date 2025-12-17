#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

/// @brief Gestion des runs avec histogrammes et statistiques
///
/// *** VERSION CORRIGÉE : Cône 60° avec normalisation correcte ***
///
/// NORMALISATION :
/// - Temps simulé : t = N_events / A_4π (TOUJOURS)
/// - Correction géométrique : f_corr = Ω_det / Ω_cône
/// - Débit réel : Ḋ_réel = Ḋ_brut × f_corr

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    
    void AddKermaEnergy(G4double edep);
    void AddKermaFluence(G4double fluenceContrib);
    void AddKermaEnergyForced(G4double forcedDeposit);
    
    void RecordEventStatistics(G4int nPrimaries, 
                               const std::vector<G4double>& primaryEnergies,
                               G4int nTransmitted,
                               G4int nAbsorbed,
                               G4double kermaDeposit);

private:
    // Données pour le calcul du kerma
    G4double fKermaTotalEnergy;
    G4double fKermaTotalEnergy2;
    G4double fKermaMass;
    G4double fKermaRadius;
    G4double fKermaPosition;
    G4int fKermaEventCount;
    
    // Méthode 2 (fluence)
    G4double fKermaFluenceTotal;
    G4int fKermaFluenceCount;
    
    // Méthode 1bis (forçage)
    G4double fKermaForcedTotal;
    G4int fKermaForcedCount;
    
    // Paramètres source
    G4double fActivity4pi;          // Activité 4π (Bq)
    G4double fConeAngle;            // Angle du cône d'émission (60°)
    
    // Paramètres géométriques
    G4double fSourcePosZ;
    G4double fDetectorPosZ;
    G4double fMeanGammasPerDecay;

    // Compteurs globaux
    G4int fTotalPrimariesGenerated;
    G4int fTotalEventsWithZeroGamma;
    G4int fTotalTransmitted;
    G4int fTotalAbsorbed;
    
    G4String fOutputFileName;
};

#endif
