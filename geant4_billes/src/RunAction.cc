#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════
// VERSION : CÔNE 60° AVEC NORMALISATION CORRECTE
// ═══════════════════════════════════════════════════════════════════════════
//
// PRINCIPE FONDAMENTAL :
// ----------------------
// 1 événement Geant4 = 1 désintégration du noyau Eu-152
// Le temps simulé est TOUJOURS : t = N_events / A_4π
//
// CORRECTION GÉOMÉTRIQUE :
// ------------------------
// Les gammas sont émis dans un cône de 60° au lieu de 4π.
// On CONCENTRE donc tous les gammas dans ce cône.
// Le débit brut est surestimé d'un facteur 4π / Ω_cône = 1 / f_cone
//
// Pour un cône de 60° :
//     f_cone = (1 - cos(60°)) / 2 = (1 - 0.5) / 2 = 0.25
//
// Donc :
//     Ḋ_réel = Ḋ_brut × f_cone = Ḋ_brut × 0.25
//
// ═══════════════════════════════════════════════════════════════════════════

RunAction::RunAction()
: G4UserRunAction(),
  fKermaTotalEnergy(0.),
  fKermaTotalEnergy2(0.),
  fKermaMass(0.),
  fKermaRadius(2.0*cm),
  fKermaPosition(20.0*cm),
  fKermaEventCount(0),
  fKermaFluenceTotal(0.),
  fKermaFluenceCount(0),
  fKermaForcedTotal(0.),
  fKermaForcedCount(0),
  fActivity4pi(44000.0),              // Activité totale 4π (Bq)
  fConeAngle(60.0*deg),               // *** CÔNE D'ÉMISSION DE 60° ***
  fSourcePosZ(2.0*cm),
  fDetectorPosZ(20.0*cm),
  fMeanGammasPerDecay(1.924),
  fTotalPrimariesGenerated(0),
  fTotalEventsWithZeroGamma(0),
  fTotalTransmitted(0),
  fTotalAbsorbed(0),
  fOutputFileName("dose_water_output")
{
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);

    // Calcul du facteur de correction pour affichage
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;

    G4cout << "\n╔════════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║  RunAction initialized - VERSION CÔNE 60°                          ║" << G4endl;
    G4cout << "║  Output file: " << fOutputFileName << ".root                              ║" << G4endl;
    G4cout << "║  Normalisation: t = N_events / A_4π (CORRECTE)                     ║" << G4endl;
    G4cout << "║  Angle du cône: " << fConeAngle/deg << "°                                            ║" << G4endl;
    G4cout << "║  Facteur correction: f_corr = f_cone = " << f_cone << "                     ║" << G4endl;
    G4cout << "╚════════════════════════════════════════════════════════════════════╝\n" << G4endl;
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

    // Reset des compteurs
    fKermaTotalEnergy = 0.;
    fKermaTotalEnergy2 = 0.;
    fKermaEventCount = 0;
    fKermaFluenceTotal = 0.;
    fKermaFluenceCount = 0;
    fKermaForcedTotal = 0.;
    fKermaForcedCount = 0;
    fTotalPrimariesGenerated = 0;
    fTotalEventsWithZeroGamma = 0;
    fTotalTransmitted = 0;
    fTotalAbsorbed = 0;

    // Calcul de la masse du détecteur d'eau
    G4double kerma_volume = (4.0/3.0) * M_PI * std::pow(fKermaRadius, 3);
    G4double water_density = 1.0 * g/cm3;
    fKermaMass = kerma_volume * water_density;

    // Création des histogrammes
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fOutputFileName);

    analysisManager->CreateH1("nGammasPerEvent",
                              "Number of primary gammas per event;N_{#gamma};Counts",
                              15, -0.5, 14.5);
    analysisManager->CreateH1("energySpectrum",
                              "Energy spectrum of generated gammas;E (keV);Counts",
                              1500, 0., 1500.);
    analysisManager->CreateH1("totalEnergyPerEvent",
                              "Total primary energy per event;E_{tot} (keV);Counts",
                              500, 0., 5000.);
    analysisManager->CreateH1("nTransmittedPerEvent",
                              "Number of transmitted gammas per event;N_{trans};Counts",
                              15, -0.5, 14.5);
    analysisManager->CreateH1("nAbsorbedPerEvent",
                              "Number of absorbed gammas per event;N_{abs};Counts",
                              15, -0.5, 14.5);
    analysisManager->CreateH1("dosePerEvent",
                              "Dose energy deposit per event;E_{dose} (keV);Counts",
                              200, 0., 100.);
    analysisManager->CreateH2("nGammas_vs_totalEnergy",
                              "Number of gammas vs Total energy;N_{#gamma};E_{tot} (keV)",
                              15, -0.5, 14.5, 100, 0., 5000.);

    // Ntuples
    analysisManager->CreateNtuple("EventData", "Event-level data");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("nPrimaries");
    analysisManager->CreateNtupleDColumn("totalEnergy");
    analysisManager->CreateNtupleIColumn("nTransmitted");
    analysisManager->CreateNtupleIColumn("nAbsorbed");
    analysisManager->CreateNtupleIColumn("nScattered");
    analysisManager->CreateNtupleIColumn("nSecondaries");
    analysisManager->CreateNtupleDColumn("doseDeposit");
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("GammaData", "Primary gamma data");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("gammaIndex");
    analysisManager->CreateNtupleDColumn("energyInitial");
    analysisManager->CreateNtupleDColumn("energyUpstream");
    analysisManager->CreateNtupleDColumn("energyDownstream");
    analysisManager->CreateNtupleDColumn("theta");
    analysisManager->CreateNtupleDColumn("phi");
    analysisManager->CreateNtupleIColumn("detectedUpstream");
    analysisManager->CreateNtupleIColumn("detectedDownstream");
    analysisManager->CreateNtupleIColumn("transmitted");
    analysisManager->FinishNtuple();

    // ═══════════════════════════════════════════════════════════════════════
    // CALCUL DES FACTEURS GÉOMÉTRIQUES
    // ═══════════════════════════════════════════════════════════════════════
    G4double distance = fDetectorPosZ - fSourcePosZ;
    G4double detectorTheta = std::atan(fKermaRadius / distance);
    
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;       // Ω_cône / 4π
    G4double f_det = (1.0 - std::cos(detectorTheta)) / 2.0;     // Ω_det / 4π
    
    // ═══════════════════════════════════════════════════════════════════════
    // FACTEUR DE CORRECTION : f_corr = f_cone
    // ═══════════════════════════════════════════════════════════════════════
    G4double f_corr = f_cone;

    G4cout << "\n╔════════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║  PARAMÈTRES GÉOMÉTRIQUES ET CORRECTION                             ║" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Distance source-détecteur : " << distance/cm << " cm" << G4endl;
    G4cout << "║  Rayon détecteur           : " << fKermaRadius/cm << " cm" << G4endl;
    G4cout << "║  Angle détecteur θ_det     : " << detectorTheta/deg << "°" << G4endl;
    G4cout << "║  Angle cône θ_cône         : " << fConeAngle/deg << "°" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Fraction cône (Ω_cône/4π) : " << f_cone*100 << " %" << G4endl;
    G4cout << "║  Fraction dét. (Ω_det/4π)  : " << f_det*100 << " %" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  *** FACTEUR DE CORRECTION ***                                     ║" << G4endl;
    G4cout << "║  f_corr = Ω_cône / 4π = " << f_corr << G4endl;
    G4cout << "║                                                                    ║" << G4endl;
    G4cout << "║  Note : cos(60°) = 0.5, donc f_cone = (1-0.5)/2 = 0.25            ║" << G4endl;
    G4cout << "║  Le débit brut sera multiplié par 0.25 pour obtenir le réel.      ║" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Fraction gammas vers détecteur : " << f_det/f_cone*100 << " % du cône" << G4endl;
    G4cout << "║  Accélération vs 4π             : " << 1.0/f_cone << "×" << G4endl;
    G4cout << "╚════════════════════════════════════════════════════════════════════╝\n" << G4endl;
}

void RunAction::RecordEventStatistics(G4int nPrimaries,
                                      const std::vector<G4double>& primaryEnergies,
                                      G4int nTransmitted,
                                      G4int nAbsorbed,
                                      G4double doseDeposit)
{
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillH1(0, nPrimaries);

    G4double totalEnergy = 0.;
    for (const auto& energy : primaryEnergies) {
        analysisManager->FillH1(1, energy/keV);
        totalEnergy += energy;
    }

    analysisManager->FillH1(2, totalEnergy/keV);
    analysisManager->FillH1(3, nTransmitted);
    analysisManager->FillH1(4, nAbsorbed);
    analysisManager->FillH1(5, doseDeposit/keV);
    analysisManager->FillH2(0, nPrimaries, totalEnergy/keV);

    fTotalPrimariesGenerated += nPrimaries;
    if (nPrimaries == 0) fTotalEventsWithZeroGamma++;
    fTotalTransmitted += nTransmitted;
    fTotalAbsorbed += nAbsorbed;
}

void RunAction::AddKermaEnergy(G4double edep)
{
    fKermaTotalEnergy += edep;
    fKermaTotalEnergy2 += edep * edep;
    if (edep > 0) fKermaEventCount++;
}

void RunAction::AddKermaFluence(G4double fluence)
{
    fKermaFluenceTotal += fluence;
    fKermaFluenceCount++;
}

void RunAction::AddKermaEnergyForced(G4double forcedDeposit)
{
    fKermaForcedTotal += forcedDeposit;
    fKermaForcedCount++;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) {
        G4cout << "\n### Run ended: No events processed.\n" << G4endl;
        return;
    }

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    // ═══════════════════════════════════════════════════════════════════════
    // STATISTIQUES DE BASE
    // ═══════════════════════════════════════════════════════════════════════
    G4double meanGammasPerEvent = (G4double)fTotalPrimariesGenerated / nofEvents;
    G4double fractionZeroGamma = (G4double)fTotalEventsWithZeroGamma / nofEvents * 100.;

    G4double transmissionRate = 0.;
    G4double absorptionRate = 0.;
    if (fTotalPrimariesGenerated > 0) {
        transmissionRate = (G4double)fTotalTransmitted / fTotalPrimariesGenerated * 100.;
        absorptionRate = (G4double)fTotalAbsorbed / fTotalPrimariesGenerated * 100.;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // CALCULS GÉOMÉTRIQUES
    // ═══════════════════════════════════════════════════════════════════════
    G4double distance = fDetectorPosZ - fSourcePosZ;
    G4double detectorTheta = std::atan(fKermaRadius / distance);
    
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;       // Ω_cône / 4π
    G4double f_det = (1.0 - std::cos(detectorTheta)) / 2.0;     // Ω_det / 4π
    
    // ═══════════════════════════════════════════════════════════════════════
    // *** FACTEUR DE CORRECTION = f_cone ***
    // Pour cône 60° : f_cone = (1 - cos(60°))/2 = 0.25
    // ═══════════════════════════════════════════════════════════════════════
    G4double f_corr = f_cone;

    // ═══════════════════════════════════════════════════════════════════════
    // NORMALISATION TEMPORELLE CORRECTE
    // t = N_events / A_4π (toujours, quelle que soit la géométrie d'émission)
    // ═══════════════════════════════════════════════════════════════════════
    G4double simulatedTime_s = (G4double)nofEvents / fActivity4pi;
    G4double simulatedTime_h = simulatedTime_s / 3600.0;

    // ═══════════════════════════════════════════════════════════════════════
    // MÉTHODE 1 : DÉPÔT D'ÉNERGIE MC (BRUT)
    // ═══════════════════════════════════════════════════════════════════════
    G4double doseTotal_Gy = (fKermaTotalEnergy / keV) * 1.602e-16 / (fKermaMass / kg);
    G4double doseRate_brut_nGyPerH = (doseTotal_Gy / simulatedTime_s) * 3600.0 * 1.0e9;

    // Erreur statistique
    G4double meanDosePerEvent = fKermaTotalEnergy / nofEvents;
    G4double variance = 0.;
    if (nofEvents > 1) {
        variance = (fKermaTotalEnergy2 / nofEvents) - std::pow(meanDosePerEvent, 2);
        variance = std::max(variance, 0.);
    }
    G4double stdDev = std::sqrt(variance);
    G4double stdError = stdDev / std::sqrt(nofEvents);
    G4double convergence = (meanDosePerEvent > 0) ? stdError / meanDosePerEvent * 100. : 0.;

    // ═══════════════════════════════════════════════════════════════════════
    // MÉTHODE 1bis : FORÇAGE D'INTERACTION (BRUT)
    // ═══════════════════════════════════════════════════════════════════════
    G4double doseForced_Gy = (fKermaForcedTotal / fKermaMass) / gray;
    G4double doseRateForced_brut_nGyPerH = (doseForced_Gy / simulatedTime_s) * 3600.0 * 1.0e9;
    G4double convergenceForced = (fKermaForcedCount > 0) ? 
                                  100.0 / std::sqrt((G4double)fKermaForcedCount) : 0.;

    // ═══════════════════════════════════════════════════════════════════════
    // MÉTHODE 2 : FLUENCE (BRUT)
    // ═══════════════════════════════════════════════════════════════════════
    G4double doseFluence_Gy = (fKermaFluenceTotal / fKermaMass) / gray;
    G4double doseRateFluence_brut_nGyPerH = (doseFluence_Gy / simulatedTime_s) * 3600.0 * 1.0e9;
    G4double convergenceFluence = (fKermaFluenceCount > 0) ?
        100.0 / std::sqrt((G4double)fKermaFluenceCount) : 0.;

    // ═══════════════════════════════════════════════════════════════════════
    // *** APPLICATION DE LA CORRECTION GÉOMÉTRIQUE ***
    // Ḋ_réel = Ḋ_brut × f_cone
    // ═══════════════════════════════════════════════════════════════════════
    G4double doseRate_corrige_nGyPerH = doseRate_brut_nGyPerH * f_corr;
    G4double doseRateForced_corrige_nGyPerH = doseRateForced_brut_nGyPerH * f_corr;
    G4double doseRateFluence_corrige_nGyPerH = doseRateFluence_brut_nGyPerH * f_corr;

    G4double doseRateError_nGyPerH = doseRate_corrige_nGyPerH * convergence / 100.;
    G4double doseRateForcedError_nGyPerH = doseRateForced_corrige_nGyPerH * convergenceForced / 100.;
    G4double doseRateFluenceError_nGyPerH = doseRateFluence_corrige_nGyPerH * convergenceFluence / 100.;

    // ═══════════════════════════════════════════════════════════════════════
    // VALEUR THÉORIQUE (calcul détaillé avec μ_en/ρ)
    // ═══════════════════════════════════════════════════════════════════════
    G4double doseRate_theo_nGyPerH = 174.8;  // nGy/h (valeur validée)

    // Estimation du nombre de gammas attendus dans le détecteur
    G4double expectedGammasInDet = fTotalPrimariesGenerated * (f_det / f_cone);

    // ═══════════════════════════════════════════════════════════════════════
    // AFFICHAGE DU RÉSUMÉ
    // ═══════════════════════════════════════════════════════════════════════
    G4cout << "\n";
    G4cout << "╔═══════════════════════════════════════════════════════════════════════════╗\n";
    G4cout << "║                         RUN SUMMARY                                       ║\n";
    G4cout << "║             *** VERSION CÔNE 60° - f_corr = 0.25 ***                      ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  Number of events processed: " << nofEvents << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  PRIMARY GAMMA GENERATION STATISTICS:                                     ║\n";
    G4cout << "║    Total gammas generated     : " << fTotalPrimariesGenerated << G4endl;
    G4cout << "║    Mean gammas per event      : " << meanGammasPerEvent << G4endl;
    G4cout << "║    Expected (theory)          : 1.924" << G4endl;
    G4cout << "║    Events with 0 gamma        : " << fTotalEventsWithZeroGamma
           << " (" << fractionZeroGamma << "%)" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  ÉMISSION ET GÉOMÉTRIE:                                                   ║\n";
    G4cout << "║    Mode émission              : CÔNE de " << fConeAngle/deg << "°" << G4endl;
    G4cout << "║    Distance source-détecteur  : " << distance/cm << " cm" << G4endl;
    G4cout << "║    Rayon détecteur            : " << fKermaRadius/cm << " cm" << G4endl;
    G4cout << "║    Angle détecteur θ_det      : " << detectorTheta/deg << "°" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** NORMALISATION CORRECTE ***                                           ║\n";
    G4cout << "║    Activité (A_4π)            : " << fActivity4pi << " Bq" << G4endl;
    G4cout << "║    Temps simulé               : " << simulatedTime_s << " s ("
           << simulatedTime_h << " h)" << G4endl;
    G4cout << "║    Formule                    : t = N_events / A_4π" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** CORRECTION GÉOMÉTRIQUE ***                                           ║\n";
    G4cout << "║    Fraction cône (Ω_cône/4π)  : " << f_cone*100 << " %" << G4endl;
    G4cout << "║    Fraction dét. (Ω_det/4π)   : " << f_det*100 << " %" << G4endl;
    G4cout << "║    Facteur correction f_corr  : " << f_corr << G4endl;
    G4cout << "║    Formule                    : Ḋ_réel = Ḋ_brut × " << f_corr << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  GAMMAS DANS LE DÉTECTEUR:                                                ║\n";
    G4cout << "║    Gammas entrants (observé)  : " << fKermaForcedCount << G4endl;
    G4cout << "║    Gammas attendus (géom.)    : " << (G4int)expectedGammasInDet << G4endl;
    G4cout << "║    Fraction du cône           : " << 100.0*fKermaForcedCount/fTotalPrimariesGenerated << " %" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 1 : DÉPÔT D'ÉNERGIE MC ***                                   ║\n";
    G4cout << "║    Énergie déposée            : " << fKermaTotalEnergy/keV << " keV" << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRate_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║    Débit CORRIGÉ (× " << f_corr << ")    : " << doseRate_corrige_nGyPerH << " ± " 
           << doseRateError_nGyPerH << " nGy/h" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 1bis : FORÇAGE D'INTERACTION ***                             ║\n";
    G4cout << "║    Gammas avec forçage        : " << fKermaForcedCount << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRateForced_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║    Débit CORRIGÉ (× " << f_corr << ")    : " << doseRateForced_corrige_nGyPerH << " ± " 
           << doseRateForcedError_nGyPerH << " nGy/h" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 2 : FLUENCE ***                                              ║\n";
    G4cout << "║    Gammas traversant          : " << fKermaFluenceCount << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRateFluence_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║    Débit CORRIGÉ (× " << f_corr << ")    : " << doseRateFluence_corrige_nGyPerH << " ± "
           << doseRateFluenceError_nGyPerH << " nGy/h" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** VALEUR THÉORIQUE (Eu-152, calcul détaillé) ***                       ║\n";
    G4cout << "║    Débit théorique            : " << doseRate_theo_nGyPerH << " nGy/h" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** COMPARAISON (Méthode 1bis corrigée vs Théorie) ***                   ║\n";
    G4double ecart = 100.0 * (doseRateForced_corrige_nGyPerH - doseRate_theo_nGyPerH) / doseRate_theo_nGyPerH;
    G4cout << "║    Écart                      : " << ecart << " %" << G4endl;
    if (std::abs(ecart) < 5) {
        G4cout << "║    → EXCELLENT ! Écart < 5%                                             ║\n";
    } else if (std::abs(ecart) < 10) {
        G4cout << "║    → BON ! Écart < 10%                                                  ║\n";
    } else {
        G4cout << "║    → ÉCART SIGNIFICATIF - Investiguer                                   ║\n";
    }
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  OUTPUT FILE: " << fOutputFileName << ".root" << G4endl;
    G4cout << "╚═══════════════════════════════════════════════════════════════════════════╝\n";
    G4cout << G4endl;
}
