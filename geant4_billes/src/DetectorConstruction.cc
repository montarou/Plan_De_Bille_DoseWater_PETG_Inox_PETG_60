#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fDetecteurRP_log(nullptr),
  fPlaqueX(10.0*cm),           // Largeur en X : 10 cm
  fPlaqueY(10.0*cm),           // Largeur en Y : 10 cm
  // ═══════════════════════════════════════════════════════════════════
  // STRUCTURE SANDWICH : W/PETG (7mm) + INOX (4mm) + W/PETG (7mm) = 18 mm
  // ═══════════════════════════════════════════════════════════════════
  fPlaqueZ(18.0*mm),           // Épaisseur totale : 18 mm
  fPlaquePosZ(4.65*cm),        // Centre à z = 4.65 cm (face avant à 3.75 cm)
  fWPETG_thickness(7.0*mm),    // Épaisseur de chaque couche W/PETG : 7 mm
  fInox_thickness(4.0*mm),     // Épaisseur de la couche Inox : 4 mm
  fPETG(nullptr),
  fTungsten(nullptr),
  fW_PETG(nullptr),
  fInox(nullptr)
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();

    // =============================================================================
    // MATÉRIAUX DE BASE
    // =============================================================================

    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

    // =============================================================================
    // TUNGSTÈNE (W) - depuis la base de données NIST
    // =============================================================================
    fTungsten = nist->FindOrBuildMaterial("G4_W");

    // =============================================================================
    // PETG
    // Approximation PETG ~ PET (C10 H8 O4) avec densité typique 1.27 g/cm3
    // =============================================================================
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elO = nist->FindOrBuildElement("O");

    fPETG = new G4Material("PETG", 1.27*g/cm3, 3, kStateSolid);
    fPETG->AddElement(elC, 10);
    fPETG->AddElement(elH,  8);
    fPETG->AddElement(elO,  4);

    // =============================================================================
    // Mélange W/PETG : 75%/25% (fractions massiques)
    // =============================================================================
    G4double rhoW = fTungsten->GetDensity();
    G4double rhoPETG = fPETG->GetDensity();
    G4double massFracW = 0.75;    // 75% en masse de W
    G4double massFracPETG = 0.25; // 25% en masse de PETG

    // Règle des mélanges: 1/ρ_mix = Σ(w_i / ρ_i)
    G4double rhoMix = 1.0 / (massFracW/rhoW + massFracPETG/rhoPETG);

    fW_PETG = new G4Material("W_PETG_75_25", rhoMix, 2, kStateSolid);
    fW_PETG->AddMaterial(fTungsten, massFracW);
    fW_PETG->AddMaterial(fPETG, massFracPETG);

    // Afficher les propriétés du mélange W/PETG
    G4cout << "\n=== Material Properties ===" << G4endl;
    G4cout << "Plaque material: W/PETG 75%/25%" << G4endl;
    G4cout << "  W density      = " << G4BestUnit(rhoW, "Volumic Mass") << G4endl;
    G4cout << "  PETG density   = " << G4BestUnit(rhoPETG, "Volumic Mass") << G4endl;
    G4cout << "  Mix density    = " << G4BestUnit(rhoMix, "Volumic Mass") << G4endl;
    G4cout << "\nDetector material: Water (H2O)" << G4endl;
    G4cout << "  Density = " << water->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
    G4cout << "===========================\n" << G4endl;

    // =============================================================================
    // INOX (Acier inoxydable 304)
    // Composition typique : Fe ~70%, Cr ~18%, Ni ~8%, Mn ~2%, Si ~1%, C ~0.08%
    // Densité : 7.93 - 8.00 g/cm3
    // =============================================================================
    fInox = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    G4cout << "\n=== Material Properties ===" << G4endl;
    G4cout << "Inox (Stainless Steel 304):" << G4endl;
    G4cout << "  Density = " << fInox->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
    G4cout << "===========================\n" << G4endl;

    // =============================================================================
    // WORLD (MONDE)
    // =============================================================================

    G4double world_size = 400*cm;
    G4Box* solidWorld = new G4Box("World",
                                  world_size/2,
                                  world_size/2,
                                  world_size/2);

    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, air, "World");

    G4VPhysicalVolume* physWorld = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     logicWorld,
                                                     "World",
                                                     0,
                                                     false,
                                                     0);

    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ===================================================================
    // ENVELOPPE
    // ===================================================================
    G4Material* envelope_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* solidEnveloppe = new G4Box("Enveloppe", 40*cm, 40*cm, 40*cm);

    G4LogicalVolume* logicEnveloppe = new G4LogicalVolume(solidEnveloppe, envelope_mat, "Enveloppe");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(),
                      logicEnveloppe,
                      "Enveloppe",
                      logicWorld,
                      false,
                      0,
                      true);

    G4VisAttributes* enveloppeVis = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1));
    enveloppeVis->SetVisibility(true);
    logicEnveloppe->SetVisAttributes(enveloppeVis);

    // ═══════════════════════════════════════════════════════════════════
    // STRUCTURE SANDWICH : W/PETG (7mm) + INOX (4mm) + W/PETG (7mm)
    // ═══════════════════════════════════════════════════════════════════
    // Dimensions : 10 cm × 10 cm × 18 mm (total)
    // Face avant : z = 3.75 cm
    // Centre global : z = 4.65 cm
    // Face arrière : z = 5.55 cm
    //
    // Structure (de la source vers le détecteur) :
    //   1. W/PETG avant  : 7 mm (z = 3.75 à 4.45 cm)
    //   2. INOX          : 4 mm (z = 4.45 à 4.85 cm)
    //   3. W/PETG arrière: 7 mm (z = 4.85 à 5.55 cm)
    // ═══════════════════════════════════════════════════════════════════

    G4double plaque_front_z = fPlaquePosZ - fPlaqueZ/2;  // 3.75 cm

    // Positions des centres de chaque couche
    G4double wpetg1_center_z = plaque_front_z + fWPETG_thickness/2;           // 3.75 + 0.35 = 4.10 cm
    G4double inox_center_z   = plaque_front_z + fWPETG_thickness + fInox_thickness/2;  // 3.75 + 0.7 + 0.2 = 4.65 cm
    G4double wpetg2_center_z = plaque_front_z + fWPETG_thickness + fInox_thickness + fWPETG_thickness/2;  // 3.75 + 0.7 + 0.4 + 0.35 = 5.20 cm

    // ───────────────────────────────────────────────────────────────
    // COUCHE 1 : W/PETG AVANT (7 mm)
    // ───────────────────────────────────────────────────────────────
    G4Box* solidWPETG1 = new G4Box("WPETG_Front",
                                    fPlaqueX/2,
                                    fPlaqueY/2,
                                    fWPETG_thickness/2);

    G4LogicalVolume* logicWPETG1 = new G4LogicalVolume(solidWPETG1,
                                                        fW_PETG,
                                                        "WPETG_FrontLog");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, wpetg1_center_z),
                      logicWPETG1,
                      "WPETG_Front",
                      logicEnveloppe,
                      false,
                      0,
                      true);

    // Attributs visuels (vert pour W/PETG)
    G4VisAttributes* wpetgVis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0, 0.5));
    wpetgVis->SetForceSolid(true);
    logicWPETG1->SetVisAttributes(wpetgVis);

    // ───────────────────────────────────────────────────────────────
    // COUCHE 2 : INOX (4 mm)
    // ───────────────────────────────────────────────────────────────
    G4Box* solidInox = new G4Box("Inox",
                                  fPlaqueX/2,
                                  fPlaqueY/2,
                                  fInox_thickness/2);

    G4LogicalVolume* logicInox = new G4LogicalVolume(solidInox,
                                                      fInox,
                                                      "InoxLog");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, inox_center_z),
                      logicInox,
                      "Inox",
                      logicEnveloppe,
                      false,
                      0,
                      true);

    // Attributs visuels (gris métallique pour Inox)
    G4VisAttributes* inoxVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.7));
    inoxVis->SetForceSolid(true);
    logicInox->SetVisAttributes(inoxVis);

    // ───────────────────────────────────────────────────────────────
    // COUCHE 3 : W/PETG ARRIÈRE (7 mm)
    // ───────────────────────────────────────────────────────────────
    G4Box* solidWPETG2 = new G4Box("WPETG_Back",
                                    fPlaqueX/2,
                                    fPlaqueY/2,
                                    fWPETG_thickness/2);

    G4LogicalVolume* logicWPETG2 = new G4LogicalVolume(solidWPETG2,
                                                        fW_PETG,
                                                        "WPETG_BackLog");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, wpetg2_center_z),
                      logicWPETG2,
                      "WPETG_Back",
                      logicEnveloppe,
                      false,
                      0,
                      true);

    logicWPETG2->SetVisAttributes(wpetgVis);

    // Face arrière du sandwich
    G4double plaque_back_z = fPlaquePosZ + fPlaqueZ/2;   // 5.55 cm

    // =============================================================================
    // PLANS DE COMPTAGE (UPSTREAM ET DOWNSTREAM)
    // =============================================================================

    G4double detector_thickness = 2.0*mm;
    G4double detector_gap = 2.0*mm;
    G4double detector_size = 12.0*cm;

    G4double upstream_detector_z = plaque_front_z - detector_gap - detector_thickness/2.0;
    G4double downstream_detector_z = plaque_back_z + detector_gap + detector_thickness/2.0;

    G4Box* solidDetector = new G4Box("Detector",
                                     detector_size/2.0,
                                     detector_size/2.0,
                                     detector_thickness/2.0);

    // ───────────────────────────────────────────────────────────────
    // UPSTREAM DETECTOR
    // ───────────────────────────────────────────────────────────────

    G4LogicalVolume* logicUpstreamDetector = new G4LogicalVolume(solidDetector, air, "UpstreamDetector");

    G4VisAttributes* upstreamVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    upstreamVisAtt->SetForceSolid(true);
    logicUpstreamDetector->SetVisAttributes(upstreamVisAtt);

    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., upstream_detector_z),
                      logicUpstreamDetector,
                      "UpstreamDetector",
                      logicEnveloppe,
                      false,
                      0,
                      false);

    // ───────────────────────────────────────────────────────────────
    // DOWNSTREAM DETECTOR
    // ───────────────────────────────────────────────────────────────

    G4LogicalVolume* logicDownstreamDetector = new G4LogicalVolume(solidDetector, air, "DownstreamDetector");

    G4VisAttributes* downstreamVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.3));
    downstreamVisAtt->SetForceSolid(true);
    logicDownstreamDetector->SetVisAttributes(downstreamVisAtt);

    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., downstream_detector_z),
                      logicDownstreamDetector,
                      "DownstreamDetector",
                      logicEnveloppe,
                      false,
                      0,
                      false);

    // =============================================================================
    // DÉTECTEUR DOSE (sphère d'EAU à z = 20 cm)
    // =============================================================================

    G4double distance_plaque_detecteur = 14.45*cm;
    G4double dose_position_z = plaque_back_z + distance_plaque_detecteur;
    G4double dose_radius = 2.0*cm;

    G4Orb* doseDetector_solid = new G4Orb("DoseDetector", dose_radius);

    G4LogicalVolume* doseDetector_log = new G4LogicalVolume(
        doseDetector_solid,
        water,
        "DoseDetectorLog");

    G4double maxStep = 0.1*mm;
    G4UserLimits* doseLimits = new G4UserLimits(maxStep);
    doseDetector_log->SetUserLimits(doseLimits);

    G4VisAttributes* doseVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.5));
    doseVis->SetForceSolid(true);
    doseDetector_log->SetVisAttributes(doseVis);

    G4ThreeVector dose_position = G4ThreeVector(0, 0, dose_position_z);

    new G4PVPlacement(
        0,
        dose_position,
        doseDetector_log,
        "DoseDetector",
        logicEnveloppe,
        false,
        0,
        true);

    // =============================================================================
    // CALCULS POUR L'AFFICHAGE
    // =============================================================================
    G4double dose_volume = (4.0/3.0) * M_PI * std::pow(dose_radius, 3);
    G4double water_density = water->GetDensity();
    G4double dose_mass = dose_volume * water_density;

    G4double wpetg_volume = fPlaqueX * fPlaqueY * fWPETG_thickness;
    G4double inox_volume = fPlaqueX * fPlaqueY * fInox_thickness;
    G4double wpetg_mass = wpetg_volume * fW_PETG->GetDensity();
    G4double inox_mass = inox_volume * fInox->GetDensity();
    G4double total_mass = 2 * wpetg_mass + inox_mass;

    G4double source_z = 2.0*cm;
    G4double distance_source_detecteur = dose_position_z - source_z;

    // =============================================================================
    // AFFICHAGE RÉCAPITULATIF
    // =============================================================================

    G4cout << "\n╔══════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║      PROGRAMME : MESURE DE DOSE DANS L'EAU                        ║" << G4endl;
    G4cout << "║      (Version SANDWICH : W/PETG + INOX + W/PETG)                  ║" << G4endl;
    G4cout << "╠══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "║ STRUCTURE SANDWICH :                                              ║" << G4endl;
    G4cout << "║ ┌────────────────────────────────────────────────────────────┐   ║" << G4endl;
    G4cout << "║ │  W/PETG avant  │     INOX      │  W/PETG arrière │         │   ║" << G4endl;
    G4cout << "║ │     7 mm       │     4 mm      │     7 mm        │ = 18 mm │   ║" << G4endl;
    G4cout << "║ └────────────────────────────────────────────────────────────┘   ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "║ • Dimensions X×Y     : " << fPlaqueX/cm << " × " << fPlaqueY/cm << " cm                          ║" << G4endl;
    G4cout << "║ • Épaisseur totale   : " << fPlaqueZ/mm << " mm                                ║" << G4endl;
    G4cout << "║ • Position centre    : z = " << fPlaquePosZ/cm << " cm                         ║" << G4endl;
    G4cout << "║ • Face avant         : z = " << plaque_front_z/cm << " cm                        ║" << G4endl;
    G4cout << "║ • Face arrière       : z = " << plaque_back_z/cm << " cm                        ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ COUCHE 1 - W/PETG AVANT :                                         ║" << G4endl;
    G4cout << "║ • Épaisseur          : " << fWPETG_thickness/mm << " mm                                 ║" << G4endl;
    G4cout << "║ • Position centre    : z = " << wpetg1_center_z/cm << " cm                         ║" << G4endl;
    G4cout << "║ • Densité            : " << fW_PETG->GetDensity()/(g/cm3) << " g/cm³                       ║" << G4endl;
    G4cout << "║ • Masse              : " << wpetg_mass/g << " g                             ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ COUCHE 2 - INOX :                                                 ║" << G4endl;
    G4cout << "║ • Épaisseur          : " << fInox_thickness/mm << " mm                                 ║" << G4endl;
    G4cout << "║ • Position centre    : z = " << inox_center_z/cm << " cm                        ║" << G4endl;
    G4cout << "║ • Densité            : " << fInox->GetDensity()/(g/cm3) << " g/cm³                            ║" << G4endl;
    G4cout << "║ • Masse              : " << inox_mass/g << " g                               ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ COUCHE 3 - W/PETG ARRIÈRE :                                       ║" << G4endl;
    G4cout << "║ • Épaisseur          : " << fWPETG_thickness/mm << " mm                                 ║" << G4endl;
    G4cout << "║ • Position centre    : z = " << wpetg2_center_z/cm << " cm                          ║" << G4endl;
    G4cout << "║ • Masse              : " << wpetg_mass/g << " g                             ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ MASSE TOTALE SANDWICH : " << total_mass/g << " g                            ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ PLANS DE COMPTAGE :                                               ║" << G4endl;
    G4cout << "║ • Upstream           : z = " << upstream_detector_z/cm << " cm                       ║" << G4endl;
    G4cout << "║ • Downstream         : z = " << downstream_detector_z/cm << " cm                       ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╟───────────────────────────────────────────────────────────────────╢" << G4endl;
    G4cout << "║ DÉTECTEUR DOSE DANS L'EAU :                                       ║" << G4endl;
    G4cout << "║ • Position centre    : z = " << dose_position_z/cm << " cm                           ║" << G4endl;
    G4cout << "║ • Distance source    : " << distance_source_detecteur/cm << " cm                           ║" << G4endl;
    G4cout << "║ • Rayon              : " << dose_radius/cm << " cm                                 ║" << G4endl;
    G4cout << "║ • Masse              : " << dose_mass/g << " g                            ║" << G4endl;
    G4cout << "║                                                                   ║" << G4endl;
    G4cout << "╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;

    return physWorld;
}
