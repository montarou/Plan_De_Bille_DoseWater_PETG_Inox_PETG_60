#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();

private:
    G4LogicalVolume* fDetecteurRP_log;
    
    // Paramètres de la plaque (accessibles si besoin)
    G4double fPlaqueX;      // Largeur en X
    G4double fPlaqueY;      // Largeur en Y  
    G4double fPlaqueZ;      // Épaisseur totale en Z (18 mm)
    G4double fPlaquePosZ;   // Position du centre en Z

    // ═══════════════════════════════════════════════════════════════════
    // STRUCTURE SANDWICH : W/PETG (7mm) + INOX (4mm) + W/PETG (7mm)
    // ═══════════════════════════════════════════════════════════════════
    G4double fWPETG_thickness;  // Épaisseur de chaque couche W/PETG (7 mm)
    G4double fInox_thickness;   // Épaisseur de la couche Inox (4 mm)

    // Matériaux
    G4Material* fPETG;          // Matériau PETG
    G4Material* fTungsten;      // Matériau Tungstène (W)
    G4Material* fW_PETG;        // Mélange W/PETG 75%/25%
    G4Material* fInox;          // Matériau Inox (acier inoxydable 304)
};

#endif
