# Microkinetic model for ammonia oxidation
# E.V. Rebrov, M.H.J.M. de Croon, J.C. Schouten
# Development of the kinetic model of platinum catalyzed ammonia oxidation in a microreactor
# Chemical Engineering Journal 90 (2002) 61–76

database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','CHON_G4','NOx2018', 'GRI-Mech3.0-N', 'NitrogenCurran','DFT_QCI_thermo'],
    reactionLibraries =['Surface/CPOX_Pt/Deutschmann2006'],
    #reactionLibraries =['Surface/CPOX_Pt/Deutschmann2006','Surface/Rebrov_Pt111','Surface/Nitrogen','Surface/Arevalo_Pt111','Surface/Kraehnert_Pt111','Surface/Schneider_Pt111','Surface/Novell_Pt111','Surface/Offermans_Pt111','Surface/Scheuer_Pt'],
    seedMechanisms = ['Surface/Schneider_Rh211','Surface/Nitrogen'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface'],
    #kineticsFamilies = [
    #'Surface_Adsorption_Single',
    #'Surface_Adsorption_vdW',
    #'Surface_Adsorption_Dissociative',
    #'Surface_Dissociation',
    #'Surface_Abstraction',
    #'Surface_EleyRideal_Addition_Multiple_Bond',
    #'Surface_Migration',
    #'Surface_Dissociation_Double_vdW',
    #'Surface_Addition_Single_vdW',
    #'Surface_Dissociation_vdW',
    #'Surface_Abstraction_vdW',
    #'Surface_Dual_Adsorption_vdW',
    #'Surface_Dissociation_Beta',
    #'Surface_Adsorption_Abstraction_vdW',
    #'Surface_Adsorption_Bidentate',
    #'Surface_Bidentate_Dissociation',
    #'Surface_DoubleBond_to_Bidentate', 
    #'Surface_vdW_to_Bidentate',
    #'Surface_Abstraction_Single_vdW',
    #'Surface_Adsorption_Dissociative_Double',
    #'default'
    #],
    kineticsEstimator = 'rate rules',
)

catalystProperties(
    metal = 'Rh211'
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=0,
    maximumOxygenAtoms=2,
    maximumNitrogenAtoms=2,
    maximumSurfaceSites=2,
    maximumRadicalElectrons=2,
)

# List of species
species(
    label='X',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
    label='O2',
    reactive=True,
    structure=adjacencyList(
"""
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)

species(
    label='N2',
    reactive=True,
    structure=SMILES("N#N"),
)

species(
    label='NO',
    reactive=True,
    structure=adjacencyList(
"""
multiplicity 2
1 N u1 p1 c0 {2,D}
2 O u0 p2 c0 {1,D}
"""),
)

species(
    label='NH3',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
)

species(
    label='N2O',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p2 c-1 {2,D}
2 N u0 p0 c+1 {1,D} {3,D}
3 O u0 p2 c0 {2,D}
"""),
)

species(
    label='NO2',
    reactive=True,
    structure=adjacencyList(
"""
multiplicity 2
1 N u0 p1 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 O u1 p2 c0 {1,S}
"""),
)

species(
    label='He',
    reactive=False,
    structure=adjacencyList(
"""
1 He u0 p1 c0
"""),
)

species(
    label='H2OX',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0
"""),
)
species(
    label='NOX',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {2,D}
2 N u0 p1 c0 {1,D} {3,S}
3 X u0 p0 c0 {2,S}
"""),
)

species(
    label='OX',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
"""),
)

species(
    label='NX',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
"""),
)

species(
    label='OHX',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
"""),
)

species(
    label='NH3X',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0
"""),
)

species(
    label='NH2X',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0 {1,S}
"""),
)

species(
    label='NHX',
    reactive=True,
    structure=adjacencyList(
"""
1 N u0 p1 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
"""),
)

species(
    label='N2OX',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {2,D}
2 N u0 p1 c0 {1,D} {3,S}
3 N u0 p1 c0 {2,S} {4,D}
4 X u0 p0 c0 {3,D}
"""),
)

species(
    label='NO2X',
    reactive=True,
    structure=adjacencyList(
"""
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 N u0 p1 c0 {1,S} {2,D}
4 X u0 p0 c0 {1,S}
"""),
)


#-------------

#temperature from 523-673K 
surfaceReactor(  
    temperature=[(300,'K'),(1500,'K')],
    initialPressure=(1.0, 'bar'),
    nSims=6,
    initialGasMoleFractions={
        "NH3": 0.002,
        "O2": 0.01,
        "He": 0.938,
        "NO":0.0,
        "H2O":0.05,
        "N2O":0.0,
        "N2":0.0,
        "NO2":0.0,
        "NX":0.0,
        "NOX":0.0,
        "NH3X":0.0,
        "NH2X":0.0,
        "NHX":0.0,
        "OX":0.0,
        "OHX":0.0,
        "H2OX":0.0,
        "N2OX":0.0,
        "NO2X":0.0,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(2.3832928e4, 'm^-1'), #Spherical area=0.0660382956m2,V of Pt=0.00002531468cm3 
    #surfaceVolumeRatio=(1.4285714e4, 'm^-1'), #A/V = 280µm*π*9mm/140µm*140µm*π*9mm = 2.8571428e4^m-1
    terminationConversion = {"NH3":0.95,},
    #terminationTime=(10, 's'),
)

simulator( #default for surface reaction atol=1e-18,rtol=1e-12
    atol=1e-22, #absolute tolerance are 1e-15 to 1e-25
    rtol=1e-12, #relative tolerance is usually 1e-4 to 1e-8
)

model( 
    toleranceKeepInEdge=0.05,
    toleranceMoveToCore=0.5, 
    toleranceInterruptSimulation=1e6, 
    maximumEdgeSpecies=5000, 
    minCoreSizeForPrune=50,
    toleranceThermoKeepSpeciesInEdge=0.5, # prune before simulation based on thermo
    minSpeciesExistIterationsForPrune=2, # prune rxns from edge that dont move into core
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=True, 
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
