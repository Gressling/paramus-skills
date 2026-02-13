# SKILL.md — PARAMUS Chemistry & Materials Science

## Server Identity

- **Name**: `io.github.gressling/paramus-chemistry`
- **Title**: PARAMUS Chemistry & Materials Science
- **Version**: 2.0.0
- **Endpoint**: `https://cloud1.paramus.ai/mcp`

---

## Skills

### chemistry.molecular_properties

**Description**: Calculate molecular properties from SMILES strings  
**Category**: Chemistry  
**Tools**: `calculate_molecular_weight`, `calculate_logp`, `calculate_tpsa`, `calculate_num_hba`, `calculate_num_hbd`, `calculate_num_rotatable_bonds`, `calculate_num_rings`, `calculate_num_aromatic_rings`, `calculate_fraction_csp3`, `calculate_formal_charge`, `calculate_molar_refractivity`, `calculate_exact_mass`, `calculate_heavy_atom_count`, `calculate_heavy_atom_molecular_weight`

**Example prompts**:
- "What is the molecular weight of aspirin?"
- "Calculate LogP for caffeine (SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
- "Get the topological polar surface area of ibuprofen"
- "How many hydrogen bond donors does paracetamol have?"

**Input**: SMILES string (e.g., `CC(=O)Oc1ccccc1C(=O)O`)  
**Output**: Numeric property value with units

---

### chemistry.molecular_conversion

**Description**: Convert between molecular representations  
**Category**: Chemistry  
**Tools**: `smiles_to_inchi`, `smiles_to_inchikey`, `inchi_to_smiles`, `smiles_to_canonical`, `smiles_to_mol_block`, `validate_smiles`

**Example prompts**:
- "Convert aspirin SMILES to InChI"
- "What is the InChIKey for CC(=O)O?"
- "Is this a valid SMILES string: C(C)O?"
- "Canonicalise this SMILES: c1ccccc1O"

**Input**: SMILES or InChI string  
**Output**: Converted representation

---

### chemistry.fingerprints_similarity

**Description**: Generate molecular fingerprints and compute similarity  
**Category**: Chemistry  
**Tools**: `generate_morgan_fingerprint`, `generate_rdkit_fingerprint`, `generate_maccs_keys`, `calculate_tanimoto_similarity`, `calculate_dice_similarity`, `calculate_fingerprint_similarity_matrix`

**Example prompts**:
- "How similar are aspirin and ibuprofen?"
- "Generate a Morgan fingerprint for ethanol"
- "Calculate Tanimoto similarity between caffeine and theobromine"

**Input**: One or two SMILES strings + fingerprint parameters  
**Output**: Fingerprint bit vector or similarity coefficient (0–1)

---

### chemistry.structure_analysis

**Description**: Analyse molecular structure and substructures  
**Category**: Chemistry  
**Tools**: `check_aromaticity`, `find_substructure`, `get_ring_info`, `get_atom_count`, `enumerate_stereoisomers`, `generate_3d_conformer`, `calculate_3d_descriptors`

**Example prompts**:
- "Is this molecule aromatic?"
- "Does aspirin contain a carboxylic acid group?"
- "How many rings does caffeine have?"
- "Generate a 3D conformer for benzene"

**Input**: SMILES string + optional substructure SMARTS  
**Output**: Structural analysis results

---

### chemistry.descriptors

**Description**: Calculate comprehensive molecular descriptors  
**Category**: Chemistry  
**Tools**: `calculate_all_descriptors`, `calculate_drug_likeness`, `check_lipinski_rule`, `check_veber_rules`, `calculate_qed`, `calculate_synthetic_accessibility`

**Example prompts**:
- "Does this molecule satisfy Lipinski's Rule of Five?"
- "What is the QED score for this drug candidate?"
- "Calculate synthetic accessibility for CC(=O)Oc1ccccc1C(=O)O"
- "Run a full descriptor calculation for paracetamol"

**Input**: SMILES string  
**Output**: Descriptor values, pass/fail assessments

---

### thermodynamics.fluid_properties

**Description**: Query thermodynamic properties of pure fluids using CoolProp  
**Category**: Thermodynamics  
**Tools**: `coolprop_fluid_properties`, `coolprop_saturation_properties`, `coolprop_phase_envelope`, `coolprop_transport_properties`, `coolprop_critical_point`, `coolprop_triple_point`, `coolprop_available_fluids`

**Example prompts**:
- "What is the density of water at 25°C and 1 atm?"
- "Get the viscosity of ethanol at 78°C"
- "What are the critical point properties of CO2?"
- "List all available fluids in CoolProp"
- "Calculate Cp and Cv for nitrogen at 300K"

**Input**: Fluid name + temperature (K) + pressure (Pa)  
**Output**: Thermodynamic properties (density, viscosity, Cp, Cv, enthalpy, entropy, etc.)

---

### thermodynamics.kinetics

**Description**: Chemical kinetics calculations using Cantera  
**Category**: Thermodynamics  
**Tools**: `cantera_equilibrium`, `cantera_reaction_rates`, `cantera_flame_speed`, `cantera_ignition_delay`

**Example prompts**:
- "What is the adiabatic flame temperature of methane in air?"
- "Calculate the equilibrium composition of H2/O2 at 2000K"

**Input**: Mechanism, species, conditions  
**Output**: Equilibrium state, reaction rates, or flame properties

---

### electrochemistry.fundamentals

**Description**: Electrochemistry calculations  
**Category**: Electrochemistry  
**Tools**: `nernst_equation`, `butler_volmer`, `calculate_conductivity`, `electrode_potential`, `faraday_electrolysis`

**Example prompts**:
- "Calculate the Nernst potential for a Zn/Cu cell"
- "What current density from Butler-Volmer at 50mV overpotential?"
- "How much copper is deposited by 10A for 1 hour?"

**Input**: Electrochemical parameters (potential, concentration, temperature, etc.)  
**Output**: Calculated electrochemical quantities

---

### materials.polymer_prediction

**Description**: Predict polymer and material properties using ML models  
**Category**: Materials Science  
**Tools**: `predict_tg`, `polymer_fingerprint`, `bigsmiles_parse`, `bigsmiles_validate`

**Example prompts**:
- "Predict the glass transition temperature of polystyrene"
- "Parse this BigSMILES string"
- "Generate a polymer fingerprint for polyethylene"

**Input**: SMILES, BigSMILES, or text description  
**Output**: Predicted properties or processed polymer representations

---

### data_science.statistics

**Description**: Statistical analysis and experimental design  
**Category**: Chemistry (cross-cutting)  
**Tools**: `doe_full_factorial`, `doe_latin_hypercube`, `pca_analysis`, `kmeans_clustering`, `regression_analysis`, `descriptive_statistics`

**Example prompts**:
- "Design a full factorial experiment with 3 factors at 2 levels"
- "Run PCA on this dataset"
- "Cluster these molecules by their descriptors"

**Input**: Data arrays or design parameters  
**Output**: Statistical results, designs, or cluster assignments

---

## Capabilities Summary

| Capability | Supported |
|---|---|
| **Tools** | 319 tools via 7 gateway meta-tools |
| **Resources** | Not provided |
| **Prompts** | Not provided |
| **Sampling** | Not used |
| **Authentication** | OAuth 2.0 + PKCE |
| **Transport** | Streamable HTTP, SSE, JSON-RPC |

## Intended Users

- **Researchers** in chemistry, materials science, and chemical engineering
- **Students** learning computational chemistry
- **Engineers** designing materials, evaluating drug candidates, or querying thermodynamic data
- **Data scientists** working with molecular datasets

## Limitations

- Tools operate on individual molecules or small datasets — not batch processing
  of millions of compounds.
- ML predictions (e.g., Tg) are model-dependent and should be validated
  experimentally for critical applications.
- CoolProp covers ~120 pure fluids — not every substance is available.
- No arbitrary code execution — all tools are pre-registered handler functions.
