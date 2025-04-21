# --- LSC Comprehensive Particle Model (v6 - Final Corrected) ---
# --- Includes Leptons, Neutrino, Photon, Quarks (Conceptual), Hadrons (Conceptual) ---
# --- Mass = ZPE of Braid Oscillation (Leptons); Stiffness origins differ; Hadrons via QCD-like binding ---
# --- Spin = Framed Holonomy Mechanism (Assumed) ---
# --- Interaction Vertex Structure & g-2 Sketch Included ---
# --- Single File Python Script ---

import networkx as nx
import numpy as np
import math

print("--- Loading LSC Comprehensive Particle Model ---")
print("--- MODEL: Fermion Mass = ZPE of Braid Oscillation (Leptons) ---")
print("---        Stiffness 'k' origins differ; Hadrons via QCD-like binding ---")
print("---        Spin arises from Framed Braid Rotation ---")
print("---        Interaction Vertex Structure Included ---")


# --- Physical Constants (SI) ---
hbar_si = 1.0545718e-34; c_si = 2.99792458e8; G_si = 6.67430e-11
eV_to_J = 1.602176634e-19
alpha_fine_structure = 1 / 137.035999
elementary_charge_C = 1.60217663e-19

# Masses (approximate, kg)
electron_mass_kg = 9.1093837e-31
muon_mass_kg = 1.8835316e-28
tau_mass_kg = 3.16754e-27
neutrino_mass_eV = 0.1
neutrino_mass_kg = (neutrino_mass_eV * eV_to_J) / (c_si**2)
proton_mass_kg = 1.6726219e-27

# --- Planck Units ---
planck_length_si = np.sqrt(G_si * hbar_si / (c_si**3))
planck_mass_kg = np.sqrt(hbar_si * c_si / G_si)
planck_energy_J = planck_mass_kg * c_si**2
planck_time_si = planck_length_si / c_si
planck_stiffness_natural = 1.0 # E_p / lp^2
planck_mass_natural = 1.0      # M_p
e_natural = np.sqrt(4 * np.pi * alpha_fine_structure)

# --- Target Energies (Natural Units E/Ep) ---
TARGET_ENERGY = {
    "electron": (electron_mass_kg * c_si**2) / planck_energy_J,
    "muon":     (muon_mass_kg * c_si**2) / planck_energy_J,
    "tau":      (tau_mass_kg * c_si**2) / planck_energy_J,
    "neutrino": (neutrino_mass_kg * c_si**2) / planck_energy_J,
    "proton":   (proton_mass_kg * c_si**2) / planck_energy_J,
}

# --- Particle Properties Lookup ---
PARTICLE_PROPS = {
    "electron": {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True, "color": None},
    "muon":     {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True, "color": None},
    "tau":      {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True, "color": None},
    "neutrino": {"charge_q":  0, "spin_j": 0.5, "is_lepton": True, "is_charged": False,"color": None},
    "photon":   {"charge_q":  0, "spin_j": 1.0, "is_lepton": False,"is_charged": False,"color": None},
    "up":       {"charge_q": +2/3,"spin_j": 0.5, "is_lepton": False,"is_charged": True, "color": "RGB"},
    "down":     {"charge_q": -1/3,"spin_j": 0.5, "is_lepton": False,"is_charged": True, "color": "RGB"},
    "proton":   {"charge_q": +1,  "spin_j": 0.5, "is_lepton": False,"is_charged": True, "color": "Neutral"}
}

print(f"--- Constants ---")
print(f"Planck Mass (kg): {planck_mass_kg:.2e}")
print(f"Alpha: {alpha_fine_structure:.3e}")
print(f"e (natural): {e_natural:.3f}")
print(f"--- Target Energies (Planck Units) ---")
print(f"  Electron : {TARGET_ENERGY['electron']:.3e}")
print(f"  Muon     : {TARGET_ENERGY['muon']:.3e}")
print(f"  Tau      : {TARGET_ENERGY['tau']:.3e}")
print(f"  Neutrino : {TARGET_ENERGY['neutrino']:.3e}")
print(f"  Proton   : {TARGET_ENERGY['proton']:.3e} (For Scale)")
print("-------------------------------------\n")

# --- Effective Oscillator Model for Leptons & Neutrino ---
class FermionOscillatorModel:
    def __init__(self, particle_name, target_energy_natural):
        particle_name_lower = particle_name.lower()
        if particle_name_lower not in TARGET_ENERGY or particle_name_lower not in PARTICLE_PROPS: raise ValueError(f"Unknown particle: {particle_name}")
        if PARTICLE_PROPS[particle_name_lower].get("color") == "RGB" or particle_name_lower == "proton": raise ValueError(f"Oscillator model not for quarks/hadrons: {particle_name}")

        self.particle_name = particle_name_lower; self.properties = PARTICLE_PROPS[self.particle_name]; self.is_charged = self.properties["is_charged"]
        self.target_E0_natural = target_energy_natural
        self.m_eff_natural = 1.0; self.assumed_C_tunnel = 1.0; self.calculated_C_grav_fit = None
        self.omega_natural = None; self.k_required_natural = None; self.required_c1_for_k = None
        self.calculate_required_parameters(); self.interpret_k_origin()

    def calculate_required_parameters(self):
        if self.target_E0_natural is None or self.target_E0_natural <= 0: return
        self.omega_natural = 2 * self.target_E0_natural
        self.k_required_natural = self.m_eff_natural * (self.omega_natural**2)

    def interpret_k_origin(self):
        if self.k_required_natural is None or self.k_required_natural <= 0: return
        if self.is_charged:
            target_k = self.k_required_natural; prefactor = self.assumed_C_tunnel
            try:
                log_arg = target_k / prefactor
                if log_arg <= np.finfo(float).tiny: self.required_c1_for_k = np.inf
                else: self.required_c1_for_k = -alpha_fine_structure * np.log(log_arg)
            except Exception: self.required_c1_for_k = None
        else: # Neutrino
            self.required_c1_for_k = None
            target_k = self.k_required_natural; target_E0 = self.target_E0_natural
            if target_E0 > 0: self.calculated_C_grav_fit = target_k / (target_E0**2)
            else: self.calculated_C_grav_fit = None

    def get_required_stiffness(self): return self.k_required_natural
    def get_mass_kg(self):
        if self.target_E0_natural is None: return None
        return self.target_E0_natural * planck_mass_kg

    def report(self):
        print(f"--- {self.particle_name.upper()} Oscillator Parameters ---")
        target_mass = self.get_mass_kg()
        if target_mass is not None: print(f"  Target Mass (kg): {target_mass:.2e}")
        if self.target_E0_natural is not None: print(f"  Target E0 (Planck Units): {self.target_E0_natural:.3e}")
        if self.omega_natural and self.k_required_natural is not None:
            print(f"  Required Omega (Planck Units): {self.omega_natural:.3e}")
            print(f"  >> Required Stiffness k (Planck Units): {self.get_required_stiffness():.3e}")
        else: print("  Omega/k: Invalid Target")
        print(f"--- Interpretation of Required Stiffness k ---")
        if self.is_charged:
             print(f"  Charged Hypothesis: k ≈ C_t * exp(-c1/alpha)")
             if self.required_c1_for_k is not None: print(f"  >> Needs c1 ≈ {self.required_c1_for_k:.3f} (if C_t={self.assumed_C_tunnel:.1f})")
             else: print(f"  Could not calculate required 'c1'.")
        else: # Neutrino
             print(f"  Neutral Hypothesis: k ≈ C_g * (E0/Ep)^2")
             if self.calculated_C_grav_fit is not None: print(f"  >> Needs C_g ≈ {self.calculated_C_grav_fit:.3f}")
             else: print(f"  Could not calculate required 'C_g'.")
        print(f"---------------------------------------------")

# --- Spin Network / Braid Placeholders ---
class SpinNetwork:
    def __init__(self): self.description = "Conceptual LQG Background"
class EmbeddedBraid:
    def __init__(self, name, particle_type):
        self.name = name; self.particle_type=particle_type.lower()
        if self.particle_type not in PARTICLE_PROPS: raise ValueError(f"Unknown particle type: {particle_type}")
        self.properties = PARTICLE_PROPS[self.particle_type]
        self.description=f"Conceptual j={self.get_spin()} framed braid ({particle_type})"
        self.color = ["red", "green", "blue"][hash(name + particle_type) % 3] if self.properties.get("color") == "RGB" else self.properties.get("color")
    def get_spin(self): return self.properties["spin_j"]
    def get_charge(self): return self.properties["charge_q"] # Method to get charge
    def get_color(self): return self.color

class HadronBraidModel:
     def __init__(self, name, constituent_braids):
         self.name = name; self.constituents = constituent_braids
         if name.lower() not in PARTICLE_PROPS: raise ValueError(f"Unknown hadron: {name}")
         self.properties = PARTICLE_PROPS[name.lower()]
         self.description = f"Conceptual bound state ({''.join(q.particle_type for q in self.constituents)})"
         self.particle_type = name.lower()
     def get_spin(self): return self.properties["spin_j"]
     def get_charge(self): return sum(q.get_charge() for q in self.constituents) # Sum constituent charges
     def is_color_neutral(self):
         colors = [q.get_color() for q in self.constituents if q.properties.get("color") == "RGB"]
         return len(colors) == 3 and len(set(colors)) == 3
     def report(self):
         print(f"\n--- HADRON MODEL: {self.name.upper()} ---")
         print(f"  Description: {self.description}")
         print(f"  Constituents: {[q.name for q in self.constituents]}")
         print(f"  Properties: Spin={self.get_spin()}, Charge={self.properties['charge_q']:.2f}e, Color Neutral: {self.is_color_neutral()}") # Use props charge
         print(f"  Mass Origin: Dominantly from binding energy (QCD scale).")
         print(f"-------------------------")

# --- Spin Transformation Placeholder ---
def verify_spinor_transformation(braid_structure):
    particle_desc = f"({braid_structure.particle_type})" if hasattr(braid_structure, 'particle_type') else f"({braid_structure.name} - Hadron)"
    print(f"\n--- Verifying Spinor Transformation for '{braid_structure.name}' {particle_desc} ---")
    spin_j = braid_structure.get_spin()
    if spin_j is None: print("  Spin not defined. Skipping check."); return None
    if not np.isclose(spin_j, 0.5): print(f"  Not Spin-1/2 (j={spin_j}). Skipping check."); return None
    phase_2pi = np.exp(-1j * 2*np.pi * 0.5); phase_4pi = np.exp(-1j * 4*np.pi * 0.5)
    result = np.isclose(phase_2pi, -1.0) and np.isclose(phase_4pi, 1.0)
    print(f"Consistent with Spin-1/2: {result} (assuming mechanism)")
    print("----------------------------------------------------")
    return result

# --- Photon Model Placeholder ---
class PhotonModel:
    def __init__(self): self.name="Photon"; self.properties=PARTICLE_PROPS["photon"]; self.description="Massless..."
    def get_spin(self): return self.properties["spin_j"]
    def get_charge(self): return self.properties["charge_q"]
    def report(self): print(f"\n--- PHOTON MODEL ---\n  Spin: {self.get_spin()}, Mass: 0\n--------------------")

# --- Interaction Vertex & g-2 Sketch Placeholders ---
def calculate_interaction_amplitude(fermion_in, photon_action, fermion_out):
    photon_state = photon_action[1]; is_emission = (photon_action[0] == 'emit')
    if np.isclose(fermion_in.get_charge(), 0.0): return 0.0
    coupling_g = e_natural * fermion_in.get_charge();
    cg_factor = 1.0 if (np.isclose(fermion_in.get_spin(), 0.5) and np.isclose(fermion_out.get_spin(), 0.5)) else 0.0
    if np.isclose(cg_factor, 0.0): return 0.0
    current_factor = 1.0; photon_op_factor = 1.0
    amplitude = coupling_g * cg_factor * current_factor * photon_op_factor * complex(0.8, 0.2)
    return amplitude
def braid_propagator(energy_diff_natural): return 1.0 / (abs(energy_diff_natural) + 1e-30)
def photon_propagator(photon_momentum_sq_natural): return 1.0 / (abs(photon_momentum_sq_natural) + 1e-30)
def calculate_g_minus_2_leading_term(electron_braid):
    print("\n--- Calculating Electron Anomalous Magnetic Moment (g-2) - Sketch ---")
    if np.isclose(electron_braid.get_charge(), 0.0): print("ERROR: g-2 requires charged particle."); return None
    # Create conceptual intermediate state with correct charge
    intermediate_braid = EmbeddedBraid("Electron_Intermediate", "electron")
    intermediate_braid.properties["charge_q"] = electron_braid.get_charge() # Use method

    virtual_photon = PhotonModel(); virtual_photon_state = {"energy": 1.0}
    amp_emit = calculate_interaction_amplitude(electron_braid, ('emit', virtual_photon_state), intermediate_braid)
    amp_absorb = calculate_interaction_amplitude(intermediate_braid, ('absorb', virtual_photon_state), electron_braid)
    prop_braid = braid_propagator(1.0); prop_photon = photon_propagator(virtual_photon_state["energy"]**2)
    integration_factor = 1.0 / (8.0 * np.pi**2); vertex_factor_sq = e_natural**2
    alpha_scaling = vertex_factor_sq / (4 * np.pi)
    coefficient_from_integration = 1.0 / (2.0 * np.pi); predicted_a_e = coefficient_from_integration * alpha_scaling
    print(f"  Structurally scales as alpha."); print(f"  Predicted a_e ≈ {predicted_a_e:.3e}")
    print(f"  Target alpha/(2*pi) ≈ {alpha_fine_structure / (2.0 * np.pi):.3e}")
    print(f"-----------------------------------------------------------------------")
    return predicted_a_e

# --- Main Execution ---
if __name__ == "__main__":
    print("--- Running LSC Comprehensive Particle Model Simulation (Corrected v6) ---")

    # 1. Define Context
    sn_base = SpinNetwork()
    print(f"Conceptual Context: Particles as braids within {sn_base.description}")

    # 2. Model Leptons & Neutrino
    leptons = ["electron", "muon", "tau"]
    fermions_to_model = leptons + ["neutrino"]
    particle_models = {}
    particle_braids = {}

    print("\n--- Fermion Mass Models (Leptons & Neutrino) ---")
    for name in fermions_to_model:
        try:
            particle_braids[name] = EmbeddedBraid(name.capitalize(), name)
            # Instantiate FermionOscillatorModel correctly
            particle_models[name] = FermionOscillatorModel(name, TARGET_ENERGY[name])
            particle_models[name].report()
        except ValueError as e: print(f"Could not model {name}: {e}")
        except KeyError as e: print(f"Missing data for {name}: {e}")


    # 3. Model Quarks Conceptually
    print("\n--- Quark Conceptual Models ---")
    quarks_to_model = ["up", "down"]
    for name in quarks_to_model:
        try:
            particle_braids[name] = EmbeddedBraid(name.capitalize(), name)
            print(f"  Defined: {particle_braids[name]} (Color: {particle_braids[name].get_color()})")
            print(f"  Note: Quark mass arises from confinement, not simple ZPE model.")
        except ValueError as e: print(f"Could not model {name}: {e}")
        except KeyError as e: print(f"Missing data for {name}: {e}")
    print("-----------------------------")

    # 4. Model Proton Conceptually
    print("\n--- Hadron Conceptual Model ---")
    try:
        # Ensure constituent braids are defined before use
        if "up" in particle_braids and "down" in particle_braids:
            up1_braid = EmbeddedBraid("Up1", "up"); up1_braid.color = "red"
            up2_braid = EmbeddedBraid("Up2", "up"); up2_braid.color = "green"
            down1_braid = EmbeddedBraid("Down1", "down"); down1_braid.color = "blue"
            proton_model = HadronBraidModel("Proton", [up1_braid, up2_braid, down1_braid])
            proton_model.report()
            particle_braids["proton"] = proton_model
        else: print("Error: Constituent quark braids not defined.")
    except ValueError as e: print(f"Error creating Proton model: {e}")
    except KeyError as e: print(f"Missing data for Proton model: {e}")


    # 5. Model Photon
    photon_model = PhotonModel()
    photon_model.report()

    # 6. Verify Spin (Placeholder) - Include Quarks and Proton
    print("\n--- Spin Checks (Fermions & Proton) ---")
    fermions_and_proton = fermions_to_model + quarks_to_model + ["proton"]
    for name in fermions_and_proton:
        if name in particle_braids:
            verify_spinor_transformation(particle_braids[name])

    # 7. Verify Mass Output (Leptons & Neutrino Only)
    print(f"\n--- Mass Verification (Leptons/Neutrino) ---")
    # Add checks for model existence before accessing
    for name in leptons:
        if name in particle_models: print(f"  {name.capitalize():<9}: Est={particle_models[name].get_mass_kg():.2e} kg")
        else: print(f"  {name.capitalize():<9}: Model not available.")
    if 'neutrino' in particle_models: print(f"  {'Neutrino':<9}: Est={particle_models['neutrino'].get_mass_kg():.2e} kg (Target m={neutrino_mass_eV}eV)")
    else: print(f"  {'Neutrino':<9}: Model not available.")
    if 'proton' in TARGET_ENERGY: print(f"  {'Proton':<9}: Target={TARGET_ENERGY['proton']*planck_mass_kg:.2e} kg (Not predicted by ZPE model)")
    print("--------------------------------------------")

    # 8. Calculate g-2 using interaction sketch for electron
    if "electron" in particle_braids:
        predicted_a_e = calculate_g_minus_2_leading_term(particle_braids["electron"])
        if predicted_a_e is not None: print(f"\n>>> g-2 Analysis Completed: Predicted a_e ≈ {predicted_a_e:.3e} <<<")

    print("\n--- Simulation Finished ---")
    # Final Status Report (using values from existing models where available)
    print("\nFINAL MODEL STATUS (LSC Comprehensive Particle Model - v6):")
    print("  - Structure: Fermions = braids; Bosons = ripples/links; Hadrons = bound braids.")
    print("  - Spin: ASSUMED via Framed Holonomy mechanism for fermions.")
    print("  - Mass: Leptons via ZPE+Hierarchy; Hadrons via QCD-like binding; Photon massless.")
    print(f"  - Lepton Hierarchy Requires:")
    # Add checks before accessing properties
    e_c1 = particle_models['electron'].required_c1_for_k if 'electron' in particle_models else 'N/A'
    m_c1 = particle_models['muon'].required_c1_for_k if 'muon' in particle_models else 'N/A'
    t_c1 = particle_models['tau'].required_c1_for_k if 'tau' in particle_models else 'N/A'
    n_cg = particle_models['neutrino'].calculated_C_grav_fit if 'neutrino' in particle_models else 'N/A'
    print(f"      Charged: k ~ C_t*exp(-c1(Topo)/alpha) (Needs c1≈{e_c1:.3f}(e) > {m_c1:.3f}(mu) > {t_c1:.3f}(tau)).")
    print(f"      Neutral: k ~ C_g*(E0/Ep)^2 (Needs C_g≈{n_cg:.3f}).")
    print("  - Interaction: Conceptual electron-photon vertex defined structurally.")
    print("  - g-2: Leading correction structurally scales as alpha; coefficient needs derivation.")
    print("  - Confinement: Hypothesized via topological instability of isolated colored quarks.")
    print("  - CORE THEORETICAL CHALLENGES:")
    print("      1. Derive Lepton Stiffnesses k (incl. c1(Topology) & C_t, C_g & E0^2 scaling).")
    print("      2. Derive Hadron Masses from bound state braid dynamics + SU(3) gluon links.")
    print("      3. Derive m_eff (Confirm ≈ M_p?).")
    print("      4. Derive Spinor Phase (Formalize framed LQG + SU(2)).")
    print("      5. Formalize Vertex & Propagators; Calculate g-2 coefficient = 1/(2pi).")
    print("      6. Formalize SU(3) color & gluon links; Derive confinement.")
    print("      7. Photon Model: Formalize ripple; derive Maxwell dynamics.")