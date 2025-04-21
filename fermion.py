# --- LSC Comprehensive Particle Model (Corrected) ---
# --- Includes Leptons (e, mu, tau), Neutrino (generic), Photon ---
# --- Mass = ZPE of Braid Oscillation (Different k origins for charged/neutral) ---
# --- Spin = Framed Holonomy Mechanism (Assumed) ---
# --- Interaction Vertex Structure & g-2 Sketch Included ---
# --- Single File Python Script ---

import networkx as nx
import numpy as np
import math

print("--- Loading LSC Comprehensive Particle Model ---")
print("--- MODEL: Fermion Mass = ZPE of Braid Oscillation ---")
print("---        Stiffness 'k' arises from Cancellation Residual ---")
print("---        Hypothesized Origin depends on charge ---")
print("---        Spin arises from Framed Braid Rotation ---")
print("---        Interaction Vertex Structure Included ---")


# --- Physical Constants (SI) ---
hbar_si = 1.0545718e-34 # J*s
c_si = 2.99792458e8    # m/s
G_si = 6.67430e-11     # m^3 kg^-1 s^-2
# Masses (approximate, kg)
electron_mass_kg = 9.1093837e-31
muon_mass_kg = 1.8835316e-28 # approx 105.66 MeV/c^2
tau_mass_kg = 3.16754e-27    # approx 1776.86 MeV/c^2
neutrino_mass_eV = 0.1 # Representative value
eV_to_J = 1.602176634e-19
neutrino_mass_kg = (neutrino_mass_eV * eV_to_J) / (c_si**2)
# Fine structure constant
alpha_fine_structure = 1 / 137.035999
elementary_charge_C = 1.60217663e-19 # Magnitude

# --- Planck Units ---
planck_length_si = np.sqrt(G_si * hbar_si / (c_si**3))
planck_mass_kg = np.sqrt(hbar_si * c_si / G_si)
planck_energy_J = planck_mass_kg * c_si**2
planck_time_si = planck_length_si / c_si
# Natural Planck Units (hbar=c=1)
planck_stiffness_natural = 1.0 # Units of E_p / lp^2
planck_mass_natural = 1.0      # Units of M_p
# Elementary charge in natural units (Heaviside-Lorentz, hbar=c=1)
e_natural = np.sqrt(4 * np.pi * alpha_fine_structure)

# --- Target Energies (Natural Units E/Ep) ---
TARGET_ENERGY = {
    "electron": (electron_mass_kg * c_si**2) / planck_energy_J,
    "muon":     (muon_mass_kg * c_si**2) / planck_energy_J,
    "tau":      (tau_mass_kg * c_si**2) / planck_energy_J,
    "neutrino": (neutrino_mass_kg * c_si**2) / planck_energy_J,
}

# --- Particle Properties Lookup ---
# Define this *before* FermionOscillatorModel uses it
PARTICLE_PROPS = {
    "electron": {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True},
    "muon":     {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True},
    "tau":      {"charge_q": -1, "spin_j": 0.5, "is_lepton": True, "is_charged": True},
    "neutrino": {"charge_q":  0, "spin_j": 0.5, "is_lepton": True, "is_charged": False},
    "photon":   {"charge_q":  0, "spin_j": 1.0, "is_lepton": False,"is_charged": False}
}

print(f"--- Constants ---")
print(f"Planck Mass (kg): {planck_mass_kg:.2e}")
print(f"Alpha: {alpha_fine_structure:.3e}")
print(f"e (natural): {e_natural:.3f}")
print(f"--- Target Energies (Planck Units) ---")
for p, E in TARGET_ENERGY.items(): print(f"  {p.capitalize():<9}: {E:.3e}")
print("-------------------------------------\n")


# --- Effective Oscillator Model for Fermions ---
class FermionOscillatorModel:
    """ Models fermion mass as ZPE of braid oscillation. """
    def __init__(self, particle_name, target_energy_natural): # Removed is_charged argument
        particle_name_lower = particle_name.lower()
        if particle_name_lower not in TARGET_ENERGY:
            raise ValueError(f"Unknown particle: {particle_name}")
        self.particle_name = particle_name_lower
        self.properties = PARTICLE_PROPS[self.particle_name]
        self.is_charged = self.properties["is_charged"] # Determine internally
        self.target_E0_natural = target_energy_natural
        # Assumptions
        self.m_eff_natural = 1.0 # Assume Planckian inertia
        self.assumed_C_tunnel = 1.0 # Prefactor for charged k interpretation
        self.calculated_C_grav_fit = None # Coefficient for neutral k interpretation
        # Calculated
        self.omega_natural = None
        self.k_required_natural = None
        self.required_c1_for_k = None # Only relevant for charged

        self.calculate_required_parameters()
        self.interpret_k_origin()

    def calculate_required_parameters(self):
        """Calculate omega and k needed to match target E0, assuming m_eff=1."""
        if self.target_E0_natural <= 0: return
        self.omega_natural = 2 * self.target_E0_natural # E0 = 0.5 * omega
        self.k_required_natural = self.m_eff_natural * (self.omega_natural**2) # k = m * omega^2

    def interpret_k_origin(self):
        """Interpret the required k using the relevant hypothesis."""
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
    def get_mass_kg(self): return self.target_E0_natural * planck_mass_kg

    def report(self):
        """Reports the calculated required parameters and their interpretation."""
        print(f"--- {self.particle_name.upper()} Oscillator Parameters ---")
        print(f"  Target Mass (kg): {self.get_mass_kg():.2e}")
        print(f"  Target E0 (Planck Units): {self.target_E0_natural:.3e}")
        if self.omega_natural and self.k_required_natural is not None:
            print(f"  Required Omega (Planck Units): {self.omega_natural:.3e}")
            print(f"  >> Required Stiffness k (Planck Units): {self.get_required_stiffness():.3e}")
        else: print("  Omega/k: Invalid")
        print(f"--- Interpretation of Required Stiffness k ---")
        if self.is_charged:
             print(f"  Charged Hypothesis: k ≈ C_tunnel * exp(-c1/alpha)")
             if self.required_c1_for_k is not None:
                 print(f"  >> Needs c1 ≈ {self.required_c1_for_k:.3f} (if C_t={self.assumed_C_tunnel:.1f})")
             else: print(f"  Could not calculate required 'c1'.")
        else: # Neutrino
             print(f"  Neutral Hypothesis: k ≈ C_grav * (E0/Ep)^2")
             if self.calculated_C_grav_fit is not None:
                 print(f"  >> Needs C_g ≈ {self.calculated_C_grav_fit:.3f}")
             else: print(f"  Could not calculate required 'C_g'.")
        print(f"---------------------------------------------")

# --- Spin Network / Braid Placeholders ---
class SpinNetwork:
    def __init__(self): self.description = "Conceptual LQG Background"
class EmbeddedBraid:
    def __init__(self, name, particle_type):
        self.name = name; self.particle_type=particle_type.lower() # Ensure lowercase for lookup
        self.properties = PARTICLE_PROPS[self.particle_type]
        self.description=f"Conceptual j={self.get_spin()} framed braid ({particle_type})"
    def get_spin(self): return self.properties["spin_j"]
    def get_charge(self): return self.properties["charge_q"]

# --- Spin Transformation Placeholder ---
def verify_spinor_transformation(braid_structure):
    """Checks if the placeholder phase function yields spinor properties."""
    print(f"\n--- Verifying Spinor Transformation for '{braid_structure.name}' ---")
    spin_j = braid_structure.get_spin()
    if not np.isclose(spin_j, 0.5): print(f"  Not Spin-1/2 (j={spin_j}). Skipping check."); return None
    phase_2pi = np.exp(-1j * 2*np.pi * 0.5); phase_4pi = np.exp(-1j * 4*np.pi * 0.5)
    result = np.isclose(phase_2pi, -1.0) and np.isclose(phase_4pi, 1.0)
    print(f"Consistent with Spin-1/2: {result} (assuming mechanism)")
    print("----------------------------------------------------")
    return result

# --- Photon Model Placeholder ---
class PhotonModel:
    """Conceptual model for the photon."""
    def __init__(self):
        self.name = "Photon"
        self.properties = PARTICLE_PROPS["photon"]
        self.description = "Massless, spin-1, chargeless propagating ripple in connection field"
    def get_spin(self): return self.properties["spin_j"]
    def get_charge(self): return self.properties["charge_q"]
    def report(self):
        print(f"\n--- PHOTON MODEL ---")
        print(f"  Description: {self.description}")
        print(f"  Spin: {self.get_spin()}, Charge: {self.get_charge()}e, Mass: 0")
        print(f"--------------------")

# --- Interaction Vertex & g-2 Sketch Placeholders ---
# (Functions: braid_current_operator, photon_field_operator, su2_coupling_factor,
# calculate_vertex_amplitude, braid_propagator, photon_propagator,
# calculate_g_minus_2_leading_term remain the same as previous correct version)
def braid_current_operator(braid_state_in, braid_state_out): return braid_state_in.get_charge()
def photon_field_operator(photon_state, creation=False): return 1.0
def su2_coupling_factor(spin_in, spin_photon, spin_out):
    if np.isclose(spin_in, 0.5) and np.isclose(spin_photon, 1.0) and np.isclose(spin_out, 0.5): return 1.0
    return 0.0
def calculate_vertex_amplitude(braid_in, photon_action, braid_out):
    photon_state = photon_action[1]; is_emission = (photon_action[0] == 'emit')
    coupling_g = e_natural * braid_in.get_charge();
    if np.isclose(coupling_g, 0.0): return 0.0
    cg_factor = su2_coupling_factor(braid_in.get_spin(), 1.0, braid_out.get_spin());
    if np.isclose(cg_factor, 0.0): return 0.0
    current_factor = braid_current_operator(braid_in, braid_out)
    photon_op_factor = photon_field_operator(photon_state, creation=is_emission)
    amplitude = coupling_g * cg_factor * current_factor * photon_op_factor * complex(0.8, 0.2)
    return amplitude
def braid_propagator(energy_diff_natural): return 1.0 / (abs(energy_diff_natural) + 1e-30)
def photon_propagator(photon_momentum_sq_natural): return 1.0 / (abs(photon_momentum_sq_natural) + 1e-30)
def calculate_g_minus_2_leading_term(electron_braid):
    print("\n--- Calculating Electron Anomalous Magnetic Moment (g-2) - Sketch ---")
    if np.isclose(electron_braid.get_charge(), 0.0): print("ERROR: g-2 requires charged particle."); return None
    intermediate_braid = EmbeddedBraid("Electron_Intermediate", "electron")
    virtual_photon = PhotonModel(); virtual_photon_state = {"energy": 1.0}
    amp_emit = calculate_vertex_amplitude(electron_braid, ('emit', virtual_photon_state), intermediate_braid)
    amp_absorb = calculate_vertex_amplitude(intermediate_braid, ('absorb', virtual_photon_state), electron_braid)
    prop_braid = braid_propagator(1.0)
    prop_photon = photon_propagator(virtual_photon_state["energy"]**2)
    integration_factor = 1.0 / (8.0 * np.pi**2) # Assumed factor
    vertex_factor_sq = e_natural**2
    alpha_scaling = vertex_factor_sq / (4 * np.pi)
    coefficient_from_integration = 1.0 / (2.0 * np.pi) # Assumed coeff
    predicted_a_e = coefficient_from_integration * alpha_scaling
    print(f"  Structurally scales as alpha.")
    print(f"  Predicted a_e ≈ {predicted_a_e:.3e} (Assuming LQG integral yields 1/(2pi) coeff)")
    print(f"  Target alpha/(2*pi) ≈ {alpha_fine_structure / (2.0 * np.pi):.3e}")
    print(f"-----------------------------------------------------------------------")
    return predicted_a_e

# --- Main Execution ---
if __name__ == "__main__":
    print("--- Running LSC Comprehensive Particle Model Simulation (Corrected) ---")

    # 1. Define Context
    sn_base = SpinNetwork()
    print(f"Conceptual Context: Particles as braids within {sn_base.description}")

    # 2. Model Leptons & Neutrino
    leptons = ["electron", "muon", "tau"]
    fermions_to_model = leptons + ["neutrino"]
    particle_models = {}
    particle_braids = {}

    for name in fermions_to_model:
        particle_braids[name] = EmbeddedBraid(name.capitalize(), name) # Use lowercase name for lookup
        # Corrected Instantiation Call:
        particle_models[name] = FermionOscillatorModel(name, TARGET_ENERGY[name])
        particle_models[name].report()

    # 3. Model Photon
    photon_model = PhotonModel()
    photon_model.report()

    # 4. Verify Spin (Placeholder)
    print("\n--- Spin Checks ---")
    for name in fermions_to_model:
        verify_spinor_transformation(particle_braids[name])

    # 5. Verify Mass Output (By Construction for model consistency)
    print(f"\n--- Mass Verification ---")
    for name in leptons: # Check charged leptons
        print(f"  {name.capitalize():<9}: Est={particle_models[name].get_mass_kg():.2e} kg, Actual={TARGET_ENERGY[name]*planck_mass_kg:.2e} kg")
    # Check neutrino
    print(f"  {'Neutrino':<9}: Est={particle_models['neutrino'].get_mass_kg():.2e} kg (Target m={neutrino_mass_eV}eV)")
    print("----------------------")

    # 6. Calculate g-2 using interaction sketch for electron
    predicted_a_e = calculate_g_minus_2_leading_term(particle_braids["electron"]) # Pass the electron braid object
    if predicted_a_e is not None:
        print(f"\n>>> g-2 Analysis Completed <<<")
        print(f"  Predicted leading term a_e ≈ {predicted_a_e:.3e}")

    print("\n--- Simulation Finished ---")
    print("\nFINAL MODEL STATUS (LSC Fermion Model + Interaction - Corrected):")
    print("  - Structure: Fermions = Dynamic, framed, j=1/2 braid excitations.")
    print("  - Spin: ASSUMED via Framed Holonomy mechanism.")
    print("  - Mass: Modeled as ZPE of braid oscillation (m = E0 = 0.5*omega).")
    print(f"  - Hierarchy requires different stiffness origins:")
    print(f"      Charged Leptons: Requires k ~ C_t*exp(-c1/alpha) (Needs c1≈{particle_models['electron'].required_c1_for_k:.3f}(e), {particle_models['muon'].required_c1_for_k:.3f}(mu), {particle_models['tau'].required_c1_for_k:.3f}(tau) if C_t=1).")
    print(f"      Neutrino (neutral): Requires k ~ C_g*(E0/Ep)^2 (Needs C_g≈{particle_models['neutrino'].calculated_C_grav_fit:.3f}).")
    print("  - Interaction: Conceptual vertex V ~ e * CG * J_braid * A_photon defined structurally.")
    print("  - g-2: Leading correction structurally scales as alpha; coefficient 1/(2pi) needs derivation.")
    print("  - CORE THEORETICAL CHALLENGES:")
    print("      1. Electron k: Derive mechanism yielding k ~ exp(-0.74/alpha) (incl. c1 & prefactor C_t).")
    print("      2. Neutrino k: Derive mechanism yielding k ~ 4.0*(E0/Ep)^2 (incl. coefficient C_g & scaling).")
    print("      3. Derive effective inertia 'm_eff' (Confirm ≈ M_p for all?).")
    print("      4. Derive spinor transformation phase rigorously (Formalize framed LQG).")
    print("      5. Formalize Vertex & Propagators; Calculate g-2 coefficient = 1/(2pi).")
    print("      6. Extend to Quarks/Hadrons (Incorporate SU(3) color).")
    print("      7. Photon Model: Formalize as geometric excitation; derive Maxwell dynamics.")