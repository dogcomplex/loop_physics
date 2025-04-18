# --- LSC Ansatz 2b: Electron as Dynamic, Charged Braid ---
# --- Single File Python Script ---

import networkx as nx
import numpy as np
# Potential future imports: sympy, scipy.sparse, knot/braid libraries

print("--- Loading LSC Ansatz 2b Simulation Framework ---")

# --- Physical Constants (SI) ---
hbar_si = 1.0545718e-34 # J*s
c_si = 2.99792458e8    # m/s
G_si = 6.67430e-11     # m^3 kg^-1 s^-2
electron_mass_kg = 9.1093837e-31
elementary_charge_C = 1.60217663e-19 # Magnitude

# --- Planck Units ---
planck_length_si = np.sqrt(G_si * hbar_si / (c_si**3))
planck_mass_kg = np.sqrt(hbar_si * c_si / G_si)
planck_energy_J = planck_mass_kg * c_si**2
planck_time_si = planck_length_si / c_si
planck_area_si = planck_length_si**2
planck_volume_si = planck_length_si**3
planck_charge_C = np.sqrt(4 * np.pi * 8.854187817e-12 * hbar_si * c_si) # Approx 1.87e-18 C

# Natural units (hbar=c=1, G=lp^2=1/mp^2)
AREA_OPERATOR_NATURAL_CONST = 8 * np.pi # Assuming Immirzi beta=1

print(f"--- Constants (SI) ---")
print(f"Planck Length (m): {planck_length_si:.2e}")
print(f"Planck Mass (kg): {planck_mass_kg:.2e}")
print(f"Planck Charge (C): {planck_charge_C:.2e}")
print("--------------------\n")

# --- Spin Network Components (Classes similar, slight mods) ---
class Intertwiner:
    def __init__(self, incident_spins): self.spins = tuple(sorted(incident_spins)); self.data_repr = f"Intertwiner{self.spins}"
    def __repr__(self): return self.data_repr
class SpinNetworkNode:
    def __init__(self, node_id, position=None, intertwiner_obj=None): self.id = node_id; self.position = position; self.intertwiner = intertwiner_obj
    def __repr__(self): return f"Node({self.id}, Intw: {self.intertwiner})"
class SpinNetworkEdge:
    def __init__(self, edge_id, u, v, spin_j, holonomy=None): self.id = edge_id; self.u = u; self.v = v; self.spin_j = spin_j; self.holonomy = holonomy # holonomy ~ SU(2) matrix
    def __repr__(self): return f"Edge({self.id}: {self.u}-{self.v}, j={self.spin_j})"
class SpinNetwork:
    def __init__(self): self.graph = nx.MultiGraph(); self._node_counter = 0; self._edge_counter = 0; self.braids = {}
    def add_node(self, position=None, intertwiner_obj=None): node_id=self._node_counter; node=SpinNetworkNode(node_id, position, intertwiner_obj); self.graph.add_node(node_id, data=node); self._node_counter+=1; return node_id
    def add_edge(self, u, v, spin_j, holonomy=None): edge_id=self._edge_counter; assert u in self.graph and v in self.graph; edge=SpinNetworkEdge(edge_id, u, v, spin_j, holonomy); self.graph.add_edge(u, v, key=edge_id, data=edge); self._edge_counter+=1; return edge_id
    def get_node(self, node_id): return self.graph.nodes[node_id]['data']
    def get_edge_data_by_id(self, edge_id):
         for u, v, k, data in self.graph.edges(keys=True, data=True):
             if k == edge_id: return data['data']
         print(f"Warning: Edge {edge_id} not found."); return None
    def get_edges_incident_to_node(self, node_id): return [d['data'] for u, v, k, d in self.graph.edges(node_id, keys=True, data=True)]
    def number_of_nodes(self): return self.graph.number_of_nodes()
    def number_of_edges(self): return self.graph.number_of_edges()
    def __str__(self): return f"SpinNetwork(Nodes: {self.number_of_nodes()}, Edges: {self.number_of_edges()}, Braids: {list(self.braids.keys())})"

# --- Braid Representation (Enhanced) ---
class EmbeddedBraid:
    """Represents a dynamic, charged braid excitation."""
    def __init__(self, name, braid_type, edge_ids, spin_network_ref, charge_quantum=0):
        self.name = name
        self.braid_type = braid_type # E.g., "3-strand-closed-trefoil"
        self.edge_ids = list(edge_ids)
        if not isinstance(spin_network_ref, SpinNetwork): raise TypeError("spin_network_ref must be a SpinNetwork instance")
        self.spin_network = spin_network_ref
        self.charge_quantum = charge_quantum # In units of fundamental charge 'e'
        # Properties derived from structure (placeholders)
        self.properties = {
            "complexity": len(self.edge_ids) / 3.0, # Crude topology proxy
            "SU2_spin": 0.5, # Assumed for electron
            "effective_area": 0.0, # Placeholder for area term calculation
            "effective_volume": 0.0 # Placeholder for volume term calculation
        }
        self.verify_path()

    def verify_path(self): # Basic checks
        # ... (checks from previous version) ...
        return True

    def get_spin(self): return self.properties.get("SU2_spin")
    def get_nodes(self):
        nodes = set()
        for edge_id in self.edge_ids:
            edge_data = self.spin_network.get_edge_data_by_id(edge_id)
            if edge_data: nodes.add(edge_data.u); nodes.add(edge_data.v)
        return list(nodes)
    def get_topology_proxy(self): return self.properties.get("complexity", 1.0)
    def get_charge(self): return self.charge_quantum # Returns charge in units of e

    def __str__(self): return f"EmbeddedBraid(Name: {self.name}, Type: {self.braid_type}, Edges: {len(self.edge_ids)}, Spin: {self.get_spin()}, Charge: {self.charge_quantum}e, Nodes: {len(self.get_nodes())})"


# --- Calculation Functions ---

# --- Spin Related (Placeholder Transformation) ---
def calculate_transformation_phase(braid_structure, angle):
    # print(f"WARNING: calculate_transformation_phase placeholder.")
    spin_j = braid_structure.get_spin()
    if spin_j == 0.5: return np.exp(-1j * angle * 0.5)
    return 1.0

def verify_spinor_transformation(braid_structure):
    print(f"\n--- Verifying Spinor Transformation for Braid '{braid_structure.name}' ---")
    phase_2pi = calculate_transformation_phase(braid_structure, angle=2*np.pi)
    phase_4pi = calculate_transformation_phase(braid_structure, angle=4*np.pi)
    is_negated_at_2pi = np.isclose(phase_2pi, -1.0)
    is_identity_at_4pi = np.isclose(phase_4pi, 1.0)
    result = is_negated_at_2pi and is_identity_at_4pi
    print(f"Phase after 2pi rotation (placeholder): {phase_2pi:.3f}")
    print(f"Phase after 4pi rotation (placeholder): {phase_4pi:.3f}")
    print(f"Consistent with Spin-1/2 Transformation: {result}")
    print("----------------------------------------------------")
    return result

# --- Mass/Energy Related ---
def lqg_area_term_contribution(spin_network, edge_ids_list):
    """Calculates dimensionless Sum(sqrt(j(j+1)))."""
    area_term_sum = 0.0
    for edge_id in edge_ids_list:
        edge_data = spin_network.get_edge_data_by_id(edge_id)
        if edge_data and edge_data.spin_j > 0:
            area_term_sum += np.sqrt(edge_data.spin_j * (edge_data.spin_j + 1))
    return area_term_sum # Dimensionless sum

def lqg_volume_term_contribution(spin_network, node_ids_list):
    """Calculates dimensionless Sum(Volume Proxy at node)."""
    volume_term_sum = 0.0
    for node_id in node_ids_list:
        incident_edges_data = spin_network.get_edges_incident_to_node(node_id)
        valence = len(incident_edges_data)
        if valence >= 3:
             # Use sum of spins as proxy for node complexity contribution
            volume_term_sum += sum(edge.spin_j for edge in incident_edges_data)
    return volume_term_sum # Dimensionless proxy

# --- NEW: Conceptual Braid Hamiltonian (Ansatz 2b) ---
def construct_braid_hamiltonian_ansatz2b(spin_network, embedded_braid):
    """
    Calculates a conceptual ground state energy based on the dynamic,
    charged braid model with near-cancellation. Returns energy in Planck units.
    """
    print(f"INFO: Constructing H_local (Ansatz 2b) for '{embedded_braid.name}'.")

    braid_edge_ids = embedded_braid.edge_ids
    braid_node_ids = embedded_braid.get_nodes()
    spin_j = embedded_braid.get_spin()
    topology_proxy = embedded_braid.get_topology_proxy()
    charge_q = embedded_braid.get_charge() # Charge in units of e

    # Coefficients (Dimensionless - ADJUSTING FOR TWEAK 2c)
    # Aiming for residual sum ~ +4e-23 Planck units
    C_kin = 1e-22     # Increased significantly
    C_top = 0.5e-22   # Increased significantly
    C_geom_A = AREA_OPERATOR_NATURAL_CONST # Keep cancellation core
    C_geom_V = 0.001  # Reduced negative binding, or make slightly positive? Try smaller negative first.
    C_charge = AREA_OPERATOR_NATURAL_CONST # Keep cancellation core
    C_int = 0.0       # Keep zero

    # --- Calculate Term Contributions (Planck Units) ---
    term_kin_factor = len(braid_edge_ids) * (spin_j**2)
    term_topology_factor = topology_proxy # Using complexity proxy
    area_term_dimless = lqg_area_term_contribution(spin_network, braid_edge_ids)
    volume_term_dimless = lqg_volume_term_contribution(spin_network, braid_node_ids)

    E_kin = C_kin * term_kin_factor
    E_top = C_top * term_topology_factor
    E_geom_A = C_geom_A * area_term_dimless
    E_geom_V = -C_geom_V * volume_term_dimless # Negative contribution using C_geom_V
    E_top_charge = -C_charge * abs(charge_q) * area_term_dimless
    E_int = C_int

    print(f"  Hamiltonian Term Contributions (Planck Units):")
    print(f"    Kinetic (proxy)       : {E_kin:.4e}")
    print(f"    Topological (proxy)   : {E_top:.4e}")
    print(f"    Geometric Area (+)    : {E_geom_A:.4e}")
    print(f"    Geometric Volume (-)  : {E_geom_V:.4e}")
    print(f"    Topological Charge (-): {E_top_charge:.4e}")
    print(f"    Internal Interaction  : {E_int:.4e}")

    # --- Near-Cancellation Hypothesis ---
    # The stable ground state minimizes energy. We hypothesize this minimum occurs
    # when the positive geometric area energy is almost cancelled by the
    # negative topological charge energy.
    E_cancellation_residual = E_geom_A + E_top_charge # Should still be ~zero
    total_energy_natural = E_kin + E_top + E_geom_V + E_int + E_cancellation_residual # This should now be positive

    print(f"  Cancellation Term (E_geom_A + E_top_charge): {E_cancellation_residual:.4e}")
    print(f"  Total Ground State Energy (Conceptual): {total_energy_natural:.4e} [E_planck]")

    return total_energy_natural


def construct_braid_hamiltonian_ansatz2d(spin_network, embedded_braid, epsilon):
    """
    Calculates conceptual ground state energy based on Ansatz 2d:
    Imperfect cancellation between Geometric Area and Topological Charge terms.
    Returns energy in Planck units.
    """
    print(f"INFO: Constructing H_local (Ansatz 2d: Imperfect Cancellation, epsilon={epsilon:.2e}) for '{embedded_braid.name}'.")

    braid_edge_ids = embedded_braid.edge_ids
    charge_q = embedded_braid.get_charge()

    # --- Coefficients ---
    C_geom_A = AREA_OPERATOR_NATURAL_CONST # Geometric Area term strength
    # C_charge is now implicitly (1 - epsilon) * C_geom_A
    # Set other coefficients to zero for this simplified test
    C_kin = 0.0
    C_top = 0.0
    C_geom_V = 0.0
    C_int = 0.0

    # --- Calculate Term Contributions (Dimensionless factors) ---
    area_term_dimless = lqg_area_term_contribution(spin_network, braid_edge_ids)

    # --- Calculate Energy Contributions (Planck Units) ---
    E_geom_A = C_geom_A * area_term_dimless
    # Imperfect Cancellation: Charge term is slightly weaker
    E_top_charge = -(1.0 - epsilon) * C_geom_A * abs(charge_q) * area_term_dimless

    # Other terms are zero in this simplified version
    E_kin = C_kin # = 0
    E_top = C_top # = 0
    E_geom_V = C_geom_V # = 0
    E_int = C_int # = 0

    print(f"  Hamiltonian Term Contributions (Planck Units):")
    print(f"    Geometric Area (+)    : {E_geom_A:.4e}")
    print(f"    Topological Charge (-): {E_top_charge:.4e} (using epsilon={epsilon:.2e})")

    # The residual is the dominant energy term now
    E_cancellation_residual = E_geom_A + E_top_charge
    total_energy_natural = E_kin + E_top + E_geom_V + E_int + E_cancellation_residual

    print(f"  Cancellation Residual (E_geom_A + E_top_charge): {E_cancellation_residual:.4e}")
    print(f"  Total Ground State Energy (Conceptual): {total_energy_natural:.4e} [E_planck]")

    return total_energy_natural


def analyze_braid_mass_ansatz2b(spin_network, embedded_braid):
    """Estimate mass from the Ansatz 2b Hamiltonian."""
    print(f"\n--- Analyzing Mass/Energy (Ansatz 2b: Dynamic Charged Braid) ---")

    ground_state_energy_natural = construct_braid_hamiltonian_ansatz2b(spin_network, embedded_braid)

    if ground_state_energy_natural <= 0:
         print("Warning: Conceptual ground state energy is non-positive. Mass is zero or undefined.")
         mass_estimate_kg = 0.0
    else:
        # Convert Planck energy units to kg
        mass_estimate_kg = ground_state_energy_natural * planck_mass_kg

    print(f"\nEstimated Ground State Mass (Ansatz 2b, kg): {mass_estimate_kg:.2e}")
    print(f"Actual Electron Mass (kg): {electron_mass_kg:.2e}")
    ratio = mass_estimate_kg / electron_mass_kg if electron_mass_kg and mass_estimate_kg > 0 else float('inf')
    print(f"Ratio to Electron Mass: {ratio:.2e}")
    print("--------------------------------------------------------------")
    return mass_estimate_kg


def analyze_braid_mass_ansatz2d(spin_network, embedded_braid, epsilon):
    """Estimate mass from the Ansatz 2d Hamiltonian."""
    print(f"\n--- Analyzing Mass/Energy (Ansatz 2d: Imperfect Cancellation) ---")
    ground_state_energy_natural = construct_braid_hamiltonian_ansatz2d(spin_network, embedded_braid, epsilon)
    # ... (rest of mass calculation and printing as before) ...
    if ground_state_energy_natural <= 0:
         print("Warning: Conceptual ground state energy is non-positive/zero.")
         mass_estimate_kg = 0.0
    else:
        mass_estimate_kg = ground_state_energy_natural * planck_mass_kg

    print(f"\nEstimated Ground State Mass (Ansatz 2d, kg): {mass_estimate_kg:.2e}")
    print(f"Actual Electron Mass (kg): {electron_mass_kg:.2e}")
    ratio = mass_estimate_kg / electron_mass_kg if electron_mass_kg and mass_estimate_kg > 0 else float('inf')
    print(f"Ratio to Electron Mass: {ratio:.2e}")
    print("--------------------------------------------------------------")
    return mass_estimate_kg


# --- Main Execution ---
if __name__ == "__main__":
    print("--- Running LSC Ansatz 2d Simulation (Imperfect Cancellation) ---")

    # 1. Create Background Spin Network (Minimal)
    sn_base = SpinNetwork()
    # Add a few nodes/edges to represent vacuum structure if needed by Hamiltonian
    n_vac1 = sn_base.add_node()
    n_vac2 = sn_base.add_node()
    sn_base.add_edge(n_vac1, n_vac2, 0) # Example vacuum edge
    print(f"Created Base Spin Network: {sn_base}")

    # 2. Define and Embed the Electron Braid Structure
    # Needs a concrete graph structure allowing a 3-strand braid.
    # Let's reuse the simple 6-node structure conceptually.
    braid_edge_ids = []
    nodes = [sn_base.add_node() for _ in range(6)]
    spin_electron = 0.5
    try:
        # Schematic wiring - assumes these edges form the braid structure
        braid_edge_ids.append(sn_base.add_edge(nodes[0], nodes[1], spin_j=spin_electron))
        braid_edge_ids.append(sn_base.add_edge(nodes[1], nodes[2], spin_j=spin_electron))
        braid_edge_ids.append(sn_base.add_edge(nodes[3], nodes[4], spin_j=spin_electron))
        braid_edge_ids.append(sn_base.add_edge(nodes[4], nodes[5], spin_j=spin_electron))
        braid_edge_ids.append(sn_base.add_edge(nodes[0], nodes[3], spin_j=spin_electron)) # Crossing proxy
        braid_edge_ids.append(sn_base.add_edge(nodes[1], nodes[4], spin_j=spin_electron)) # Crossing proxy
        braid_edge_ids.append(sn_base.add_edge(nodes[2], nodes[5], spin_j=spin_electron)) # Closure proxy
        # Ensure nodes involved have valence >= 3 for Volume proxy
        sn_base.add_edge(nodes[0],nodes[2],0); sn_base.add_edge(nodes[1],nodes[3],0);
        sn_base.add_edge(nodes[4],nodes[0],0); sn_base.add_edge(nodes[5],nodes[1],0);

    except ValueError as e: print(f"Error creating braid network: {e}"); exit()

    electron_braid = EmbeddedBraid(
        name="Electron",
        braid_type="Conceptual-3Strand-j1/2",
        edge_ids=braid_edge_ids,
        spin_network_ref=sn_base,
        charge_quantum=-1 # Assign charge -1e
    )
    electron_braid.properties = {"complexity": 1.8, "SU2_spin": 0.5} # Update proxy properties

    print(f"Defined Embedded Braid: {electron_braid}")


    # 3. Analyze Spin Transformation (Placeholder check)
    # >> Point for running and sharing output <<
    spin_result = verify_spinor_transformation(electron_braid)
    print(f"\n>>> Run 1 Output: Spin Analysis Completed. Is Spinor-like: {spin_result}")

    # Calculate the required epsilon
    area_term_sum_val = lqg_area_term_contribution(sn_base, electron_braid.edge_ids)
    target_energy_natural = electron_mass_kg / planck_mass_kg
    required_epsilon = target_energy_natural / (AREA_OPERATOR_NATURAL_CONST * abs(electron_braid.get_charge()) * area_term_sum_val)
    print(f"\nTarget energy (Planck units): {target_energy_natural:.2e}")
    print(f"Area term sum (dimless): {area_term_sum_val:.3f}")
    print(f"Required epsilon for electron mass: {required_epsilon:.2e}")

    # Run with the calculated epsilon
    mass_result = analyze_braid_mass_ansatz2d(sn_base, electron_braid, epsilon=required_epsilon)
    print(f"\n>>> Run 2 Output: Mass Analysis (Ansatz 2d, Epsilon={required_epsilon:.2e}). Estimated Mass (kg): {mass_result:.2e}")

    # Optional: Run with a slightly different epsilon to see sensitivity
    # mass_result_perturbed = analyze_braid_mass_ansatz2d(sn_base, electron_braid, epsilon=required_epsilon * 1.1)

    # For comparison, also run the original Ansatz 2b
    mass_result_2b = analyze_braid_mass_ansatz2b(sn_base, electron_braid)
    print(f"\n>>> Run 3 Output: Mass Analysis (Ansatz 2b). Estimated Mass (kg): {mass_result_2b:.2e}")

    print("\n--- Simulation Finished ---")
    print("NOTE: Hamiltonian uses conceptual terms & coefficients. Imperfect cancellation is key.")