# --- LSC Stability Filter (Relaxed - Including Links) ---
# --- Searching for Candidate Topologies (Knots & Links) for Leptons ---

import numpy as np
import math
import sys
import time

# Attempt import
try:
    import spherogram
    print("INFO: Spherogram library found.")
except ImportError:
    print("ERROR: Spherogram library not found. Please install: pip install spherogram")
    sys.exit(1)

print("\n--- Loading LSC Relaxed Stability Filter (Including Links) ---")
print("--- Applying basic topological filters to find candidate knots/links for leptons ---")

# --- Physical Constants & Targets ---
alpha_fine_structure = 1 / 137.035999
TARGET_ENERGY_ELECTRON = 4.18540e-23 # Planck Units
TARGET_ENERGY_MUON = 8.65429e-21
TARGET_ENERGY_TAU = 1.45532e-19
ASSUMED_C_TUNNEL = 1.0

# Function to calculate required c1 for a given energy E0
def get_required_c1(E0_natural, C_tunnel=1.0):
    k_required = (2 * E0_natural)**2
    if k_required <= 0: return None
    try:
        log_arg = k_required / C_tunnel
        if log_arg <= np.finfo(float).tiny: return np.inf
        c1 = -alpha_fine_structure * np.log(log_arg)
        return c1 if c1 > 0 else None
    except Exception: return None

C1_REQUIRED = {"electron": get_required_c1(TARGET_ENERGY_ELECTRON),
               "muon": get_required_c1(TARGET_ENERGY_MUON),
               "tau": get_required_c1(TARGET_ENERGY_TAU)}
C1_REQUIRED = {k: v for k, v in C1_REQUIRED.items() if v is not None}

print("\n--- Required c1 Values (if C_tunnel=1.0) ---")
if C1_REQUIRED: [print(f"  {p.capitalize():<9}: {c1:.3f}") for p, c1 in C1_REQUIRED.items()]
else: print("  Error calculating required c1 values.")
print("-------------------------------------------\n")

# --- Stability Criteria (Relaxed Further) ---
N_STRANDS = 3
MAX_BRAID_LENGTH = 7 # Keep slightly longer length

# Filter 1: Non-trivial
def filter_non_trivial(link_obj):
    try: return link_obj.identify() != [(0, 1)] # Check against Unknot ID
    except: return False

# Filter 2: Topological Stability (Relaxed: Allow simple links)
def filter_stable_candidate(link_obj):
    # Allow prime knots OR simple links (e.g., 2 components, low crossings)
    try:
        ident = link_obj.identify()
        num_components = len(ident)
        crossing_num = link_obj.crossing_number()

        if num_components == 1 and ident != [(0,1)]: # Prime Knot
            return True
        elif num_components == 2 and crossing_num <= 4: # Allow simple 2-component links (e.g., Hopf L2a1, Solomon L4a1)
            return True
        # Add criteria for other simple links if desired
        return False
    except Exception as e:
        # print(f"Warning: Stability check failed: {e}")
        return False

# Filter 3: Spinor Compatibility (Assumed possible for j=1/2 structures passing others)
# Filter 4: Mass/k>0 (Filtering based on c1(topo) removed)

# --- Braid Generation (Remains same) ---
def generate_braid_words(n_strands, max_length):
    # ... (function code from previous version) ...
    if n_strands <= 1: yield from []; return
    generators = [f"s{i+1}" for i in range(n_strands - 1)]
    inv_generators = [f"s{i+1}^-1" for i in range(n_strands - 1)]
    all_gens = generators + inv_generators
    queue = [""]
    for length in range(1, max_length + 1):
        next_queue = []; count = 0
        # print(f"INFO: Generating braids of length {length}...") # Verbose
        for word in queue:
            for gen in all_gens:
                last_gen = word.split('*')[-1] if word else None
                if last_gen and \
                   gen.startswith(last_gen.split('^')[0]) and \
                   gen.endswith("^-1") != last_gen.endswith("^-1"): continue
                new_word = f"{word}*{gen}" if word else gen
                yield new_word
                next_queue.append(new_word); count +=1
        queue = next_queue
        # print(f"      Generated {count} braids at length {length}.") # Verbose
        if not queue: break

# --- Main Filter Execution ---
if __name__ == "__main__":
    print("--- Running Relaxed Stability Filter (Including Simple Links) ---")

    start_time = time.time()
    candidate_structures = {} # Store valid candidates: link_name -> data
    processed_links = set(['0_1']) # Track links we found stable candidates for

    braid_generator = generate_braid_words(N_STRANDS, MAX_BRAID_LENGTH)
    count_processed = 0

    print("\nFiltering candidates...")
    for word in braid_generator:
        count_processed += 1
        if count_processed % 1000 == 0: print(f"  Processed {count_processed} braids...")
        try:
            braid_obj = spherogram.Braid(word)
            link_obj = braid_obj.closing()

            # --- Calculate Topology ---
            # Identify might be slow, do basic checks first if possible
            num_components = len(link_obj.components) # Faster check?
            if num_components == 1 and link_obj.is_trivial(): # Fast unknot check
                 knot_name_str = "Unknot"
            else:
                # Use identify for proper classification
                knot_name_tuple = tuple(link_obj.identify())
                knot_name_str = spherogram.knot_db.name_from_tuple(knot_name_tuple) if knot_name_tuple else f"Link_{len(knot_name_tuple)}comp"

            # Skip if already processed this link type
            if knot_name_str in processed_links: continue

            # --- Apply Filters ---
            # Filter 1: Non-trivial
            if knot_name_str == "Unknot": continue

            # Filter 2: Stable Candidate (Prime knot OR Simple Link)
            # Simplify before final stability check if needed
            link_obj_simplified = link_obj.simplified()
            if not filter_stable_candidate(link_obj_simplified): continue

            # Recalculate properties from simplified object
            crossing_num = link_obj_simplified.crossing_number()
            signature_val = link_obj_simplified.signature() # Example invariant
            num_components = len(link_obj_simplified.identify()) # Use identified components

            # Store the first valid candidate found for this link type
            # Use crossing number as primary sort key, then maybe signature for tie-breaking
            sort_key = (crossing_num, abs(signature_val))

            if knot_name_str not in candidate_structures or sort_key < candidate_structures[knot_name_str]['sort_key']:
                 candidate_structures[knot_name_str] = {
                     'complexity_Nc': crossing_num,
                     'signature': signature_val,
                     'components': num_components,
                     'origin_braid_word': word,
                     'sort_key': sort_key
                 }
                 processed_links.add(knot_name_str)
                 print(f"  Found Candidate: {knot_name_str} (Nc={crossing_num}, Sig={signature_val}, Comp={num_components}) from braid '{word}'")

        except Exception as e:
            # Catch errors from braid processing or identification
            # print(f"  Skipping braid '{word}' due to error: {e}")
            pass

    end_time = time.time()
    print(f"\nFinished processing {count_processed} braids in {end_time - start_time:.2f} seconds.")

    # --- Analyze Survivors ---
    print("\n--- Analysis of Filtered Candidate Structures (Prime Knots & Simple Links) ---")
    if not candidate_structures:
        print("No potentially stable candidates found matching relaxed criteria.")
    else:
        # Sort by complexity (crossing number, then |signature|)
        sorted_candidates = sorted(candidate_structures.items(), key=lambda item: item[1]['sort_key'])
        print("Potentially stable knot/link types found (sorted by Nc, |Sig|):")
        for link_name, data in sorted_candidates:
            print(f"  - {link_name}: Nc={data['complexity_Nc']}, Signature={data['signature']}, Components={data['components']} (Origin: {data['origin_braid_word']})")

        # --- Check Lepton Assignment Consistency ---
        print("\n--- Checking Assignment to Lepton Hierarchy based on Complexity (Nc) ---")
        target_c1_list = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)

        # Select only KNOTS (1 component) from candidates for lepton assignment
        knot_candidates = [(name, data) for name, data in sorted_candidates if data['components'] == 1]

        consistent = True
        if len(knot_candidates) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct stable KNOT types found ({len(knot_candidates)}) to map to {len(target_c1_list)} lepton generations.")
             consistent = False
        else:
            print("Hypothetical Assignment (Increasing Knot Nc -> Increasing Mass / Decreasing c1):")
            assigned_complexities = []
            for i in range(len(target_c1_list)):
                 lepton, req_c1 = target_c1_list[i]
                 knot_name, data = knot_candidates[i]
                 assigned_complexities.append(data['complexity_Nc'])
                 print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={data['complexity_Nc']}), Needs c1≈{req_c1:.3f}")

            # Check if complexity strictly increases
            if not all(assigned_complexities[i] < assigned_complexities[i+1] for i in range(len(assigned_complexities)-1)):
                 print("  >> Warning: Complexity (Nc) of simplest knots does not strictly increase.")
            else:
                 print("  >> Observation: Complexity (Nc) of simplest knots found increases correctly.")

            # We cannot check c1 match without a function c1(Topology)
            print("  >> Assessment: Found stable knot candidates whose complexity ordering is *potentially*")
            print("     consistent with lepton hierarchy. Derivation of c1(Topology) needed.")


    print("\n--- Stability Filter Finished ---")