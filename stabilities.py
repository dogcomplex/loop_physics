# --- LSC Stability Filter (v9 - Deep Search for Specific Simple Knots) ---
# --- Searching for Candidate Topologies (Knots & Links) for Leptons ---

import numpy as np
import math
import sys
import time
import signal
import collections
import re # For parsing braid strings if needed later

# Attempt import
try:
    import spherogram
    print("INFO: Spherogram library found.")
    # SnapPy is essential for knot identification in spherogram
    import snappy
    print("INFO: SnapPy library found (needed for knot identification).")
except ImportError as e:
    print(f"ERROR: Required library not found. {e}")
    print("Please install: pip install spherogram snappy")
    print("Note: SnapPy may require additional system dependencies.")
    sys.exit(1)
except RuntimeError as e: print(f"Warning: Spherogram backend issue? {e}")


print("\n--- Loading LSC Deep Stability Filter ---")
print("--- Searching for simplest braid representations of low-crossing prime knots ---")

# --- Physical Constants & Targets ---
alpha_fine_structure = 1 / 137.035999
planck_mass_kg = 2.176434e-8
TARGET_ENERGY_ELECTRON = 4.18540e-23
TARGET_ENERGY_MUON = 8.65429e-21
TARGET_ENERGY_TAU = 1.45532e-19
ASSUMED_C_TUNNEL = 1.0 # Prefactor for k interpretation
def get_required_c1(E0_natural, C_tunnel=1.0):
    """Calculates the c1 needed in k ~ C*exp(-c1/alpha) to get target energy E0."""
    k_required = (2 * E0_natural)**2 # Assumes m_eff=1
    if k_required <= 0: return None
    try:
        log_arg = k_required / C_tunnel
        if log_arg <= np.finfo(float).tiny: return np.inf # Avoid log(0)
        c1 = -alpha_fine_structure * np.log(log_arg)
        return c1 if c1 > 0 else None # c1 should be positive for exp suppression
    except Exception: return None

C1_REQUIRED = {"electron": get_required_c1(TARGET_ENERGY_ELECTRON),
               "muon": get_required_c1(TARGET_ENERGY_MUON),
               "tau": get_required_c1(TARGET_ENERGY_TAU)}
C1_REQUIRED = {k: v for k, v in C1_REQUIRED.items() if v is not None} # Filter out None

print("\n--- Required c1 Values (if C_tunnel=1.0) ---")
if C1_REQUIRED:
    # Sort by required c1 descending (maps to increasing mass)
    C1_REQUIRED_SORTED = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)
    for p, c1 in C1_REQUIRED_SORTED: print(f"  {p.capitalize():<9}: {c1:.3f}")
else: print("  Error calculating required c1 values.")
print("-------------------------------------------\n")


# --- Search Parameters ---
TARGET_KNOTS_INFO = { # Prime knots up to 7 crossings
    # Name : (Nc, |Signature|, Expected DT Code if needed)
    "3_1": (3, 2), # Trefoil (chiral)
    "4_1": (4, 0), # Figure Eight (achiral)
    "5_1": (5, 4), # Cinquefoil (chiral)
    "5_2": (5, 2), # Three-twist (chiral)
    "6_1": (6, 2), # Stevedore knot (chiral)
    "6_2": (6, 0), # (achiral)
    "6_3": (6, 0), # (achiral)
    "7_1": (7, 6), # Heptafoil (chiral)
    "7_2": (7, 4), # (chiral)
    # ... add more knots if needed ...
}
MAX_CANDIDATE_CROSSINGS = 8 # Only interested in knots up to this complexity for leptons

# Known braid representations (as tuples)
KNOWN_BRAIDS = {
    # 3-strand representations
    "3_1": (1, 1, 1),
    "4_1": (1, -2, 1, -2), # Standard short representation
    "5_1": (1, 1, 1, 1, 1),
    "5_2": (1, -2, 1, -2, 1),
    "6_1": (-1, 2, -1, 2, -1, 2), # One common rep
    "6_2": (1, -2, -2, 1, -2, -2),
    "6_3": (-1, -1, 2, -1, -1, 2),
    "7_1": (1, 1, 1, 1, 1, 1, 1),
    # 4-strand representations (can be shorter)
    "4_1_4s": (1, 2, 1, 2), # Figure eight is same on 4 strands
    "5_2_4s": (1, 2, 1, -3, -2, -1), # Example 4-strand rep
    "6_1_4s": (1, -2, 3, 1, -2, -3), # Example 4-strand rep
    # Note: Finding the *absolute minimal* braid index and length is hard.
    # These are just known examples.
}

MAX_BRAID_LENGTH = 12 # Search depth
N_STRANDS_TO_CHECK = [3, 4] # Check 3 and 4 strands

# Helper function to get SnapPy manifold for identification
def get_knot_manifold(name):
    try: return snappy.Manifold(name) # SnapPy understands common knot names
    except Exception as e:
         print(f"Warning: Failed to create SnapPy manifold for {name}: {e}")
         return None

# Helper to format braid word for printing
def format_braid_word(braid_tuple):
    return ' '.join([f"s{abs(g)}{'^-1' if g < 0 else ''}" for g in braid_tuple])

# --- Braid Generation ---
def generate_braid_words_optimized(n_strands, max_length):
    # ... (generator function remains the same as v8) ...
    if n_strands <= 1: yield from []; return
    generators = list(range(1, n_strands))
    inv_generators = [-g for g in generators]
    all_gens = generators + inv_generators
    queue = collections.deque([(g,) for g in all_gens])
    yield from queue
    current_len = 1
    total_yielded = len(queue)
    while queue and current_len < max_length:
        current_len += 1; level_size = len(queue)
        # print(f"INFO: Generating braids of length {current_len} ({n_strands} strands)...") # Verbose
        start_level_time = time.time(); count_yielded_this_level = 0
        for _ in range(level_size):
            word_tuple = queue.popleft()
            last_gen_val = word_tuple[-1] if word_tuple else None
            for gen_val in all_gens:
                if last_gen_val and gen_val == -last_gen_val: continue
                new_word_tuple = word_tuple + (gen_val,)
                yield new_word_tuple; total_yielded += 1
                queue.append(new_word_tuple); count_yielded_this_level += 1
        end_level_time = time.time()
        rate = count_yielded_this_level / (end_level_time - start_level_time + 1e-9)
        # print(f"      Len {current_len}: Gen {count_yielded_this_level} braids ({rate:.1f}/s). Total: {total_yielded}") # Verbose
        if not queue: break


# --- Timeout Handler ---
class TimeoutException(Exception): pass
def timeout_handler(signum, frame): print("\nERROR: Processing timed out!"); raise TimeoutException

# --- Main Filter Execution Function ---
def find_simplest_braids(n_strands, max_braid_length, target_knots_info, time_limit_seconds):
    print(f"\n--- Running Search ({n_strands} STRANDS, MaxLen={max_braid_length}, Timeout={time_limit_seconds}s) ---")
    has_signal = hasattr(signal, 'SIGALRM'); timer_set = False
    if has_signal:
        try:
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(time_limit_seconds); timer_set = True
        except ValueError: print("Warning: Cannot set SIGALRM timer in this environment (e.g., Windows thread).")

    start_time = time.time()
    # Store best candidate found: knot_name -> {'braid_tuple': tuple, 'length': L}
    best_braid_repr = {}
    count_processed = 0; count_identified = 0; count_candidates_found = 0
    timed_out = False; errors_encountered = 0
    last_print_time = time.time()

    # Pre-load known braids for this strand count
    for name, b_tuple in KNOWN_BRAIDS.items():
        if name in target_knots_info:
            max_gen = max(abs(g) for g in b_tuple) if b_tuple else 0
            req_strands = max_gen + 1
            if req_strands == n_strands:
                 best_braid_repr[name] = {'braid_tuple': b_tuple, 'length': len(b_tuple), 'source': 'known'}
                 print(f"  Added known {name} braid (L={len(b_tuple)}): {format_braid_word(b_tuple)}")

    try:
        target_manifolds = {name: get_knot_manifold(name) for name in target_knots_info}
        target_manifolds = {name: mf for name, mf in target_manifolds.items() if mf is not None} # Filter failed lookups

        braid_generator = generate_braid_words_optimized(n_strands, max_braid_length)
        print("\nSearching generated braids...")

        for braid_tuple in braid_generator:
            count_processed += 1
            # Progress reporting
            current_time = time.time()
            if current_time - last_print_time > 5.0: # Print every 5 seconds
                rate = count_processed / (current_time - start_time + 1e-9)
                print(f"  Processed {count_processed} braids... ({rate:.1f} braids/sec)")
                last_print_time = current_time

            try:
                if not braid_tuple: continue
                # Skip if already found shorter representation for all targets
                if len(best_braid_repr) == len(target_knots_info) and \
                   all(data['length'] < len(braid_tuple) for data in best_braid_repr.values()):
                    print(f"INFO: Found shorter representations for all targets. Stopping search early at length {len(braid_tuple)}.")
                    break

                # Create Spherogram objects (can be slow)
                braid_obj = spherogram.Braid(braid_tuple)
                link_obj = braid_obj.closing()

                # --- Quick Filters ---
                if len(link_obj.components) != 1: continue # Only knots
                if link_obj.crossing_number() > MAX_CANDIDATE_CROSSINGS: continue

                # --- Identify using SnapPy (can be slow/error prone) ---
                manifold = None
                try:
                    manifold = link_obj.exterior()
                    # Optional: Check volume/invariants before expensive isometry check?
                    # if manifold.volume() < 0.1: continue # Skip small volume manifolds (likely unknot)
                except Exception as e_mf:
                    # print(f"Warning: Could not get manifold for {format_braid_word(braid_tuple)}: {e_mf}")
                    errors_encountered += 1
                    continue # Skip if manifold creation fails

                # --- Check against Targets ---
                for knot_name, target_manifold in target_manifolds.items():
                    # Skip if already found a shorter or equal length representation
                    if knot_name in best_braid_repr and best_braid_repr[knot_name]['length'] <= len(braid_tuple):
                        continue

                    try:
                        if target_manifold and manifold.is_isometric_to(target_manifold):
                            current_length = len(braid_tuple)
                            best_braid_repr[knot_name] = {
                                'braid_tuple': braid_tuple,
                                'length': current_length,
                                'source': 'search'
                            }
                            count_candidates_found += 1
                            print(f"  ** Found {knot_name}: Length={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                            # Don't break here, maybe a different target matches too? (Unlikely for knots)

                    except snappy.UnsuitableCaseError: pass # is_isometric_to can fail
                    except Exception as e_iso:
                        # print(f"Warning: Isometry check failed for {knot_name} vs braid {braid_tuple}: {e_iso}")
                        errors_encountered += 1
                        pass

            except Exception as e_inner:
                 # print(f"  Skipping braid {format_braid_word(braid_tuple)} due to error: {type(e_inner).__name__}")
                 errors_encountered += 1; pass

    except TimeoutException: timed_out = True
    finally:
        if has_signal and timer_set: signal.alarm(0) # Disable alarm

    end_time = time.time()
    rate = count_processed/(end_time - start_time + 1e-9)
    print(f"\nFinished search ({n_strands} strands, up to maxlen={max_braid_length}) in {end_time - start_time:.2f} seconds.")
    print(f"Processed {count_processed} braids at {rate:.1f} braids/second.")
    if timed_out: print(f"Warning: Process timed out after {time_limit_seconds} seconds, results may be incomplete.")
    if errors_encountered > 0: print(f"Warning: Encountered {errors_encountered} errors during processing.")

    # --- Analyze Survivors ---
    print(f"\n--- Analysis: Simplest Braid Reps Found ({n_strands} Strands) ---")
    if not best_braid_repr: print("No target knots found.")
    else:
        sorted_found_knots = sorted(best_braid_repr.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        print("Simplest braid representations found/known for target knots:")
        for knot_name in sorted_found_knots:
             data = best_braid_repr[knot_name]
             nc, sig = TARGET_KNOTS_INFO[knot_name]
             source = data.get('source', 'search')
             print(f"  - {knot_name:<4} (Nc={nc}, Sig={sig}): Min Braid Length={data['length']} [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")

        # Lepton Assignment Check
        print("\n--- Checking Lepton Assignment Consistency ---")
        target_c1_list = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)
        # Use only found knots for assignment check
        found_knot_list_sorted = [(name, TARGET_KNOTS_INFO[name][0]) for name in sorted_found_knots]

        consistent = True
        if len(found_knot_list_sorted) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct target KNOT types found ({len(found_knot_list_sorted)}) to map to {len(target_c1_list)} lepton generations.")
             consistent = False
        else:
            print("Hypothetical Assignment (Increasing Knot Nc -> Increasing Mass / Decreasing c1):")
            complexities = []
            for i in range(len(target_c1_list)):
                 lepton, req_c1 = target_c1_list[i]
                 # Assign based on complexity order of *found* knots
                 knot_name, nc = found_knot_list_sorted[i]
                 complexities.append(nc)
                 print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")
            # Check if Nc is non-decreasing (allow ties)
            if not all(complexities[i] <= complexities[i+1] for i in range(len(complexities)-1)):
                 print("  >> Warning: Complexity (Nc) ordering inconsistent with lepton hierarchy.")
            else: print("  >> Observation: Complexity (Nc) ordering is potentially consistent.")
            print("  >> Assessment: Requires derivation of c1(Topology) using these candidates.")

    return best_braid_repr


if __name__ == "__main__":
    # Set a long timeout
    timeout_very_long = 28800 # seconds = 8 hours

    all_results = {}

    for n_strands in N_STRANDS_TO_CHECK:
        max_len = MAX_BRAID_LENGTH if n_strands==3 else max(5, MAX_BRAID_LENGTH - 3) # Adjust max_len based on strands
        print("\n" + "="*20 + f" {n_strands}-STRAND DEEP SEARCH (MaxLen={max_len}, Timeout={timeout_very_long}s per run) " + "="*20)
        results = find_simplest_braids(n_strands=n_strands, max_braid_length=max_len,
                                             target_knots_info=TARGET_KNOTS_INFO, time_limit_seconds=timeout_very_long)
        all_results[n_strands] = results

    print("\n" + "="*40)
    print("--- FINAL SUMMARY ACROSS ALL STRANDS ---")
    print("Simplest representations found:")
    combined_best = {}
    for n_strands, results in all_results.items():
         for name, data in results.items():
             if name not in combined_best or data['length'] < combined_best[name]['length']:
                 combined_best[name] = data
                 combined_best[name]['strands'] = n_strands # Record which strand count was best

    if not combined_best: print("No target knots found in any run.")
    else:
        sorted_final = sorted(combined_best.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        for name in sorted_final:
            data = combined_best[name]
            nc, sig = TARGET_KNOTS_INFO[name]
            source = data.get('source', 'search')
            print(f"  - {name:<4} (Nc={nc}): Min Braid Length={data['length']} ({data['strands']} strands) [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")

    print("\n>>> Does the complexity ordering of found knots (e.g., 3_1, 4_1, 5_1/5_2) match lepton hierarchy?")
    print(">>> Final assessment requires theoretical derivation of c1(Topology).")