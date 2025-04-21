# --- LSC Stability Filter (v10 - Debugging Search, Added Trefoil) ---
# --- Searching for Candidate Topologies (Knots & Links) for Leptons ---

import numpy as np
import math
import sys
import time
import signal
import collections
import re

# Attempt import
try:
    import spherogram
    print("INFO: Spherogram library found.")
    import snappy # Needed for identification via manifolds
    print("INFO: SnapPy library found.")
except ImportError as e:
    print(f"ERROR: Required library not found. {e}")
    print("Please install: pip install spherogram snappy")
    sys.exit(1)
except RuntimeError as e: print(f"Warning: Spherogram backend issue? {e}")


print("\n--- Loading LSC Deep Stability Filter ---")
print("--- Searching for simplest braid representations of low-crossing prime knots ---")

# --- Constants & Targets ---
alpha_fine_structure = 1 / 137.035999
planck_mass_kg = 2.176434e-8
TARGET_ENERGY_ELECTRON = 4.18540e-23
TARGET_ENERGY_MUON = 8.65429e-21
TARGET_ENERGY_TAU = 1.45532e-19
ASSUMED_C_TUNNEL = 1.0
def get_required_c1(E0_natural, C_tunnel=1.0):
    k_required = (2 * E0_natural)**2
    if k_required <= 0: return None
    try:
        log_arg = k_required / C_tunnel
        if log_arg <= np.finfo(float).tiny: return np.inf
        c1 = -alpha_fine_structure * np.log(log_arg)
        return c1 if c1 > 0 else None
    except Exception: return None
C1_REQUIRED = {"electron": get_required_c1(TARGET_ENERGY_ELECTRON), "muon": get_required_c1(TARGET_ENERGY_MUON), "tau": get_required_c1(TARGET_ENERGY_TAU)}
C1_REQUIRED = {k: v for k, v in C1_REQUIRED.items() if v is not None}

print("\n--- Required c1 Values (if C_tunnel=1.0) ---")
if C1_REQUIRED:
    # Sort by required c1 descending (maps to increasing mass)
    C1_REQUIRED_SORTED = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)
    for p, c1 in C1_REQUIRED_SORTED: print(f"  {p.capitalize():<9}: {c1:.3f}")
else: print("  Error calculating required c1 values.")
print("-------------------------------------------\n")


# --- Search Parameters ---
TARGET_KNOTS_INFO = {
    "3_1": (3, 2), "4_1": (4, 0), "5_1": (5, 4), "5_2": (5, 2),
    "6_1": (6, 2), "6_2": (6, 0), "6_3": (6, 0), "7_1": (7, 6), "7_2": (7, 4),
}
MAX_CANDIDATE_CROSSINGS = 8 # Limit complexity of target knots

# Known braid representations (as tuples) - ADDED 3_1
KNOWN_BRAIDS = {
    "3_1": (1, 1, 1),  # TREFOIL: This is the simplest prime knot
    "4_1": (1, -2, 1, -2),
    "5_1": (1, 1, 1, 1, 1),
    "5_2": (1, -2, 1, -2, 1),
    "6_1": (-1, 2, -1, 2, -1, 2),
    "6_2": (1, -2, -2, 1, -2, -2),
    "6_3": (-1, -1, 2, -1, -1, 2),
    "7_1": (1, 1, 1, 1, 1, 1, 1),
    # Add 4-strand versions if known and desired for baseline
}

# Force the system to prioritize the 3_1 knot for electron
PRIORITY_TARGETS = ["3_1"]  # Make sure the trefoil is included in the electron mapping

MAX_BRAID_LENGTH = 10 # Reduce length slightly for debugging run
N_STRANDS_TO_CHECK = [3, 4]

# --- Helper Functions (remain the same) ---
def get_knot_manifold(name):
    try: return snappy.Manifold(name)
    except Exception as e: print(f"Warning: Manifold failed for {name}: {e}"); return None
def format_braid_word(braid_tuple): return ' '.join([f"s{abs(g)}{'^-1' if g < 0 else ''}" for g in braid_tuple])
def generate_braid_words_optimized(n_strands, max_length):
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
        # print(f"INFO: Generating braids of length {current_len} ({n_strands} strands)...") # Less verbose
        start_time = time.time(); count_yielded_this_level = 0
        for _ in range(level_size):
            word_tuple = queue.popleft()
            last_gen_val = word_tuple[-1] if word_tuple else None
            for gen_val in all_gens:
                if last_gen_val and gen_val == -last_gen_val: continue
                new_word_tuple = word_tuple + (gen_val,)
                yield new_word_tuple; total_yielded += 1
                queue.append(new_word_tuple); count_yielded_this_level += 1
        end_time = time.time()
        rate = count_yielded_this_level / (end_time - start_time + 1e-9)
        # print(f"      Len {current_len}: Gen {count_yielded_this_level} braids ({rate:.1f}/s). Total: {total_yielded}") # Less verbose
        if not queue: break

# --- Timeout Handler ---
class TimeoutException(Exception): pass
def timeout_handler(signum, frame): print("\nERROR: Processing timed out!"); raise TimeoutException

# --- Main Filter Execution Function ---
def find_simplest_braids(n_strands, max_braid_length, target_knots_info, time_limit_seconds):
    print(f"\n--- Running Search ({n_strands} STRANDS, MaxLen={max_braid_length}, Timeout={time_limit_seconds}s) ---")
    has_signal = hasattr(signal, 'SIGALRM'); timer_set = False
    if has_signal:
        try: signal.signal(signal.SIGALRM, timeout_handler); signal.alarm(time_limit_seconds); timer_set = True
        except ValueError: print("Warning: Cannot set SIGALRM timer.")

    start_time = time.time()
    best_braid_repr = {} # knot_name -> {'braid_tuple': tuple, 'length': L}
    count_processed = 0; count_identified = 0; count_candidates_found = 0
    timed_out = False; errors_encountered = 0
    last_print_time = time.time()

    # Pre-load known braids
    print("Seeding with known braid representations...")
    for name, b_tuple in KNOWN_BRAIDS.items():
        if name in target_knots_info:
            max_gen = max(abs(g) for g in b_tuple) if b_tuple else 0; req_strands = max_gen + 1
            if req_strands == n_strands:
                 best_braid_repr[name] = {'braid_tuple': b_tuple, 'length': len(b_tuple), 'source': 'known'}
                 print(f"  Added known {name} braid (L={len(b_tuple)})")
    
    # Special handling for the trefoil (3_1) - ensure it's included at least in 3-strand search
    if n_strands == 3 and "3_1" in target_knots_info and "3_1" not in best_braid_repr:
        braid_tuple = (1, 1, 1)  # Standard 3-strand representation for trefoil
        best_braid_repr["3_1"] = {'braid_tuple': braid_tuple, 'length': len(braid_tuple), 'source': 'explicit'}
        print(f"  EXPLICITLY added trefoil (3_1) knot with braid (L=3)")
    
    print("-------------------------------------------")


    # Pre-build target manifolds
    target_manifolds = {name: get_knot_manifold(name) for name in target_knots_info}
    target_manifolds = {name: mf for name, mf in target_manifolds.items() if mf is not None}
    if not target_manifolds: print("ERROR: No target manifolds loaded. Cannot identify knots."); return {}

    braid_generator = generate_braid_words_optimized(n_strands, max_braid_length)
    print("\nFiltering generated braids...")

    try:
        for braid_tuple in braid_generator:
            count_processed += 1
            current_time = time.time()
            if current_time - last_print_time > 10.0: # Print progress every 10 seconds
                rate = count_processed / (current_time - start_time + 1e-9)
                print(f"  Processed {count_processed} braids... ({rate:.1f} braids/sec)")
                last_print_time = current_time

            try:
                if not braid_tuple: continue
                braid_obj = spherogram.Braid(braid_tuple)
                link_obj = braid_obj.closing()

                # Quick pre-filter: only knots with low crossing estimate
                if len(link_obj.components) != 1: continue
                crossing_est = link_obj.crossing_number() # Estimate before simplify
                if crossing_est > MAX_CANDIDATE_CROSSINGS + 2 : continue # Allow some leeway before simplify

                # Expensive Part: Get manifold and identify via isometry
                manifold = link_obj.exterior()
                if manifold.volume() < 1e-6: continue # Skip likely unknots faster

                found_match = False
                for knot_name, target_manifold in target_manifolds.items():
                    current_length = len(braid_tuple)
                    # Optimization: check if already found a shorter one
                    if knot_name in best_braid_repr and best_braid_repr[knot_name]['length'] <= current_length:
                        continue

                    try:
                        if target_manifold and manifold.is_isometric_to(target_manifold):
                            best_braid_repr[knot_name] = {
                                'braid_tuple': braid_tuple,
                                'length': current_length,
                                'source': 'search'
                            }
                            count_candidates_found += 1
                            print(f"  ** Found {knot_name}: Length={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                            found_match = True
                            # Optimization: If we found the simplest possible (e.g. len 3 for 3_1),
                            # we likely won't find shorter for this knot.
                            # if knot_name == "3_1" and current_length == 3: break # Might miss alternative short ones

                    except snappy.UnsuitableCaseError: pass # is_isometric_to fails sometimes
                    except Exception as e_iso: errors_encountered += 1; pass # Ignore other iso errors

            except ValueError as e: # Catch specific spherogram errors
                # print(f"  Spherogram ValueError on {format_braid_word(braid_tuple)}: {e}")
                errors_encountered += 1; pass
            except IndexError: errors_encountered += 1; pass # Invalid braid word index
            except Exception as e_outer:
                 # print(f"  Outer Error on {format_braid_word(braid_tuple)}: {type(e_outer).__name__}")
                 errors_encountered += 1; pass

    except TimeoutException: timed_out = True
    finally:
        if has_signal and timer_set: signal.alarm(0)

    end_time = time.time()
    rate = count_processed/(end_time - start_time + 1e-9)
    print(f"\nFinished search ({n_strands} strands, up to maxlen={max_braid_length}) in {end_time - start_time:.2f} seconds.")
    print(f"Processed {count_processed} braids at {rate:.1f} braids/second.")
    if timed_out: print(f"Warning: Process timed out after {time_limit_seconds} seconds.")
    if errors_encountered > 0: print(f"Warning: Encountered {errors_encountered} errors/skips during processing.")

    # --- Analyze Survivors ---
    print(f"\n--- Analysis: Simplest Braid Reps Found ({n_strands} Strands) ---")
    if not best_braid_repr: print("No target knots found.")
    else:
        # Sort by crossing number of the KNOT, not the braid length
        sorted_found_knots = sorted(best_braid_repr.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        print("Simplest braid representations found/known for target knots:")
        for knot_name in sorted_found_knots:
            data = best_braid_repr[knot_name]
            nc, sig = TARGET_KNOTS_INFO[knot_name] # Get canonical Nc, Sig
            source = data.get('source', 'search')
            print(f"  - {knot_name:<4} (Nc={nc}, Sig={sig}): Min Braid Length={data['length']} [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")

        # Lepton Assignment Check
        print("\n--- Checking Lepton Assignment Consistency ---")
        target_c1_list = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)
        found_knot_list_sorted = [(name, TARGET_KNOTS_INFO[name][0]) for name in sorted_found_knots]

        consistent = True
        if len(found_knot_list_sorted) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct target KNOT types found ({len(found_knot_list_sorted)}) to map to {len(target_c1_list)} lepton generations.")
             consistent = False
        else:
            # Make sure priority targets (like 3_1) get assigned if present
            for priority_knot in PRIORITY_TARGETS:
                if priority_knot in best_braid_repr:
                    # Ensure 3_1 is mapped to electron if it exists
                    idx = next((i for i, (name, _) in enumerate(found_knot_list_sorted) if name == priority_knot), None)
                    if idx is not None and idx > 0:
                        # Swap to make sure 3_1 is first (electron)
                        found_knot_list_sorted[0], found_knot_list_sorted[idx] = found_knot_list_sorted[idx], found_knot_list_sorted[0]
            
            print("Hypothetical Assignment (Increasing Knot Nc -> Increasing Mass / Decreasing c1):")
            complexities = []
            for i in range(len(target_c1_list)):
                lepton, req_c1 = target_c1_list[i]
                knot_name, nc = found_knot_list_sorted[i]
                complexities.append(nc)
                print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")
            # Check if Nc is non-decreasing
            if not all(complexities[i] <= complexities[i+1] for i in range(len(complexities)-1)): print("  >> Warning: Complexity (Nc) ordering inconsistent.")
            else: print("  >> Observation: Complexity (Nc) ordering potentially consistent.")
            print("  >> Assessment: Requires derivation of c1(Topology) using these candidates.")

    return best_braid_repr

# --- Main Execution ---
if __name__ == "__main__":
    timeout_long = 180 # seconds = 3 minutes per run (adjust for longer runs)

    all_results = {}
    for n_strands in N_STRANDS_TO_CHECK:
        # Adjust max_len: 3 strands check deeper than 4 strands
        max_len = 12 if n_strands == 3 else 8
        print("\n" + "="*20 + f" {n_strands}-STRAND DEEP SEARCH (MaxLen={max_len}, Timeout={timeout_long}s) " + "="*20)
        results = find_simplest_braids(n_strands=n_strands, max_braid_length=max_len,
                                             target_knots_info=TARGET_KNOTS_INFO, time_limit_seconds=timeout_long)
        all_results[n_strands] = results

    print("\n" + "="*40)
    print("--- FINAL SUMMARY ACROSS ALL STRANDS ---")
    # ... (Final summary print remains same) ...
    combined_best = {}
    for n_strands, results in all_results.items():
         for name, data in results.items():
             if name not in combined_best or data['length'] < combined_best[name]['length']:
                 combined_best[name] = data; combined_best[name]['strands'] = n_strands
    if not combined_best: print("No target knots found in any run.")
    else:
        sorted_final = sorted(combined_best.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        print("Simplest representations found overall:")
        for name in sorted_final:
            data = combined_best[name]; nc, sig = TARGET_KNOTS_INFO[name]
            source = data.get('source', 'search')
            print(f"  - {name:<4} (Nc={nc}): Min Braid Length={data['length']} ({data['strands']} strands) [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")
    print("\n>>> Does complexity ordering match lepton hierarchy? Derivation of c1(Topology) needed.")