# --- LSC Stability Filter (v10 - Explicit Trefoil, Robust Search) ---
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
    import snappy # Needed for identification via manifolds
    print("INFO: SnapPy library found.")
except ImportError as e:
    print(f"ERROR: Required library not found. {e}")
    print("Please install: pip install spherogram snappy")
    print("Note: SnapPy may require additional system dependencies.")
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
TARGET_KNOTS_INFO = { # Prime knots up to 7 crossings
    # Name : (Nc, |Signature|)
    "3_1": (3, 2), # Trefoil (chiral)
    "4_1": (4, 0), # Figure Eight (achiral)
    "5_1": (5, 4), # Cinquefoil (chiral)
    "5_2": (5, 2), # Three-twist (chiral)
    "6_1": (6, 2), # Stevedore knot (chiral)
    "6_2": (6, 0), # (achiral)
    "6_3": (6, 0), # (achiral)
    "7_1": (7, 6), # Heptafoil (chiral)
    "7_2": (7, 4), # (chiral)
}
MAX_CANDIDATE_CROSSINGS = 8 # Limit complexity of target knots

# Known braid representations (as tuples) - INCLUDING 3_1
KNOWN_BRAIDS = {
    "3_1": (1, 1, 1), # Explicitly add simplest trefoil rep
    "4_1": (1, -2, 1, -2),
    "5_1": (1, 1, 1, 1, 1),
    "5_2": (1, -2, 1, -2, 1),
    "6_1": (-1, 2, -1, 2, -1, 2),
    "6_2": (1, -2, -2, 1, -2, -2),
    "6_3": (-1, -1, 2, -1, -1, 2),
    "7_1": (1, 1, 1, 1, 1, 1, 1),
}

MAX_BRAID_LENGTH = 12 # Keep reasonably deep search
N_STRANDS_TO_CHECK = [2, 3, 4]

# --- Helper Functions ---
def get_knot_manifold(name):
    try: return snappy.Manifold(name)
    except Exception as e: print(f"Warning: Manifold failed for {name}: {e}"); return None
def format_braid_word(braid_tuple): return ' '.join([f"s{abs(g)}{'^-1' if g < 0 else ''}" for g in braid_tuple])
def generate_braid_words_optimized(n_strands, max_length):
    # ... (generator code unchanged) ...
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
    best_braid_repr = {} # knot_name -> {'braid_tuple': tuple, 'length': L, 'source': 'known'/'search'}
    count_processed = 0; count_candidates_found = 0
    timed_out = False; errors_encountered = 0
    last_print_time = time.time()
    error_types = collections.Counter()

    # Pre-load known braids for this strand count
    print("Seeding with known braid representations...")
    known_added_count = 0
    for name, b_tuple in KNOWN_BRAIDS.items():
        if name in target_knots_info:
            max_gen = max(abs(g) for g in b_tuple) if b_tuple else 0; req_strands = max_gen + 1
            if req_strands == n_strands:
                 best_braid_repr[name] = {'braid_tuple': b_tuple, 'length': len(b_tuple), 'source': 'known'}
                 print(f"  Added known {name} braid (L={len(b_tuple)})")
                 known_added_count += 1
    print(f"Seeded {known_added_count} known braids for {n_strands} strands.")
    print("-------------------------------------------")

    # Pre-build target manifolds for efficiency
    target_manifolds = {name: get_knot_manifold(name) for name in target_knots_info}
    target_manifolds = {name: mf for name, mf in target_manifolds.items() if mf is not None}
    if not target_manifolds: print("ERROR: No target manifolds loaded. Cannot identify knots."); return {}

    # Debug: Print loaded manifolds
    print(f"Loaded {len(target_manifolds)} target manifolds: {', '.join(target_manifolds.keys())}")

    braid_generator = generate_braid_words_optimized(n_strands, max_braid_length)
    print("\nFiltering generated braids...")

    try:
        for braid_tuple in braid_generator:
            count_processed += 1
            current_time = time.time()
            if current_time - last_print_time > 10.0: # Progress update
                rate = count_processed / (current_time - start_time + 1e-9)
                print(f"  Processed {count_processed} braids... ({rate:.1f} braids/sec)")
                print(f"  Error types so far: {dict(error_types)}")
                last_print_time = current_time

            try:
                if not braid_tuple: continue

                # Optimization: Check if this braid is longer than shortest found for all targets
                if len(best_braid_repr) == len(target_knots_info) and \
                   all(data['length'] <= len(braid_tuple) for data in best_braid_repr.values()):
                   continue # Skip if cannot improve any result

                try:
                    # Skip optimization patterns based on braid structure
                    # 1. Braids that sum to 0 are often trivial or equivalent to shorter ones
                    if sum(braid_tuple) == 0:
                        continue
                        
                    # 2. Skip repetitive patterns like (1,1,-1,-1) or (1,-1,1,-1)
                    if len(braid_tuple) >= 4:
                        # Pattern: generator followed by its inverse
                        has_cancellation = False
                        for i in range(len(braid_tuple)-1):
                            if braid_tuple[i] == -braid_tuple[i+1]:
                                has_cancellation = True
                                break
                        if has_cancellation:
                            continue
                    
                    # Handle single-element tuples specially (they cause 'int' is not iterable errors)
                    if len(braid_tuple) == 1:
                        # Create a list from the single value for proper handling
                        braid_obj = spherogram.ClosedBraid([braid_tuple[0]])
                    else:
                        braid_obj = spherogram.ClosedBraid(*braid_tuple)
                except Exception as e:
                    error_types['braid_creation'] += 1
                    errors_encountered += 1
                    if count_processed <= 5: print(f"ERROR on Braid creation: {e} for braid {braid_tuple}")
                    continue

                try:
                    link_obj = braid_obj  # ClosedBraid is already a Link object, no closing needed
                except Exception as e:
                    error_types['braid_closing'] += 1
                    errors_encountered += 1
                    if count_processed <= 5: print(f"ERROR on Braid closing: {e}")
                    continue

                # Quick Pre-filters
                if len(link_obj.link_components) != 1: 
                    continue
                
                # Additional pre-filtering to save time on expensive exterior calculation
                try:
                    # Skip if crossing number is clearly too high (adds overhead but saves time overall)
                    if hasattr(link_obj, 'crossing_number'):
                        crossings = link_obj.crossing_number()
                        if crossings > MAX_CANDIDATE_CROSSINGS + 2: 
                            continue
                except Exception:
                    pass  # If crossing_number fails, just proceed
                
                # Expensive Part: Identify via isometry
                try:
                    manifold = link_obj.exterior()
                except Exception as e:
                    error_types['exterior_creation'] += 1
                    errors_encountered += 1
                    if count_processed <= 5: print(f"ERROR on exterior creation: {e}")
                    continue
                
                # Optional Faster Check: if manifold.volume() < 0.1: continue

                found_match_this_braid = False
                for knot_name, target_manifold in target_manifolds.items():
                    current_length = len(braid_tuple)
                    # Check if it's a target and if we found a shorter rep
                    if knot_name in target_knots_info and \
                       (knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']):
                        try:
                            if target_manifold and manifold.is_isometric_to(target_manifold):
                                # Found a match for this target knot!
                                best_braid_repr[knot_name] = {
                                    'braid_tuple': braid_tuple,
                                    'length': current_length,
                                    'source': 'search'
                                }
                                count_candidates_found += 1
                                print(f"  ** Found {knot_name}: New Shortest Length={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                                found_match_this_braid = True
                                # Don't break, could match multiple targets (unlikely for knots)
                        except Exception as e_iso: 
                            # Many errors are expected for non-matching manifolds
                            error_message = str(e_iso)
                            error_types['isometry_check'] += 1
                            errors_encountered += 1
                            
                            # Only print errors for the first few braids or rare error types
                            common_error = "The SnapPea kernel was not able to determine if the manifolds are isometric."
                            if count_processed <= 5 and error_message == common_error:
                                # Print just once for common errors during initial processing
                                if error_types['isometry_check'] == 1:
                                    print(f"INFO: Suppressing common expected error: '{common_error}'")
                            elif error_message != common_error:
                                # Always print uncommon errors, they might be important
                                print(f"UNCOMMON ERROR: {error_message}")
                            pass

            except ValueError as e: 
                error_types['value_error'] += 1
                errors_encountered += 1
                if count_processed <= 5: print(f"ERROR ValueError: {e}")
                pass
            except IndexError as e: 
                error_types['index_error'] += 1
                errors_encountered += 1
                if count_processed <= 5: print(f"ERROR IndexError: {e}")
                pass
            except Exception as e_outer: 
                error_types['other_exception'] += 1
                errors_encountered += 1
                if count_processed <= 5: print(f"ERROR outer exception: {e_outer}")
                pass

    except TimeoutException: timed_out = True
    finally:
        if has_signal and timer_set: signal.alarm(0)

    end_time = time.time(); elapsed_time = end_time - start_time
    rate = count_processed/(elapsed_time + 1e-9)
    print(f"\nFinished search ({n_strands} strands, up to maxlen={max_braid_length}) in {elapsed_time:.2f} seconds.")
    print(f"Processed {count_processed} braids at {rate:.1f} braids/second.")
    if timed_out: print(f"Warning: Process timed out after {time_limit_seconds} seconds.")
    if errors_encountered > 0: 
        print(f"Warning: Encountered {errors_encountered} errors/skips during processing.")
        print(f"Error distribution: {dict(error_types)}")

    # --- Analyze Survivors ---
    print(f"\n--- Analysis: Simplest Braid Reps Found ({n_strands} Strands) ---")
    if not best_braid_repr: print("No target knots found.")
    else:
        # Sort by canonical crossing number from TARGET_KNOTS_INFO
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
        # Assign found knots based on THEIR complexity order
        found_knot_list_sorted_by_Nc = sorted(
            [(name, TARGET_KNOTS_INFO[name][0]) for name in sorted_found_knots],
            key=lambda item: item[1]
        )

        consistent = True
        if len(found_knot_list_sorted_by_Nc) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct target KNOT types found ({len(found_knot_list_sorted_by_Nc)}) to map to {len(target_c1_list)} lepton generations.")
             consistent = False
        else:
            print("Hypothetical Assignment (Increasing Knot Nc -> Increasing Mass / Decreasing c1):")
            complexities = []
            for i in range(len(target_c1_list)):
                 lepton, req_c1 = target_c1_list[i]
                 knot_name, nc = found_knot_list_sorted_by_Nc[i] # Assign by complexity order
                 complexities.append(nc)
                 print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")
            # Check if Nc is non-decreasing
            if not all(complexities[i] <= complexities[i+1] for i in range(len(complexities)-1)): print("  >> Warning: Complexity (Nc) ordering inconsistent.")
            else: print("  >> Observation: Complexity (Nc) ordering potentially consistent.")
            print("  >> Assessment: Requires derivation of c1(Topology) using these candidates.")

    return best_braid_repr


if __name__ == "__main__":
    timeout_long = 180 # Reset to 3 minutes for this corrected run

    all_results = {}
    for n_strands in N_STRANDS_TO_CHECK:
        # Adjust max_len based on strands (e.g., 10 for n=3, 7 for n=4)
        if n_strands == 2:
            max_len = 12  # Longer for 2 strands since there's only one generator
        elif n_strands == 3:
            max_len = 10
        else:
            max_len = 7
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
    print(">>> Next step: Theoretical derivation of c1(Topology) relationship.")