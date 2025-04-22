# --- LSC Stability Filter (v12 - Hybrid Identification) ---
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

# Check if running inside Sage (required for Jones polynomial)
try:
    import sage.all
    print("INFO: Running inside SageMath environment.")
    SAGE_AVAILABLE = True
except ImportError:
    print("WARNING: Not running inside SageMath. Jones polynomial calculation will fail.")
    SAGE_AVAILABLE = False


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

# Define which knots are expected to be non-hyperbolic (volume ~ 0)
# Torus knots (p,q) like 3_1 (3,2), 5_1 (5,2), 7_1 (7,2) are common examples
NON_HYPERBOLIC_TARGETS = {"3_1", "5_1", "7_1"}
HYPERBOLIC_TARGETS = {k for k in TARGET_KNOTS_INFO if k not in NON_HYPERBOLIC_TARGETS}
VOLUME_THRESHOLD = 0.1 # Threshold to distinguish hyperbolic/non-hyperbolic

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
def generate_braid_words_optimized(n_strands, max_length, debug=False):
    if n_strands <= 1: yield from []; return
    generators = list(range(1, n_strands))
    inv_generators = [-g for g in generators]
    all_gens = generators + inv_generators
    
    # --- Debug for trefoil ---
    trefoil_found = False
    # --- End debug for trefoil ---
    
    queue = collections.deque([(g,) for g in all_gens])
    
    # Check if we're debugging 2-strand braids that might contain the trefoil
    if debug and n_strands == 2:
        print("\nDEBUG: Starting braid generation for 2 strands...")
        print(f"DEBUG: Initial generators: {all_gens}")
        print(f"DEBUG: Initial queue: {list(queue)}")
    
    # Process length 1 braids
    for braid in queue:
        if debug and n_strands == 2:
            print(f"DEBUG: Yielding Length 1 braid: {braid}")
        yield braid
    
    current_len = 1
    total_yielded = len(queue)
    
    while queue and current_len < max_length:
        current_len += 1
        level_size = len(queue)
        start_level_time = time.time()
        count_yielded_this_level = 0
        
        if debug and n_strands == 2 and current_len <= 3:
            print(f"\nDEBUG: Processing braids of length {current_len}...")
        
        for _ in range(level_size):
            word_tuple = queue.popleft()
            last_gen_val = word_tuple[-1] if word_tuple else None
            
            for gen_val in all_gens:
                # Skip if adding inverse of last generator (would cancel)
                if last_gen_val and gen_val == -last_gen_val: continue
                
                new_word_tuple = word_tuple + (gen_val,)
                
                # --- Debug check for trefoil braid (1,1,1) ---
                if new_word_tuple == (1,1,1) and not trefoil_found:
                    trefoil_found = True
                    print(f"\nDEBUG: TREFOIL BRAID (1,1,1) GENERATED at position {total_yielded + count_yielded_this_level + 1}")
                # --- End debug check ---
                
                if debug and n_strands == 2 and current_len <= 3:
                    print(f"DEBUG: Yielding Length {current_len} braid: {new_word_tuple}")
                
                yield new_word_tuple
                total_yielded += 1
                queue.append(new_word_tuple)
                count_yielded_this_level += 1
        
        end_level_time = time.time()
        rate = count_yielded_this_level / (end_level_time - start_level_time + 1e-9)
        
        if debug and n_strands == 2:
            print(f"DEBUG: Generated {count_yielded_this_level} braids of length {current_len} ({rate:.1f}/s)")
            print(f"DEBUG: Total braids so far: {total_yielded}")
        
        if not queue: break

# --- Timeout Handler ---
class TimeoutException(Exception): pass
def timeout_handler(signum, frame): print("\nERROR: Processing timed out!"); raise TimeoutException

# --- Function to get Jones polynomial safely (Requires Sage) ---
def get_jones_polynomial(link_obj, braid_tuple, error_state):
    if not SAGE_AVAILABLE:
        return None # Cannot calculate outside Sage
    try:
        poly = link_obj.jones_polynomial()
        return poly # Return the raw Sage polynomial object
    except spherogram.sage_helper.SageNotAvailable as e:
        print(f"      ERROR: Sage not available for Jones polynomial: {e}")
        return None
    except AssertionError as e_assert:
        # Print specific info for AssertionErrors, now using passed tuple
        # Limit the number of times this is printed per run
        if error_state['printed_count'] < error_state['max_prints']:
             print(f"      ASSERTION ERROR calculating Jones polynomial for braid {braid_tuple}: {e_assert}")
             error_state['printed_count'] += 1
             if error_state['printed_count'] == error_state['max_prints'] and not error_state['suppressed_msg_shown']:
                 print("      (Suppressing further polynomial assertion errors for this run)")
                 error_state['suppressed_msg_shown'] = True
        # Always return None on assertion error
        return None # Treat assertion errors as calculation failure
    except Exception as e:
        # Handle potential errors during polynomial calculation
        print(f"      ERROR calculating Jones polynomial: {e} (Type: {type(e)})")
        return None

# --- Main Filter Execution Function ---
def find_simplest_braids(n_strands, max_braid_length, target_knots_info, time_limit_seconds):
    print(f"\n--- Running Search ({n_strands} STRANDS, MaxLen={max_braid_length}, Timeout={time_limit_seconds}s) ---")
    print(f"--- Using HYBRID Identification (Manifold Vol > {VOLUME_THRESHOLD} ? SnapPy : Jones Poly) ---")
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

    # State for managing assertion error printing
    poly_assertion_state = {'printed_count': 0, 'max_prints': 5, 'suppressed_msg_shown': False}

    # Pre-load known braids for this strand count
    print("\nSeeding with known braid representations...")
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

    # --- Pre-calculate Target Manifolds (Hyperbolic) & Polynomials (Non-Hyperbolic) --- 
    target_manifolds = {}
    target_polynomials = {}
    target_links = {} # Store Link objects for non-hyperbolic targets
    manifold_errors = 0
    poly_errors = 0

    print("\nPre-calculating for HYPERBOLIC targets (Manifolds)...")
    for name in HYPERBOLIC_TARGETS:
        if name in target_knots_info:
            mf = get_knot_manifold(name)
            if mf is not None:
                target_manifolds[name] = mf
            else:
                print(f"  WARNING: Failed to get manifold for hyperbolic target {name}.")
                manifold_errors += 1
    print(f"Pre-calculated {len(target_manifolds)} manifolds for {len(HYPERBOLIC_TARGETS)} hyperbolic targets.")

    print("\nPre-calculating for NON-HYPERBOLIC targets (Jones Polynomials)...")
    if not SAGE_AVAILABLE:
        print("  SKIPPING polynomial calculation: Not running in Sage.")
    else:
        # Need dummy braid tuple for pre-calculation call (not used in print on success)
        dummy_tuple = ('target',) 
        for name in NON_HYPERBOLIC_TARGETS:
            if name in target_knots_info:
                try:
                    link = spherogram.Link(name)
                    target_links[name] = link # Store link for mirror check later
                    poly = get_jones_polynomial(link, dummy_tuple, poly_assertion_state) # Pass state
                    if poly is not None:
                        target_polynomials[name] = poly
                        # print(f"  {name}: {poly}") # Optional debug print
                    else:
                        print(f"  WARNING: Failed to calculate Jones polynomial for non-hyperbolic target {name}.")
                        poly_errors += 1
                except Exception as e:
                    print(f"  ERROR: Failed to create Link object for non-hyperbolic target {name}: {e}")
                    poly_errors += 1
        print(f"Pre-calculated {len(target_polynomials)} polynomials for {len(NON_HYPERBOLIC_TARGETS)} non-hyperbolic targets.")

    if manifold_errors > 0:
        print(f"WARNING: Encountered {manifold_errors} errors during target manifold calculation.")
    if poly_errors > 0:
        print(f"WARNING: Encountered {poly_errors} errors during target polynomial calculation.")
    if not target_manifolds and not target_polynomials:
        print("ERROR: No target manifolds or polynomials loaded/calculated. Cannot identify knots.")
        return {}
    print("-------------------------------------------")

    # Enable debug tracking for 2-strand braids to identify trefoil issues
    enable_debug = (n_strands == 2)
    braid_generator = generate_braid_words_optimized(n_strands, max_braid_length, debug=enable_debug)
    print("\nFiltering generated braids using Hybrid Method...")

    # Add specific trefoil debug for 2-strand case
    if n_strands == 2:
        print("DEBUG: Specifically tracking trefoil knot (braid: 1,1,1) processing")

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
                
                # Special debugging for trefoil braid
                is_trefoil_braid = (braid_tuple == (1, 1, 1))
                if is_trefoil_braid and n_strands == 2:
                    # Keep this debug print if helpful for generation tracking
                    print(f"\nDEBUG: GENERATED/PROCESSING TREFOIL BRAID {braid_tuple}")

                # Optimization: Check if this braid is longer than shortest found for all targets
                # Needs to check against combined targets
                all_target_names = HYPERBOLIC_TARGETS.union(NON_HYPERBOLIC_TARGETS)
                if len(best_braid_repr) == len(all_target_names) and \
                   all(data['length'] <= len(braid_tuple) for name, data in best_braid_repr.items() if name in all_target_names):
                    if is_trefoil_braid and n_strands == 2:
                        # Keep this debug print if helpful for generation tracking
                        print(f"\nDEBUG: GENERATED/PROCESSING TREFOIL BRAID {braid_tuple}")
                    else:
                        continue # Skip if cannot improve any result

                try:
                    # Skip optimization patterns based on braid structure
                    # Debug check for trefoil braid
                    is_trefoil_debug = (braid_tuple == (1, 1, 1) and n_strands == 2)
                    
                    # 1. Braids that sum to 0 are often trivial or equivalent to shorter ones
                    if sum(braid_tuple) == 0:
                        if is_trefoil_debug:
                            print(f"DEBUG: Trefoil would be filtered by sum=0 check, but sum={sum(braid_tuple)}")
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
                            if is_trefoil_debug:
                                print("DEBUG: Trefoil would be filtered by cancellation check")
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

                # REMOVED SPECIAL CASE for trefoil (1, 1, 1) - rely on hybrid method now
                # if braid_tuple == (1, 1, 1):
                #     ... removed ...

                try:
                    link_obj = braid_obj  # ClosedBraid is already a Link object
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
                current_length = len(braid_tuple)

                # --- Hybrid Identification --- 
                manifold = None
                volume = -1.0 # Default to invalid volume
                try:
                    # Still need the manifold for volume check
                    manifold = link_obj.exterior()
                    # Calculate volume here (might raise error for some links)
                    volume = manifold.volume()
                except Exception as e:
                    # Errors can happen here for complex/degenerate links
                    error_types['exterior_vol_error'] += 1
                    errors_encountered += 1
                    # Don't print every time, can be noisy
                    # if count_processed <= 10 or errors_encountered % 100 == 0:
                    #    print(f"Warning: Exterior/Volume calculation failed: {e}")
                    continue # Cannot identify if exterior/volume fails

                # --- Branch based on Volume --- 
                if volume > VOLUME_THRESHOLD:
                    # --- Hyperbolic Path: Use Snappy Isometry --- 
                    identification_method = "Isometry" 
                    if not target_manifolds:
                         # error_types['no_hyperbolic_targets'] += 1 # Redundant?
                         pass # No hyperbolic targets to check against
                    else:
                        for knot_name, target_manifold in target_manifolds.items():
                            # Check if it's a potential target and shorter rep needed
                            if knot_name in target_knots_info and \
                               (knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']):
                                try:
                                    if manifold.is_isometric_to(target_manifold):
                                        # Found a match!
                                        best_braid_repr[knot_name] = {
                                            'braid_tuple': braid_tuple,
                                            'length': current_length,
                                            'source': 'search_iso'
                                        }
                                        count_candidates_found += 1
                                        print(f"  ** Found HYPERBOLIC {knot_name}: L={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                                        found_match_this_braid = True
                                        # Don't break, continue checking (unlikely match multiple)
                                except Exception as e_iso:
                                    # Expected errors for non-matching manifolds
                                    error_types['isometry_check'] += 1
                                    # Suppress common error message (already printed at start)
                                    common_error = "The SnapPea kernel was not able to determine if the manifolds are isometric."
                                    if str(e_iso) != common_error:
                                         errors_encountered += 1 # Only count uncommon errors?
                                         print(f"UNCOMMON ISOMETRY ERROR for {knot_name}: {e_iso}")
                                    pass # Continue checking other hyperbolic targets
                
                else: # Volume <= VOLUME_THRESHOLD
                    # --- Non-Hyperbolic Path: Use Jones Polynomial (if Sage available) --- 
                    identification_method = "Polynomial"
                    if not SAGE_AVAILABLE or not target_polynomials:
                        error_types['poly_check_skipped'] += 1
                        pass # Cannot check or no non-hyperbolic targets
                    else:
                        # Calculate polynomial for the candidate braid
                        candidate_poly = get_jones_polynomial(link_obj, braid_tuple, poly_assertion_state) # Pass state
                        if candidate_poly is None:
                            # Error counter is incremented based on where None comes from
                            # If it was assertion error, get_jones_poly handles print limit
                            error_types['poly_calc_failed'] += 1 # Generic failure count
                            errors_encountered += 1
                            # if count_processed <= 10: print(f"Failed poly calc for braid: {braid_tuple}")
                            continue # Cannot compare if calculation failed

                        for knot_name, target_poly in target_polynomials.items():
                            # Check if it's a potential target and shorter rep needed
                            if knot_name in target_knots_info and \
                               (knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']):
                                try:
                                    # Compare candidate polynomial with target
                                    is_match = (candidate_poly == target_poly)
                                    
                                    # Also check mirror image if direct match fails
                                    mirror_poly = None
                                    if not is_match:
                                        target_link = target_links.get(knot_name)
                                        if target_link:
                                            try:
                                                mirror_link = target_link.mirror()
                                                # Pass dummy tuple for mirror pre-calc check
                                                mirror_poly = get_jones_polynomial(mirror_link, dummy_tuple, poly_assertion_state) # Pass state
                                                if mirror_poly is not None:
                                                    is_match = (candidate_poly == mirror_poly)
                                                    # if is_match: print(f"  (Identified mirror image of {knot_name})")
                                            except Exception:
                                                error_types['mirror_poly_calc_failed'] += 1
                                                # Continue without mirror check if it fails

                                    if is_match:
                                        # Found a match!
                                        best_braid_repr[knot_name] = {
                                            'braid_tuple': braid_tuple,
                                            'length': current_length,
                                            'source': 'search_poly' + ('_mirror' if candidate_poly == mirror_poly else '')
                                        }
                                        count_candidates_found += 1
                                        print(f"  ** Found NON-HYPERBOLIC {knot_name}: L={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                                        found_match_this_braid = True
                                        # Don't break, continue checking (unlikely match multiple)

                                except Exception as e_poly_comp:
                                    error_types['poly_comparison_error'] += 1
                                    errors_encountered += 1
                                    print(f"ERROR during poly comparison for {knot_name}: {e_poly_comp}")
                                    pass # Continue checking other non-hyperbolic targets

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
            # Enable debug mode for 2-strand search to track the trefoil
            debug_mode = True
        elif n_strands == 3:
            max_len = 10
            debug_mode = False
        else:
            max_len = 7
            debug_mode = False
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

    # --- Debug Trefoil Detection ---
    print("\n" + "="*40)
    print("--- DEBUGGING TREFOIL DETECTION (Hybrid Method) ---")
    
    def debug_trefoil_detection(error_state):
        print("Testing explicit trefoil (3_1) recognition:")
        print("-" * 30)
        knot_name = "3_1"
        trefoil_braid = (1, 1, 1)
        volume = -1
        manifold = None
        poly = None

        print("Attempting Manifold/Volume/Isometry...")
        try:
            braid_obj = spherogram.ClosedBraid(*trefoil_braid)
            manifold = braid_obj.exterior()
            volume = manifold.volume()
            print(f"  Volume: {volume}")
            if volume > VOLUME_THRESHOLD:
                print("  Volume > threshold, attempting isometry check (EXPECTING FAILURE or ERROR)...")
                target_manifold = get_knot_manifold(knot_name)
                if target_manifold:
                    try:
                        result = manifold.is_isometric_to(target_manifold)
                        print(f"  Isometry check result: {result} (UNEXPECTED if True)")
                    except Exception as e:
                        print(f"  Isometry check expected error: {e}")
                else:
                    print("  Could not get target manifold for isometry check.")
            else:
                print(f"  Volume <= threshold ({VOLUME_THRESHOLD}), skipping isometry check (correct path).")
        except Exception as e:
            print(f"  Error during manifold/volume processing: {e}")
        
        print("\nAttempting Polynomial Comparison (requires Sage)... ")
        if not SAGE_AVAILABLE:
             print("  SKIPPING: Not running in Sage.")
        else:
            # Need dummy tuple for calls inside debug
            dummy_tuple = ('target',) 
            braid_obj = None # Ensure braid_obj is defined before try
            try:
                 braid_obj = spherogram.ClosedBraid(*trefoil_braid)
            except Exception as e:
                print(f"  Error creating braid object for poly test: {e}")

            try:
                # Get target polynomial
                target_poly = None
                target_link = None
                try:
                    target_link = spherogram.Link(knot_name)
                    target_poly = get_jones_polynomial(target_link, dummy_tuple, error_state)
                    if target_poly:
                         print(f"  Target Polynomial (3_1): {target_poly}")
                    else:
                         print("  Failed to get target polynomial for 3_1.")
                except Exception as e:
                    print(f"  Error getting target link/poly: {e}")

                # Get polynomial from braid
                if braid_obj and target_poly:
                    candidate_poly = get_jones_polynomial(braid_obj, trefoil_braid, error_state)
                    if candidate_poly:
                        print(f"  Candidate Polynomial ({trefoil_braid}): {candidate_poly}")
                        match = (candidate_poly == target_poly)
                        print(f"  Polynomials match: {match} (EXPECTING TRUE)")

                        # Check mirror
                        if not match:
                            print("  Checking mirror polynomial...")
                            mirror_poly = None
                            try:
                                mirror_link = target_link.mirror()
                                mirror_poly = get_jones_polynomial(mirror_link, dummy_tuple, error_state)
                                if mirror_poly:
                                     mirror_match = (candidate_poly == mirror_poly)
                                     print(f"  Mirror Polynomial (3_1_mirror): {mirror_poly}")
                                     print(f"  Candidate matches mirror: {mirror_match}")
                            except Exception as e_mirror:
                                print(f"  Error getting mirror poly: {e_mirror}")
                    else:
                         print("  Failed to calculate candidate polynomial from braid (1,1,1)." )
            except Exception as e:
                print(f"  Error during polynomial processing: {e}")

        print("\nTrefoil hybrid detection tests complete.")
        print("-" * 30)
    
    # Run the debug function
    # Define a state dictionary for the debug function's calls
    debug_poly_assertion_state = {'printed_count': 0, 'max_prints': 5, 'suppressed_msg_shown': False}
    debug_trefoil_detection(debug_poly_assertion_state) # Pass the state