import numpy as np
import math
import sys
import time
import signal
import collections
import re
# --- LSC Stability Filter (v15 -> v16: Expand Invariant Fits -> v17: Refactor Fits, Add Multi-Var) ---\n# --- Searching for Candidate Topologies & Correlating Invariants ---\n\nimport numpy as np\nimport math\nimport sys\nimport time\nimport signal\nimport collections\nimport re # For parsing braid strings if needed later
# import sympy # REMOVED - Using Sage polynomials directly

# Attempt import
try:
    import spherogram
    print("INFO: Spherogram library found.")
    import snappy
    print("INFO: SnapPy library found.")
    # Check for Sage environment for polynomial calculations
    import sage.all
    # from sage.symbolic.integration.integral import definite_integral # Test import - remove if not needed
    print("INFO: Running inside SageMath environment.")
    SAGE_AVAILABLE = True
    # Define polynomial variables for Sage (used in calculate_invariants)
    # Use q_var to avoid conflict with q as loop/value variable
    R = sage.all.LaurentPolynomialRing(sage.all.QQ, 'q_var')
    q_var = R.gen()
    T = sage.all.PolynomialRing(sage.all.ZZ, 't')
    t = T.gen()
    # Define roots of unity
    omega3 = sage.all.exp(2 * sage.all.pi * sage.all.I / 3)
    omega4 = sage.all.exp(2 * sage.all.pi * sage.all.I / 4) # == I
    omega5 = sage.all.exp(2 * sage.all.pi * sage.all.I / 5)
    omega6 = sage.all.exp(2 * sage.all.pi * sage.all.I / 6)
except ImportError:
    # Distinguish between missing Sage and other import errors
    try:
        # If spherogram/snappy are already imported, the missing one is Sage
        if 'spherogram' in sys.modules and 'snappy' in sys.modules:
             print("WARNING: Not running inside SageMath. Polynomial calculations will be skipped.")
        else:
             # Otherwise, it might be spherogram or snappy
             print(f"ERROR: Required library not found. Please install spherogram and snappy.")
    except Exception: # Fallback
        print(f"ERROR: Required library not found. Please install spherogram and snappy.")
    SAGE_AVAILABLE = False
except Exception as e:
    print(f"Warning: Issue detecting/setting up Sage environment: {e}")
    SAGE_AVAILABLE = False


print("\n--- Loading LSC Hybrid Search & Invariant Correlator ---")
print("--- v14: Added test for hypothesis: c1 ~ -log |Jones(K; q)| ---")
print("--- v15: Added test for refined hypothesis: c1 ~ -alpha * log |J(K; q)| ---")
print("--- v16: Expanded individual invariant correlation tests. ---")
print("--- v17: Refactored fits, added multi-variable fit capability. ---")

# --- Constants & Targets ---
alpha_fine_structure = 1 / 137.035999
planck_mass_kg = 2.176434e-8
TARGET_ENERGY_ELECTRON = 4.18540e-23
TARGET_ENERGY_MUON = 8.65429e-21
TARGET_ENERGY_TAU = 1.45532e-19
ASSUMED_C_TUNNEL = 1.0
def get_required_c1(E0_natural, C_tunnel=1.0):
    # --- Debug Print ---
    # print(f"DEBUG: get_required_c1 called with E0={E0_natural:.4e}")
    k_required = (2 * E0_natural)**2
    if k_required <= 0: return None
    try:
        log_arg = k_required / C_tunnel
        # --- Debug Print ---
        # print(f"DEBUG: k_required={k_required:.4e}, log_arg={log_arg:.4e}")
        if log_arg <= np.finfo(float).tiny:
            print(f"WARN: log_arg ({log_arg:.4e}) too small for log, returning inf for E0={E0_natural:.4e}")
            return np.inf # Or perhaps None is better?
        c1 = -alpha_fine_structure * np.log(log_arg)
        # --- Debug Print ---
        # print(f"DEBUG: Calculated c1={c1:.4f}")
        return c1 if c1 > 0 else None
    except Exception as e:
        print(f"ERROR in get_required_c1 for E0={E0_natural:.4e}: {e}") # Print exception
        return None

# C1 required based *purely* on energy targets (for hierarchy check)
C1_REQUIRED_TARGETS = {
    "electron": get_required_c1(TARGET_ENERGY_ELECTRON),
    "muon": get_required_c1(TARGET_ENERGY_MUON),
    "tau": get_required_c1(TARGET_ENERGY_TAU),
}
C1_REQUIRED_TARGETS = {k: v for k, v in C1_REQUIRED_TARGETS.items() if v is not None}

print("\n--- Target c1 Values (Based on Lepton Energy) ---")
if C1_REQUIRED_TARGETS:
    C1_TARGETS_SORTED = sorted(C1_REQUIRED_TARGETS.items(), key=lambda item: item[1], reverse=True)
    for p, c1 in C1_TARGETS_SORTED: print(f"  {p.capitalize():<9}: {c1:.3f}")
else: print("  Error calculating required c1 values from energy.")
print("----------------------------------------------------\n")


# --- Search Parameters (for find_simplest_braids) ---
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
NON_HYPERBOLIC_TARGETS = {"3_1", "5_1", "7_1"}
HYPERBOLIC_TARGETS = {k for k in TARGET_KNOTS_INFO if k not in NON_HYPERBOLIC_TARGETS}
VOLUME_THRESHOLD = 1e-6 # Threshold for checking volume (close to zero)

# Known braid representations (as tuples)
KNOWN_BRAIDS = {
    "3_1": (1, 1, 1),
    "4_1": (1, -2, 1, -2),
    "5_1": (1, 1, 1, 1, 1),
    "5_2": (1, -2, 1, -2, 1),
    "6_1": (-1, 2, -1, 2, -1, 2),
    "6_2": (1, -2, -2, 1, -2, -2),
    "6_3": (-1, -1, 2, -1, -1, 2),
    "7_1": (1, 1, 1, 1, 1, 1, 1),
}

MAX_BRAID_LENGTH = 12 # Default, adjusted per strand count later
N_STRANDS_TO_CHECK = [2, 3, 4]

# --- sklearn Import ---
try:
    from sklearn.linear_model import LinearRegression
    SKLEARN_AVAILABLE = True
    print("INFO: sklearn library found, multi-variable regression enabled.")
except ImportError:
    print("WARNING: sklearn library not found. Multi-variable regression will be skipped.")
    SKLEARN_AVAILABLE = False

# --- Helper Functions ---
def get_knot_manifold(name):
    try: return snappy.Manifold(name)
    except Exception as e: print(f"Warning: Manifold failed for {name}: {e}"); return None
def format_braid_word(braid_tuple): return ' '.join([f"s{abs(g)}{'^-1' if g < 0 else ''}" for g in braid_tuple])

# --- Braid Generation (kept from previous version) ---
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
        
        # if debug and n_strands == 2 and current_len <= 3:
        #     print(f"\nDEBUG: Processing braids of length {current_len}...")
        
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
                    # print(f"\nDEBUG: TREFOIL BRAID (1,1,1) GENERATED at position {total_yielded + count_yielded_this_level + 1}")
                # --- End debug check ---
                
                # if debug and n_strands == 2 and current_len <= 3:
                #     print(f"DEBUG: Yielding Length {current_len} braid: {new_word_tuple}")
                
                yield new_word_tuple
                total_yielded += 1
                queue.append(new_word_tuple)
                count_yielded_this_level += 1
        
        # end_level_time = time.time()
        # rate = count_yielded_this_level / (end_level_time - start_level_time + 1e-9)
        
        # if debug and n_strands == 2:
        #     print(f"DEBUG: Generated {count_yielded_this_level} braids of length {current_len} ({rate:.1f}/s)")
        #     print(f"DEBUG: Total braids so far: {total_yielded}")
        
        if not queue: break


# --- Timeout Handler ---
class TimeoutException(Exception): pass
def timeout_handler(signum, frame): print("\nERROR: Processing timed out!"); raise TimeoutException

# --- Function to get Jones polynomial safely (Requires Sage) ---
def get_jones_polynomial(link_obj, braid_tuple, error_state):
    if not SAGE_AVAILABLE:
        return None # Cannot calculate outside Sage
    try:
        # Use the Sage polynomial variable 'q_var' defined globally
        poly = link_obj.jones_polynomial(variable=q_var)
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
        # Print more details for general exceptions during polynomial calculation
        print(f"      ERROR calculating Jones polynomial for braid {braid_tuple}: {e} (Type: {type(e)})")
        return None

# --- Function to get Alexander polynomial safely (Requires Sage) ---
def get_alexander_polynomial(link_obj, braid_tuple, error_state):
    if not SAGE_AVAILABLE:
        return None
    try:
        # Use the Sage polynomial variable 't' defined globally
        poly = link_obj.alexander_polynomial(variable=t)
        return poly
    except spherogram.sage_helper.SageNotAvailable as e:
        print(f"      ERROR: Sage not available for Alexander polynomial: {e}")
        return None
    # Add assertion handling similar to Jones if needed?
    except Exception as e:
        print(f"      ERROR calculating Alexander polynomial for braid {braid_tuple}: {e} (Type: {type(e)})")
        return None

# --- Main Braid Finding Function (Modified Analysis Section) ---
def find_simplest_braids(n_strands, max_braid_length, target_knots_info, time_limit_seconds):
    print(f"\n--- Running Search ({n_strands} STRANDS, MaxLen={max_braid_length}, Timeout={time_limit_seconds}s) ---")
    print(f"--- Using HYBRID Identification (Manifold Vol > {VOLUME_THRESHOLD} ? SnapPy : Jones Poly) ---")
    has_signal = hasattr(signal, 'SIGALRM'); timer_set = False
    if has_signal:
        try: signal.signal(signal.SIGALRM, timeout_handler); signal.alarm(time_limit_seconds); timer_set = True
        except ValueError: print("Warning: Cannot set SIGALRM timer.")

    start_time = time.time()
    best_braid_repr = {} # knot_name -> {'braid_tuple': tuple, 'length': L, 'source': 'known'/'search_iso'/'search_poly'}
    count_processed = 0; count_candidates_found = 0
    timed_out = False; errors_encountered = 0
    last_print_time = time.time()
    error_types = collections.Counter()

    # State for managing assertion error printing
    poly_assertion_state = {'printed_count': 0, 'max_prints': 5, 'suppressed_msg_shown': False}

    # --- Pre-calculation (Manifolds & Polynomials) ---
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
        dummy_tuple = ('target',) # Dummy for error reporting during pre-calc
        for name in NON_HYPERBOLIC_TARGETS:
            if name in target_knots_info:
                try:
                    link = spherogram.Link(name)
                    target_links[name] = link # Store link for mirror check later
                    poly = get_jones_polynomial(link, dummy_tuple, poly_assertion_state) # Pass state
                    if poly is not None:
                        target_polynomials[name] = poly
                    else:
                        print(f"  WARNING: Failed to calculate Jones polynomial for non-hyperbolic target {name}.")
                        poly_errors += 1
                except Exception as e:
                    print(f"  ERROR: Failed to create Link object for non-hyperbolic target {name}: {e}")
                    poly_errors += 1
        print(f"Pre-calculated {len(target_polynomials)} polynomials for {len(NON_HYPERBOLIC_TARGETS)} non-hyperbolic targets.")

    if manifold_errors > 0 or poly_errors > 0:
        print(f"WARNING: Encountered {manifold_errors} manifold / {poly_errors} poly errors during pre-calculation.")
    if not target_manifolds and not target_polynomials:
        print("ERROR: No target manifolds or polynomials loaded/calculated. Cannot identify knots.")
        return {}
    print("-------------------------------------------\n")

    # --- Braid Generation & Filtering Loop ---
    enable_debug = (n_strands == 2) # Less debug printing needed now
    braid_generator = generate_braid_words_optimized(n_strands, max_braid_length, debug=False)
    print("\nFiltering generated braids using Hybrid Method...")

    try:
        for braid_tuple in braid_generator:
            count_processed += 1
            # Progress update
            current_time = time.time()
            if current_time - last_print_time > 10.0:
                rate = count_processed / (current_time - start_time + 1e-9)
                print(f"  Processed {count_processed} braids... ({rate:.1f} braids/sec)")
                # print(f"  Error types so far: {dict(error_types)}") # Can be noisy
                last_print_time = current_time

            # --- Basic Braid Filters ---
            if not braid_tuple: continue

            # Optimization: Check if this braid is longer than shortest found
            all_target_names = HYPERBOLIC_TARGETS.union(NON_HYPERBOLIC_TARGETS)
            if len(best_braid_repr) >= len(all_target_names) and \
               all(data['length'] <= len(braid_tuple) for data in best_braid_repr.values()):
               continue # Skip if cannot improve any result

            # Additional structure filters (sum=0, cancellation)
            if sum(braid_tuple) == 0:
                error_types['sum_zero_filter'] += 1
                continue
            if len(braid_tuple) >= 2:
                has_cancellation = False
                for i in range(len(braid_tuple)-1):
                    if braid_tuple[i] == -braid_tuple[i+1]:
                        has_cancellation = True; break
                if has_cancellation:
                    error_types['cancellation_filter'] += 1
                    continue

            # --- Create Link Object ---
            try:
                # Handle single-element tuples specially
                braid_list = [braid_tuple[0]] if len(braid_tuple) == 1 else list(braid_tuple)
                link_obj = spherogram.ClosedBraid(braid_list)
            except Exception as e:
                error_types['braid_creation'] += 1; errors_encountered += 1
                # if count_processed <= 5: print(f"ERROR on Braid creation: {e} for braid {braid_tuple}")
                continue

            # --- Knot/Link Pre-filters ---
            try:
                num_components = len(link_obj.link_components)
                if num_components != 1:
                    error_types['not_a_knot'] += 1
                    continue # Only interested in knots for now
            except Exception as e:
                error_types['link_components_error'] += 1; errors_encountered += 1
                # if count_processed <= 5: print(f"ERROR checking link components: {e} for braid {braid_tuple}")
                continue

            try:
                # Optional crossing number check (can be slow/error-prone)
                if hasattr(link_obj, 'crossing_number'):
                    crossings = link_obj.crossing_number()
                    if crossings > MAX_CANDIDATE_CROSSINGS + 4:
                        error_types['crossing_num_filter'] += 1
                        continue
            except Exception:
                error_types['crossing_num_error'] += 1
                pass # Proceed if crossing number fails

            # --- Hybrid Identification ---
            manifold = None; volume = -1.0
            try:
                manifold = link_obj.exterior()
                volume = manifold.volume()
            except Exception as e:
                error_types['exterior_vol_error'] += 1; errors_encountered += 1
                continue # Cannot identify if exterior/volume fails

            found_match_this_braid = False
            current_length = len(braid_tuple)

            # --- Branch based on Volume ---
            if volume > VOLUME_THRESHOLD:
                # --- Hyperbolic Path: Use Snappy Isometry ---
                if target_manifolds:
                    for knot_name, target_manifold in target_manifolds.items():
                        if knot_name in target_knots_info and \
                           (knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']):
                            try:
                                if manifold.is_isometric_to(target_manifold):
                                    best_braid_repr[knot_name] = {
                                        'braid_tuple': braid_tuple, 'length': current_length, 'source': 'search_iso'
                                    }
                                    count_candidates_found += 1
                                    print(f"  ** Found HYPERBOLIC {knot_name}: L={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                                    found_match_this_braid = True
                            except Exception as e_iso:
                                error_types['isometry_check'] += 1
                                common_error = "The SnapPea kernel was not able to determine if the manifolds are isometric."
                                if str(e_iso) != common_error:
                                     errors_encountered += 1
                                     # print(f"UNCOMMON ISOMETRY ERROR for {knot_name}: {e_iso}")
                                pass
            else: # Volume <= VOLUME_THRESHOLD
                # --- Non-Hyperbolic Path: Use Jones Polynomial ---
                if SAGE_AVAILABLE and target_polynomials:
                    candidate_poly = get_jones_polynomial(link_obj, braid_tuple, poly_assertion_state)
                    if candidate_poly is not None:
                        for knot_name, target_poly in target_polynomials.items():
                             if knot_name in target_knots_info and \
                                (knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']):
                                try:
                                    is_match = (candidate_poly == target_poly)
                                    mirror_poly = None
                                    if not is_match:
                                        # Check mirror (more efficient pre-calc needed if done often)
                                        target_link = target_links.get(knot_name)
                                        if target_link:
                                            try:
                                                mirror_link = target_link.mirror()
                                                mirror_poly = get_jones_polynomial(mirror_link, ('mirror', knot_name), poly_assertion_state)
                                                if mirror_poly is not None: is_match = (candidate_poly == mirror_poly)
                                            except Exception: error_types['mirror_poly_calc_failed'] += 1
                                    if is_match:
                                        source = 'search_poly' + ('_mirror' if mirror_poly and candidate_poly == mirror_poly else '')
                                        best_braid_repr[knot_name] = {
                                            'braid_tuple': braid_tuple, 'length': current_length, 'source': source
                                        }
                                        count_candidates_found += 1
                                        print(f"  ** Found NON-HYPERBOLIC {knot_name}: L={current_length} (Braid: {format_braid_word(braid_tuple)}) {'(Mirror)' if '_mirror' in source else ''} **")
                                        found_match_this_braid = True
                                except Exception as e_poly_comp:
                                    error_types['poly_comparison_error'] += 1; errors_encountered += 1
                                    # print(f"ERROR during poly comparison for {knot_name}: {e_poly_comp}")
                                    pass
                    else: # candidate_poly is None
                         error_types['poly_calc_failed'] += 1
                         # errors_encountered incremented inside get_jones_polynomial
                         pass
                else: # Sage not available or no non-hyperbolic targets
                     error_types['poly_check_skipped'] += 1
                     pass

    except TimeoutException: timed_out = True
    except Exception as e_outer:
        print(f"\nMAJOR ERROR in braid loop: {e_outer}")
        error_types['outer_loop_exception'] += 1
        errors_encountered += 1
        # Consider re-raising or exiting depending on severity
        # raise e_outer
    finally:
        if has_signal and timer_set: signal.alarm(0)

    # --- Post-Search Reporting ---
    end_time = time.time(); elapsed_time = end_time - start_time
    rate = count_processed / (elapsed_time + 1e-9)
    print(f"\nFinished search ({n_strands} strands, up to maxlen={max_braid_length}) in {elapsed_time:.2f} seconds.")
    print(f"Processed {count_processed} braids at {rate:.1f} braids/second.")
    if timed_out: print(f"Warning: Process timed out after {time_limit_seconds} seconds.")
    if errors_encountered > 0:
        print(f"Warning: Encountered {errors_encountered} errors/skips during processing.")
        print(f"Filtered by: Sum=0({error_types['sum_zero_filter']}), Cancel({error_types['cancellation_filter']}), Crossing#({error_types['crossing_num_filter']}), NotKnot({error_types['not_a_knot']})")
        print(f"Errors: Braid({error_types['braid_creation']}), Components({error_types['link_components_error']}), CrossingCalc({error_types['crossing_num_error']}), Ext/Vol({error_types['exterior_vol_error']}), Iso({error_types['isometry_check']}), PolyCalc({error_types['poly_calc_failed']}), PolyComp({error_types['poly_comparison_error']}), Mirror({error_types['mirror_poly_calc_failed']}), Other({error_types['outer_loop_exception']})")

    # --- Analysis Within Function (Hierarchy Check) ---
    print(f"\n--- Analysis: Simplest Braid Reps Found ({n_strands} Strands) ---")
    if not best_braid_repr: print("No target knots found.")
    else:
        known_achiral = {"4_1", "6_2", "6_3"}
        sorted_found_knots = sorted(best_braid_repr.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        print("Simplest braid representations found/known for target knots:")
        for knot_name in sorted_found_knots:
             data = best_braid_repr[knot_name]
             nc, sig = TARGET_KNOTS_INFO[knot_name]
             source = data.get('source', 'search')
             is_achiral = (knot_name in known_achiral)
             chiral_str = "(Achiral)" if is_achiral else "(Chiral)"
             print(f"  - {knot_name:<4} (Nc={nc}, Sig={sig}) {chiral_str}: Min Braid Length={data['length']} [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")

        print("\n--- Checking Lepton Assignment Consistency (Using Chirality Filter for this strand count) ---")
        target_c1_list = sorted(C1_REQUIRED_TARGETS.items(), key=lambda item: item[1], reverse=True) # Use targets based on energy

        chiral_knot_candidates = [(name, data) for name, data in best_braid_repr.items() if name in TARGET_KNOTS_INFO and name not in known_achiral]
        chiral_knot_list_sorted_by_Nc = sorted(chiral_knot_candidates, key=lambda item: TARGET_KNOTS_INFO[item[0]][0])
        print(f"Found {len(chiral_knot_list_sorted_by_Nc)} potentially stable chiral knot candidates for this strand count.")

        consistent = True
        if len(chiral_knot_list_sorted_by_Nc) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct stable CHIRAL KNOT types found ({len(chiral_knot_list_sorted_by_Nc)}) to map to {len(target_c1_list)} lepton generations for this strand count.")
             consistent = False
        else:
            print("Hypothetical Assignment (Increasing Chiral Knot Nc -> Increasing Mass / Decreasing c1):")
            complexities = []
            assigned_knots_local = {}
            for i in range(len(target_c1_list)):
                 lepton, req_c1 = target_c1_list[i]
                 knot_name, data = chiral_knot_list_sorted_by_Nc[i]
                 nc = TARGET_KNOTS_INFO[knot_name][0]
                 complexities.append(nc)
                 assigned_knots_local[lepton] = knot_name
                 print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")

            if not all(complexities[i] <= complexities[i+1] for i in range(len(complexities)-1)):
                 print("  >> Warning: Complexity (Nc) ordering inconsistent with lepton hierarchy.")
                 consistent = False
            else:
                 print("  >> Observation: Complexity (Nc) ordering of simplest chiral knots IS consistent for this strand count.")
            # ... (Rest of success/failure printing for local assignment) ...
            if consistent:
                 print(f"  >> SUCCESS: Plausible assignment found for this strand count:")
                 print(f"     Electron -> {assigned_knots_local.get('electron', 'N/A')}")
                 print(f"     Muon     -> {assigned_knots_local.get('muon', 'N/A')}")
                 print(f"     Tau      -> {assigned_knots_local.get('tau', 'N/A')}")
                 print(f"  >> Assessment: Requires derivation of c1(Topology) yielding correct values for these knots.")
            else:
                 print("  >> FAILURE: Consistent assignment failed with current candidates/ordering for this strand count.")

    return best_braid_repr


# --- Debug Function Definition (Moved to top level) ---
# Copied from previous working version
def debug_trefoil_detection(error_state):
    print("Testing explicit trefoil (3_1) recognition:")
    print("-" * 30)
    knot_name = "3_1"
    trefoil_braid = (1, 1, 1)
    volume = -1
    manifold = None
    poly = None
    braid_obj = None # Initialize braid_obj

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
        # braid_obj might have been created above

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
            if braid_obj is None: # Try creating if manifold part failed
                 try: braid_obj = spherogram.ClosedBraid(*trefoil_braid)
                 except: print("  Failed to create braid_obj for poly test")

            if braid_obj and target_poly:
                candidate_poly = get_jones_polynomial(braid_obj, trefoil_braid, error_state)
                if candidate_poly:
                    print(f"  Candidate Polynomial ({trefoil_braid}): {candidate_poly}")
                    match = (candidate_poly == target_poly)
                    # Note: Expecting FALSE due to Spherogram convention where Link("3_1") != Braid(1,1,1)
                    print(f"  Polynomials match: {match} (EXPECTING FALSE due to convention)")

                    # Check mirror
                    if not match and target_link: # Need target_link for mirror
                        print("  Checking mirror polynomial...")
                        mirror_poly = None
                        try:
                            mirror_link = target_link.mirror()
                            mirror_poly = get_jones_polynomial(mirror_link, dummy_tuple, error_state)
                            if mirror_poly:
                                 mirror_match = (candidate_poly == mirror_poly)
                                 print(f"  Mirror Polynomial (3_1_mirror): {mirror_poly}")
                                 # Note: Expecting TRUE as Braid(1,1,1) should match Mirror(Link("3_1"))
                                 print(f"  Candidate matches mirror: {mirror_match} (EXPECTING TRUE)")
                        except Exception as e_mirror:
                            print(f"  Error getting mirror poly: {e_mirror}")
                else:
                     print("  Failed to calculate candidate polynomial from braid (1,1,1)." )
        except Exception as e:
            print(f"  Error during polynomial processing: {e}")

    print("\nTrefoil hybrid detection tests complete.")
    print("-" * 30)


# --- Invariant Calculation Function (Cleaned & Enhanced for Log|J(q)|) ---
def calculate_invariants(knot_name):
    """Calculate key invariants for a given knot using Spherogram."""
    print(f"Calculating invariants for {knot_name}...")
    # Initialize data with defaults
    data = {'name': knot_name, 'Nc': TARGET_KNOTS_INFO.get(knot_name, (None,None))[0], # Get Nc from dict
            'signature': None, 'determinant': None,
            # 'arf': None, # Removed - not available directly
            'jones_at_5th_root': None, # Keep this |J(ω5)| for comparison
            'alexander_at_minus1': None,
            'alexander_at_minus2': None,
            'log_abs_determinant': None, # NEW
            'log_abs_jones_w3': None, # NEW: log |J(ω3)|
            'log_abs_jones_w4': None, # NEW: log |J(ω4)|
            'log_abs_jones_w5': None, # NEW: log |J(ω5)|
            'log_abs_jones_w6': None, # NEW: log |J(ω6)|
            # 'is_chiral': None # Removed - not available directly
           }
    link = None # Define link outside try block

    try:
        link = spherogram.Link(knot_name)
        # --- Attempt to generate diagram data (still useful for other methods?) ---
        # try:
        #     # Simplify might compute necessary diagram info
        #     print(f"  Attempting link.simplify() for {knot_name}...")
        #     link.simplify()
        #     print(f"  ...simplify() done.")
        # except Exception as e_simp:
        #     print(f"  Warning: link.simplify() failed for {knot_name}: {e_simp}")

        # Calculate invariants known to work
        try: data['signature'] = link.signature()
        except AttributeError: print(f"  Warn: Cannot get signature for {knot_name}")
        except Exception as e: print(f"  Error getting signature for {knot_name}: {e}")

        # --- Determinant Calculation ---
        det_value = None # Temporary variable
        try:
            det_value = abs(link.determinant())
            data['determinant'] = det_value # Store original determinant
        except AttributeError:
            print(f"  Warn: Cannot get determinant for {knot_name}")
        except Exception as e:
            print(f"  Error getting determinant for {knot_name}: {e}")

        # --- Log Determinant Calculation (only if determinant was successful) ---
        if det_value is not None:
            if det_value > 1e-9:
                try: data['log_abs_determinant'] = float(np.log(det_value))
                except Exception as e_logdet: print(f"  Warn: log(Det) failed: {e_logdet}")
            else:
                 print(f"  Warn: Determinant ({det_value}) <= 1e-9, skipping log.")
        # else: Determinant calculation failed, log automatically skipped

        # --- Polynomial Calculations ---
        if SAGE_AVAILABLE:
            try:
                # Jones polynomial calculation (raw polynomial)
                jones_poly_q = link.jones_polynomial(variable=q_var) # Use q_var

                # Evaluate at omega5 = exp(2*pi*i / 5)
                jones_eval_w5 = jones_poly_q(omega5)
                jones_abs_w5 = abs(jones_eval_w5) if hasattr(jones_eval_w5,'abs') else jones_eval_w5
                data['jones_at_5th_root'] = float(jones_abs_w5) # Keep original |J(ω5)|
                if jones_abs_w5 > 1e-9: # Avoid log(0)
                    data['log_abs_jones_w5'] = float(np.log(jones_abs_w5))
                else: print(f"  Warn: |J(ω5)| is near zero for {knot_name}")

                # Evaluate at omega4 = exp(2*pi*i/4) = i
                jones_eval_w4 = jones_poly_q(omega4)
                jones_abs_w4 = abs(jones_eval_w4) if hasattr(jones_eval_w4,'abs') else jones_eval_w4
                if jones_abs_w4 > 1e-9: # Avoid log(0)
                    data['log_abs_jones_w4'] = float(np.log(jones_abs_w4))
                else: print(f"  Warn: |J(ω4)| is near zero for {knot_name}")

                # Evaluate at omega3 = exp(2*pi*i/3)
                jones_eval_w3 = jones_poly_q(omega3)
                jones_abs_w3 = abs(jones_eval_w3) if hasattr(jones_eval_w3,'abs') else jones_eval_w3
                if jones_abs_w3 > 1e-9: # Avoid log(0)
                    data['log_abs_jones_w3'] = float(np.log(jones_abs_w3))
                else: print(f"  Warn: |J(ω3)| is near zero for {knot_name}")

                # Evaluate at omega6 = exp(2*pi*i/6)
                jones_eval_w6 = jones_poly_q(omega6)
                jones_abs_w6 = abs(jones_eval_w6) if hasattr(jones_eval_w6,'abs') else jones_eval_w6
                if jones_abs_w6 > 1e-9: # Avoid log(0)
                    data['log_abs_jones_w6'] = float(np.log(jones_abs_w6))
                else: print(f"  Warn: |J(ω6)| is near zero for {knot_name}")

            except Exception as e_jones:
                 print(f"  Warning: Jones calculation failed for {knot_name}: {e_jones}")

            try:
                # Alexander polynomial at t=-1 and t=-2 - requires Sage
                alex_poly_t = link.alexander_polynomial() # Call without variable
                if alex_poly_t is not None and not alex_poly_t.is_zero():
                     try: data['alexander_at_minus1'] = abs(alex_poly_t(-1))
                     except Exception as e_eval1: print(f"  Warn: Alex(-1) eval failed: {e_eval1}")

                     try: data['alexander_at_minus2'] = abs(alex_poly_t(-2))
                     except Exception as e_eval2: print(f"  Warn: Alex(-2) eval failed: {e_eval2}")
                else:
                     print(f"  Warn: Alexander polynomial is zero or None for {knot_name}, cannot evaluate.")
            except Exception as e_alex:
                 print(f"  Warning: Alexander calculation failed for {knot_name}: {e_alex}")

    except Exception as e_link:
        print(f"Error creating/processing Link object for {knot_name}: {e_link}")

    # Format output string
    j5_str = f"{data['jones_at_5th_root']:.3f}" if data['jones_at_5th_root'] is not None else "N/A"
    a1_str = f"{data['alexander_at_minus1']}" if data['alexander_at_minus1'] is not None else "N/A" # Integer usually
    a2_str = f"{data['alexander_at_minus2']}" if data['alexander_at_minus2'] is not None else "N/A" # Integer usually
    ljw3_str = f"{data['log_abs_jones_w3']:.3f}" if data['log_abs_jones_w3'] is not None else "N/A"
    ljw4_str = f"{data['log_abs_jones_w4']:.3f}" if data['log_abs_jones_w4'] is not None else "N/A"
    ljw5_str = f"{data['log_abs_jones_w5']:.3f}" if data['log_abs_jones_w5'] is not None else "N/A"
    ljw6_str = f"{data['log_abs_jones_w6']:.3f}" if data['log_abs_jones_w6'] is not None else "N/A"
    ldet_str = f"{data['log_abs_determinant']:.3f}" if data['log_abs_determinant'] is not None else "N/A"

    # Updated print string to include log|J(q)| values
    print(f"  -> Nc={data['Nc']}, Sig={data['signature']}, Det={data['determinant']}, log|Det|={ldet_str}, "
          f"Δ(-1)={a1_str}, Δ(-2)={a2_str}, "
          f"log|J(ω3)|={ljw3_str}, log|J(ω4)|={ljw4_str}, log|J(ω5)|={ljw5_str}, log|J(ω6)|={ljw6_str}")
    return data


# --- Main Execution ---
if __name__ == "__main__":
    timeout_long = 180 # 3 minutes per strand count

    # --- Braid Search Phase ---
    print("\n" + "="*20 + " PHASE 1: BRAID SEARCH " + "="*20)
    all_results = {}
    for n_strands in N_STRANDS_TO_CHECK:
        # Adjust max_len based on strands
        if n_strands == 2: max_len = 12
        elif n_strands == 3: max_len = 10
        else: max_len = 7
        print("\n" + "="*10 + f" {n_strands}-STRAND DEEP SEARCH (MaxLen={max_len}, Timeout={timeout_long}s) " + "="*10)
        results = find_simplest_braids(n_strands=n_strands, max_braid_length=max_len,
                                      target_knots_info=TARGET_KNOTS_INFO, time_limit_seconds=timeout_long)
        all_results[n_strands] = results

    # --- Final Braid Summary & Hierarchy Check ---
    print("\n" + "="*20 + " PHASE 2: FINAL BRAID SUMMARY & HIERARCHY " + "="*20)
    combined_best = {}
    for n_strands, results in all_results.items():
         for name, data in results.items():
             if name not in combined_best or data['length'] < combined_best[name]['length']:
                 combined_best[name] = data; combined_best[name]['strands'] = n_strands

    final_assigned_knots = {} # Initialize final assignment dictionary
    if not combined_best:
        print("No target knots found in any run.")
    else:
        sorted_final = sorted(combined_best.keys(), key=lambda k: TARGET_KNOTS_INFO[k][0])
        print("Simplest representations found overall:")
        known_achiral = {"4_1", "6_2", "6_3"}
        for name in sorted_final:
            data = combined_best[name]; nc, sig = TARGET_KNOTS_INFO[name]
            source = data.get('source', 'search')
            chiral_str = "(Achiral)" if name in known_achiral else "(Chiral)"
            print(f"  - {name:<4} (Nc={nc}, Sig={sig}) {chiral_str}: Min Braid Length={data['length']} ({data['strands']} strands) [{source.upper()}] (Braid: {format_braid_word(data['braid_tuple'])})")

        print("\n--- FINAL Lepton Assignment Check (Based on Simplest Found Chiral Knots) ---")
        target_c1_list = sorted(C1_REQUIRED_TARGETS.items(), key=lambda item: item[1], reverse=True)

        final_chiral_knot_candidates = [(name, data) for name, data in combined_best.items() if name in TARGET_KNOTS_INFO and name not in known_achiral]
        final_chiral_list_sorted_by_Nc = sorted(final_chiral_knot_candidates, key=lambda item: TARGET_KNOTS_INFO[item[0]][0])

        print(f"Found {len(final_chiral_list_sorted_by_Nc)} distinct chiral knots overall.")

        final_consistent = True
        if len(final_chiral_list_sorted_by_Nc) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct stable CHIRAL KNOT types found overall ({len(final_chiral_list_sorted_by_Nc)}) to map to {len(target_c1_list)} lepton generations.")
             final_consistent = False
        else:
            print("Final Hypothetical Assignment (Increasing Found Chiral Knot Nc -> Increasing Mass / Decreasing c1):")
            final_complexities = []
            for i in range(len(target_c1_list)):
                 lepton, req_c1 = target_c1_list[i]
                 knot_name, data = final_chiral_list_sorted_by_Nc[i]
                 nc = TARGET_KNOTS_INFO[knot_name][0]
                 final_complexities.append(nc)
                 final_assigned_knots[lepton] = knot_name
                 print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")

            if not all(final_complexities[i] <= final_complexities[i+1] for i in range(len(final_complexities)-1)):
                 print("  >> FINAL Warning: Complexity (Nc) ordering inconsistent with lepton hierarchy.")
                 final_consistent = False
            else:
                 print("  >> FINAL Observation: Complexity (Nc) ordering of overall simplest chiral knots IS consistent.")

            if final_consistent:
                 print(f"  >> FINAL SUCCESS: Plausible assignment based on complexity ordering found:")
                 # ... (print assignments) ...
                 print(f"     Electron -> {final_assigned_knots.get('electron', 'N/A')}")
                 print(f"     Muon     -> {final_assigned_knots.get('muon', 'N/A')}")
                 print(f"     Tau      -> {final_assigned_knots.get('tau', 'N/A')}")
                 print(f"  >> Assessment: Requires derivation of c1(Topology) yielding correct values for these knots.")
            else:
                 print("  >> FINAL FAILURE: Consistent assignment failed with overall best candidates based on complexity.")

    # --- Invariant Calculation Phase ---
    print("\n" + "="*20 + " PHASE 3: INVARIANT CALCULATION " + "="*20)
    KNOTS_TO_ANALYZE = ["3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3", "7_1", "7_2"] # Expand list?
    knot_invariant_data = collections.OrderedDict()
    print("Calculating invariants for selected knots...")
    for name in KNOTS_TO_ANALYZE:
        # Check if knot was found in braid search to potentially reuse Link object? No, recalculate for simplicity.
        knot_invariant_data[name] = calculate_invariants(name)
    print("\n--- Invariant Calculation Finished ---")

    # --- Correlation Search Phase ---
    print("\n" + "="*20 + " PHASE 4: CORRELATION SEARCH " + "="*20)

    # --- Subsection 4.1: Fit based on FIXED Hypothesis (3_1, 5_1, 5_2) ---
    # --- Subsection 4.1: Fit based on FIXED Hypothesis (3_1, 5_1, 5_2) ---
    print("\n--- 4.1: Correlation Test (Fixed Hypothesis: 3_1(e), 5_1(mu), 5_2(tau)) ---")
    C1_REQUIRED_FOR_FIXED_FIT = {
        "3_1": C1_REQUIRED_TARGETS.get("electron"),
        "5_1": C1_REQUIRED_TARGETS.get("muon"),
        "5_2": C1_REQUIRED_TARGETS.get("tau"),
    }
    C1_REQUIRED_FOR_FIXED_FIT = {k: v for k, v in C1_REQUIRED_FOR_FIXED_FIT.items() if v is not None}

    leptons_for_fixed_fit = ["3_1", "5_1", "5_2"]
    fit_data_points_fixed = {name: knot_invariant_data.get(name) for name in leptons_for_fixed_fit}

    # Check if all required base data is present for fixed fit
    if any(data is None for data in fit_data_points_fixed.values()) or len(C1_REQUIRED_FOR_FIXED_FIT) != 3 :
        print("ERROR: Missing base invariant data or c1 targets for knots in the fixed assignment (3_1, 5_1, 5_2). Cannot perform fixed fits.")
        fixed_fit_attempted = False
    else:
        fixed_fit_attempted = True
        target_c1_values_fixed = [C1_REQUIRED_FOR_FIXED_FIT[name] for name in leptons_for_fixed_fit]
        target_c1_array_fixed = np.array(target_c1_values_fixed)

        # --- Perform Individual Linear Fits: c1 = A * Inv + B ---
        print("\nAttempting individual linear fits (c1 vs Invariant) using fixed hypothesis (3_1, 5_1, 5_2):")
        invariants_to_fit = {
            "|Signature|": 'signature',
            "log|Det|": 'log_abs_determinant',
            "|Δ(-1)|": 'alexander_at_minus1',
            "|Δ(-2)|": 'alexander_at_minus2',
            "-log|J(ω3)|": 'log_abs_jones_w3',
            "-log|J(ω4)|": 'log_abs_jones_w4',
            "-log|J(ω5)|": 'log_abs_jones_w5',
            "-log|J(ω6)|": 'log_abs_jones_w6',
        }

        offset_threshold = 0.05 # Threshold for intercept check

        for inv_name, data_key in invariants_to_fit.items():
            print(f"\n--- Fit vs {inv_name} ---")
            x_values = []
            y_values = target_c1_array_fixed # Use the pre-calculated target c1 array
            valid_invariant_set = True # Flag for this invariant type

            # --- Loop through each knot for the CURRENT invariant ---
            for i, name in enumerate(leptons_for_fixed_fit):
                data_point = fit_data_points_fixed[name]
                raw_value = data_point.get(data_key)
                processed_value = None

                # Check for missing value first
                if raw_value is None:
                    print(f"  Skipping fit: Missing '{data_key}' for {name}.")
                    valid_invariant_set = False
                    break # Stop processing this invariant type

                try:
                    # --- Try processing the single value ---
                    # Ensure numeric type before transformations
                    current_val_numeric = None
                    if not isinstance(raw_value, (int, float, np.number)):
                        try:
                            current_val_numeric = float(raw_value)
                        except (ValueError, TypeError):
                             print(f"  Skipping fit: Non-numeric value '{raw_value}' for {inv_name} in {name}.")
                             valid_invariant_set = False
                             break # Stop processing this invariant type
                    else:
                        current_val_numeric = raw_value # Already numeric

                    # Apply transformations based on inv_name
                    if inv_name.startswith("|Sig") or inv_name.startswith("|Δ"):
                        processed_value = abs(current_val_numeric)
                    elif inv_name.startswith("-log|J"):
                        # raw_value here is log|J(ω)|. Check it's valid (> -inf)
                        if current_val_numeric > -np.inf:
                           processed_value = -current_val_numeric
                        else:
                            print(f"  Skipping fit: Invalid log input ({current_val_numeric}) for {inv_name} in {name}.")
                            valid_invariant_set = False
                            break # Stop processing this invariant type
                    elif inv_name == "log|Det|":
                         # raw_value here is log|Det|. Check it's valid (> -inf)
                         if current_val_numeric > -np.inf:
                            processed_value = current_val_numeric
                         else:
                            print(f"  Skipping fit: Invalid log input ({current_val_numeric}) for {inv_name} in {name}.")
                            valid_invariant_set = False
                            break # Stop processing this invariant type
                    else:
                        processed_value = current_val_numeric # Default case

                    # Final check if processing resulted in None unexpectedly
                    if processed_value is None:
                        print(f"  Skipping fit: Processed value became None for {inv_name} in {name}.")
                        valid_invariant_set = False
                        break # Stop processing this invariant type

                    x_values.append(float(processed_value)) # Ensure float for numpy array
                    # --- End processing single value ---

                except Exception as e_proc:
                    # Catch any unexpected error during processing THIS knot's value
                    print(f"  Skipping fit: Error processing value for {inv_name} in {name}: {e_proc}")
                    valid_invariant_set = False
                    break # Stop processing this invariant type
            # --- End loop through knots for the CURRENT invariant ---

            # Check if the loop completed successfully for all knots for this invariant
            if not valid_invariant_set or len(x_values) != 3:
                if valid_invariant_set: # Only print if loop finished but count is wrong
                     print(f"  Fit cannot be performed: Expected 3 valid data points, got {len(x_values)}.")
                # else: (reason for skipping was already printed)
                continue # Move to the next invariant

            # Perform the fit (Only if all 3 points were valid)
            x_array = np.array(x_values)
            y_array = y_values # Already an array
            try:
                coeffs = np.polyfit(x_array, y_array, 1)
                A_fit, B_fit = coeffs[0], coeffs[1]
                predictions = A_fit * x_array + B_fit
                # Avoid division by zero if all y values are the same
                denom = np.sum((y_array - np.mean(y_array))**2)
                if denom < 1e-12:
                    r_squared = 1.0 if np.sum((y_array - predictions)**2) < 1e-12 else 0.0
                else:
                    residuals = y_array - predictions
                    r_squared = 1 - np.sum(residuals**2) / denom

                is_offset_zero = abs(B_fit) < offset_threshold

                print(f"  Fit c1 ≈ {A_fit:.3f} * ({inv_name}) + {B_fit:.3f}")
                print(f"    Predictions: e={predictions[0]:.3f}, mu={predictions[1]:.3f}, tau={predictions[2]:.3f}")
                print(f"    Targets    : e={y_array[0]:.3f}, mu={y_array[1]:.3f}, tau={y_array[2]:.3f}")
                print(f"    Goodness of Fit (R²): {r_squared:.3f}")
                print(f"    Intercept near zero (abs < {offset_threshold}): {is_offset_zero}")

            except Exception as e_fit:
                print(f"  Fit failed for {inv_name}: {e_fit}")

    # --- Subsection 4.2: Fit c1 ~ -log|J(q)| based on DYNAMIC Assignment ---
    print("\n--- 4.2: Correlation Test (Dynamic Assignment vs -log|J(q)|) ---")
    # Use knots assigned in Phase 2: final_assigned_knots dictionary
    assigned_leptons = list(final_assigned_knots.keys()) # e.g., ['electron', 'muon', 'tau']
    assigned_knot_names = list(final_assigned_knots.values()) # e.g., ['3_1', '5_1', '5_2']

    log_fit_data_valid = True
    target_c1_dynamic = []
    log_jones_w3_dynamic = []
    log_jones_w5_dynamic = []

    # Check if we have 3 assigned knots
    if len(assigned_leptons) != 3 or any(name == 'N/A' for name in assigned_knot_names):
        print("ERROR: Need exactly 3 assigned knots from Phase 2 for dynamic fit. Skipping.")
        log_fit_data_valid = False
    else:
        # Gather data for assigned knots
        for i, lepton in enumerate(assigned_leptons):
            knot_name = assigned_knot_names[i]
            req_c1 = C1_REQUIRED_TARGETS.get(lepton)
            knot_data = knot_invariant_data.get(knot_name)

            if req_c1 is None:
                print(f"ERROR: Missing required c1 for assigned {lepton} ({knot_name}). Skipping fit.")
                log_fit_data_valid = False; break
            if knot_data is None:
                print(f"ERROR: Missing invariant data for assigned knot {knot_name}. Skipping fit.")
                log_fit_data_valid = False; break

            target_c1_dynamic.append(req_c1)
            if knot_data['log_abs_jones_w3'] is None:
                 print(f"Warning: Missing log|J(ω3)| for assigned knot {knot_name}. Cannot perform ω3 fit.")
                 # Set flag or handle incomplete data during fitting
            else: log_jones_w3_dynamic.append(knot_data['log_abs_jones_w3'])
            if knot_data['log_abs_jones_w5'] is None:
                 print(f"Warning: Missing log|J(ω5)| for assigned knot {knot_name}. Cannot perform ω5 fit.")
            else: log_jones_w5_dynamic.append(knot_data['log_abs_jones_w5'])

    # Perform fits if data is valid and complete
    r_squared_w3_dyn = None # Initialize for status report
    r_squared_w5_dyn = None

    if not log_fit_data_valid:
        print("Skipping dynamic correlation fits due to missing data.")
    else:
        target_c1_array_dyn = np.array(target_c1_dynamic)
        print(f"Using dynamically assigned knots: {final_assigned_knots}")
        print(f"Target c1 values: {target_c1_array_dyn}")

        # Fit vs -log|J(ω3)|
        if len(log_jones_w3_dynamic) == 3:
            print("\nAttempting fit: c1 = A * (-log|J(ω3)|) + B (Dynamic Assignment)")
            x_ljw3_dyn = -np.array(log_jones_w3_dynamic) # Use negative log
            y_c1_dyn = target_c1_array_dyn
            try:
                coeffs_w3_dyn = np.polyfit(x_ljw3_dyn, y_c1_dyn, 1) # Linear fit
                A_w3_dyn, B_w3_dyn = coeffs_w3_dyn[0], coeffs_w3_dyn[1]
                print(f"  Fit c1 ≈ {A_w3_dyn:.3f} * (-log|J(ω3)|) + {B_w3_dyn:.3f}")
                c1_pred_w3_dyn = A_w3_dyn * x_ljw3_dyn + B_w3_dyn
                print(f"    Predictions: e={c1_pred_w3_dyn[0]:.3f}, mu={c1_pred_w3_dyn[1]:.3f}, tau={c1_pred_w3_dyn[2]:.3f}")
                # Calculate R^2
                residuals_w3_dyn = y_c1_dyn - c1_pred_w3_dyn
                r_squared_w3_dyn = 1 - np.sum(residuals_w3_dyn**2) / np.sum((y_c1_dyn - np.mean(y_c1_dyn))**2)
                print(f"    Goodness of Fit (R²): {r_squared_w3_dyn:.3f}")
            except Exception as e: print(f"  Fit failed for log|J(ω3)|: {e}")
        else: print("\nSkipping fit vs log|J(ω3)| due to incomplete data for assigned knots.")

        # Fit vs -log|J(ω5)|
        if len(log_jones_w5_dynamic) == 3:
            print("\nAttempting fit: c1 = A * (-log|J(ω5)|) + B (Dynamic Assignment)")
            x_ljw5_dyn = -np.array(log_jones_w5_dynamic) # Use negative log
            y_c1_dyn = target_c1_array_dyn
            try:
                coeffs_w5_dyn = np.polyfit(x_ljw5_dyn, y_c1_dyn, 1) # Linear fit
                A_w5_dyn, B_w5_dyn = coeffs_w5_dyn[0], coeffs_w5_dyn[1]
                print(f"  Fit c1 ≈ {A_w5_dyn:.3f} * (-log|J(ω5)|) + {B_w5_dyn:.3f}")
                c1_pred_w5_dyn = A_w5_dyn * x_ljw5_dyn + B_w5_dyn
                print(f"    Predictions: e={c1_pred_w5_dyn[0]:.3f}, mu={c1_pred_w5_dyn[1]:.3f}, tau={c1_pred_w5_dyn[2]:.3f}")
                # Calculate R^2
                residuals_w5_dyn = y_c1_dyn - c1_pred_w5_dyn
                r_squared_w5_dyn = 1 - np.sum(residuals_w5_dyn**2) / np.sum((y_c1_dyn - np.mean(y_c1_dyn))**2)
                print(f"    Goodness of Fit (R²): {r_squared_w5_dyn:.3f}")
            except Exception as e: print(f"  Fit failed for log|J(ω5)|: {e}")
        else: print("\nSkipping fit vs log|J(ω5)| due to incomplete data for assigned knots.")

    # --- Subsection 4.3: Test Refined Hypothesis (Ansatz 3e: c1 ~ -alpha * log|J(q)|) ---
    print("\n--- 4.3: Correlation Test (Refined Hypothesis: c1 ≈ C_F * [-alpha*log|J(ω5)|]) ---")

    if not log_fit_data_valid:
        print("Skipping refined hypothesis test due to missing dynamic fit data.")
        refined_fit_r_squared = None
        refined_fit_slope = None
        refined_fit_offset = None
    elif len(log_jones_w5_dynamic) != 3:
        print("Skipping refined hypothesis test due to incomplete log|J(ω5)| data for assigned knots.")
        refined_fit_r_squared = None
        refined_fit_slope = None
        refined_fit_offset = None
    else:
        print(f"Using dynamically assigned knots: {final_assigned_knots}")
        # Independent variable X = -alpha * log|J(ω5)|
        x_refined = -alpha_fine_structure * np.array(log_jones_w5_dynamic)
        y_c1_dyn = target_c1_array_dyn # Dependent variable is still required c1

        print(f"Independent variable X = -alpha * log|J(ω5)| values: {x_refined}")
        print(f"Dependent variable Y = Required c1 values: {y_c1_dyn}")

        print("\nAttempting fit: c1 = C_F * X + Offset")
        try:
            coeffs_refined = np.polyfit(x_refined, y_c1_dyn, 1) # Linear fit
            refined_fit_slope, refined_fit_offset = coeffs_refined[0], coeffs_refined[1] # C_F = slope
            print(f"  Fit c1 ≈ {refined_fit_slope:.3f} * X + {refined_fit_offset:.3f}")
            c1_pred_refined = refined_fit_slope * x_refined + refined_fit_offset
            print(f"    Predictions: e={c1_pred_refined[0]:.3f}, mu={c1_pred_refined[1]:.3f}, tau={c1_pred_refined[2]:.3f}")
            print(f"    Targets    : e={y_c1_dyn[0]:.3f}, mu={y_c1_dyn[1]:.3f}, tau={y_c1_dyn[2]:.3f}")
            # Calculate R^2
            residuals_refined = y_c1_dyn - c1_pred_refined
            refined_fit_r_squared = 1 - np.sum(residuals_refined**2) / np.sum((y_c1_dyn - np.mean(y_c1_dyn))**2)
            print(f"    Goodness of Fit (R²): {refined_fit_r_squared:.3f}")
            # Check if offset is near zero
            offset_threshold = 0.05 # Arbitrary threshold for "near zero"
            is_offset_zero = abs(refined_fit_offset) < offset_threshold
            print(f"    Intercept (Offset) near zero (abs < {offset_threshold}): {is_offset_zero}")
        except Exception as e:
             print(f"  Fit failed for refined hypothesis: {e}")
             refined_fit_r_squared = None
             refined_fit_slope = None
             refined_fit_offset = None

    print("\n--- Correlation Search Finished ---")
    print("NOTE: Correlation fits are speculative numerology without theoretical derivation.")
    print("      More sophisticated fitting or theoretical work is needed.")

    # --- FINAL MODEL STATUS Printout (incorporating results) ---
    print("\n" + "="*20 + " FINAL MODEL STATUS " + "="*20)
    # Use a version number consistent with changes
    print("(LSC Comprehensive Particle Model - v15 w/ Hybrid ID & Refined Log|J| Correlation Test):")
    # Final assignment based on braid search complexity ordering
    electron_knot = final_assigned_knots.get('electron','N/A')
    muon_knot = final_assigned_knots.get('muon','N/A')
    tau_knot = final_assigned_knots.get('tau','N/A')
    # Ensure final_consistent is defined (might not be if braid search failed early)
    final_consistent = final_consistent if 'final_consistent' in locals() else False

    print(f"  - Lepton Hierarchy (Braid Search): Assigning simplest chiral knots by Nc ")
    print(f"     ({electron_knot}(e), {muon_knot}(mu), {tau_knot}(tau))")
    print(f"     Consistent Complexity Ordering Found: {final_consistent}")
    # Print required c1 values for this assignment
    c1_e = f"c1({electron_knot})≈{C1_REQUIRED_TARGETS.get('electron',0):.3f}" if electron_knot != 'N/A' else "c1(?)≈N/A"
    c1_mu = f"c1({muon_knot})≈{C1_REQUIRED_TARGETS.get('muon',0):.3f}" if muon_knot != 'N/A' else "c1(?)≈N/A"
    c1_tau = f"c1({tau_knot})≈{C1_REQUIRED_TARGETS.get('tau',0):.3f}" if tau_knot != 'N/A' else "c1(?)≈N/A"
    print(f"     Required c1 values: {c1_e}, {c1_mu}, {c1_tau}")

    # Report on Correlation Search based on FIXED hypothesis (3_1, 5_1, 5_2)
    print(f"  - Correlation Search (Fixed Hypothesis: 3_1(e), 5_1(mu), 5_2(tau)):")
    # Check if fits were attempted
    if fixed_fit_attempted: # Use the flag
        print(f"     Individual linear fits (c1 vs Invariant) performed for various invariants (see Phase 4.1 output).")
        # Specific previous fits removed, summary now refers to the loop output.
    else:
        print("     Fixed Fit attempt skipped due to missing invariant data or c1 targets.")
    # print("     Status: Simple correlations inconclusive.") # Moved status to new section

    # Report on Correlation Search based on DYNAMIC assignment vs -log|J(q)|
    print(f"  - Correlation Search (Dynamic Assignment vs -log|J(q)|):")
    if log_fit_data_valid:
         fit_desc = "Fits performed:"
         if r_squared_w3_dyn is not None: fit_desc += f" vs -log|J(ω3)| (R²={r_squared_w3_dyn:.3f})"
         else: fit_desc += " vs -log|J(ω3)| (Skipped/Failed)"
         if r_squared_w5_dyn is not None: fit_desc += f"; vs -log|J(ω5)| (R²={r_squared_w5_dyn:.3f})"
         else: fit_desc += "; vs -log|J(ω5)| (Skipped/Failed)"
         print(f"     {fit_desc}")
    else:
         print("     Log|J(q)| Fit attempt skipped due to missing assigned knots or invariant data.")
    print(f"     Overall Status: Correlations remain speculative.")

    # Report on Refined Hypothesis Test (Ansatz 3e)
    print(f"  - Correlation Search (Refined Hypothesis c1 ~ C_F*[-alpha*log|J(ω5)|]):")
    if refined_fit_r_squared is not None and refined_fit_slope is not None and refined_fit_offset is not None:
        print(f"     Fit: c1 ≈ {refined_fit_slope:.3f} * X + {refined_fit_offset:.3f} (R²={refined_fit_r_squared:.3f})")
        print(f"     Intercept near zero: {abs(refined_fit_offset) < offset_threshold}")
    else:
        print("     Refined fit skipped or failed.")

    # Placeholder for neutrino C_g value if not defined elsewhere
    neutrino_Cg_fit = 0.0 # Needs to be calculated/fitted elsewhere if desired
    print(f"  - Neutral Neutrino: Requires k ~ C_g*(E0/Ep)^2 (Needs C_g≈{neutrino_Cg_fit:.3f}).")

    print("  - CORE THEORETICAL CHALLENGES:")
    print("      1. Derive c1(Knot Topology) function yielding correct values for assigned chiral knots (from braid search).")
    print("         -> Test relationship with invariants like |Signature|, Det, or log|J(q)|.")
    print("         -> Test refined hypothesis c1 ≈ C_F * [-alpha*log|J(ω5)|]. If promising, derive C_F.")
    print("      2. Test FIXED c1 correlation hypothesis (3_1, 5_1, 5_2) with more sophisticated invariants/fits.") # Kept for reference
    print("      3. Confirm Achiral Knots (e.g., 4_1, 6_2) Map to Neutral Particles (Neutrinos?).")
    print("      4. Derive Gravitational Coupling C_g(Knot Topology) for Achiral Knots.")
    # Add other core challenges back if relevant
    print("      5. Formalize Spinor Phase, Vertex/Propagators, g-2, SU(3)/Confinement, Photon model...")


    # --- Debug Trefoil Detection (kept for verification) ---
    print("\n" + "="*20 + " DEBUGGING TREFOIL " + "="*20)
    # Define a state dictionary for the debug function's calls
    debug_poly_assertion_state = {'printed_count': 0, 'max_prints': 5, 'suppressed_msg_shown': False}
    # Ensure debug_trefoil_detection accepts the state dict argument
    try:
        # Check if function signature accepts argument before calling
        import inspect
        sig = inspect.signature(debug_trefoil_detection)
        if len(sig.parameters) == 1:
            debug_trefoil_detection(debug_poly_assertion_state) # Pass the state
        else:
             debug_trefoil_detection() # Call without if it doesn't accept
    except NameError:
         print("Debug function 'debug_trefoil_detection' not defined.") # Should not happen now
    except Exception as e_debug:
         print(f"Error running debug function: {e_debug}")

    # --- Subsection 4.4: Multi-Variable Fit Example (Dynamic Assignment) ---
    print("\n--- 4.4: Multi-Variable Fit Test (c1 vs log|Det| and |Signature|) ---")

    if not SKLEARN_AVAILABLE:
        print("  Skipping multi-variable fit: sklearn library not available.")
    elif not log_fit_data_valid:
        print("  Skipping multi-variable fit due to missing dynamic assignment data.")
    else:
        mv_model = None # Define model variable outside try block for status check
        print(f"Using dynamically assigned knots: {final_assigned_knots}")
        multi_var_data_valid = True
        X_multi_var = [] # List of lists for features
        y_multi_var = [] # Target variable (c1)

        # Define features to include (keys from knot_invariant_data)
        # Use log|Det| and |Signature|
        feature_keys = ['log_abs_determinant', 'signature']
        feature_names = ['log|Det|', '|Signature|'] # Names for reporting

        for i, lepton in enumerate(assigned_leptons):
            knot_name = assigned_knot_names[i]
            knot_data = knot_invariant_data.get(knot_name)
            req_c1 = C1_REQUIRED_TARGETS.get(lepton)

            if req_c1 is None or knot_data is None:
                multi_var_data_valid = False; break # Should be caught by log_fit_data_valid check

            features_for_knot = []
            for k, key in enumerate(feature_keys):
                raw_value = knot_data.get(key)
                if raw_value is None:
                    print(f"  Skipping multi-var fit: Missing feature '{key}' for {knot_name}.")
                    multi_var_data_valid = False; break
                # Apply transformations (ensure numeric first)
                try:
                    numeric_val = float(raw_value) # Assume log/abs done in calculate_invariants
                    if key == 'log_abs_determinant':
                        features_for_knot.append(numeric_val) # Already log
                    elif key == 'signature':
                        features_for_knot.append(abs(numeric_val))
                    else:
                        features_for_knot.append(numeric_val) # Add other transformations if needed
                except Exception as e:
                    print(f"  Skipping multi-var fit: Error processing feature '{key}' for {knot_name}: {e}")
                    multi_var_data_valid = False; break
            if not multi_var_data_valid: break

            X_multi_var.append(features_for_knot)
            y_multi_var.append(req_c1)

        if multi_var_data_valid and len(X_multi_var) == 3:
            X_multi_array = np.array(X_multi_var)
            y_multi_array = np.array(y_multi_var)

            try:
                mv_model = LinearRegression()
                mv_model.fit(X_multi_array, y_multi_array)
                coeffs_multi = mv_model.coef_
                intercept_multi = mv_model.intercept_
                r_squared_multi = mv_model.score(X_multi_array, y_multi_array)

                coeff_str = " + ".join([f"{c:.3f}*{fn}" for c, fn in zip(coeffs_multi, feature_names)])
                print(f"  Multi-Var Fit: c1 ≈ {coeff_str} + {intercept_multi:.3f}")
                print(f"    Goodness of Fit (R²): {r_squared_multi:.4f}")
                print(f"    Intercept: {intercept_multi:.4f}")
                is_mv_offset_zero = abs(intercept_multi) < offset_threshold
                print(f"    Intercept near zero (abs < {offset_threshold}): {is_mv_offset_zero}")

            except Exception as e_mv_fit:
                print(f"  Multi-variable fit failed: {e_mv_fit}")
        elif multi_var_data_valid:
            print(f"  Skipping multi-var fit: Expected 3 data points, got {len(X_multi_var)}. Please check input data.")
        # else: reason for skipping already printed

    print("\n--- Correlation Search Finished ---")
    print("NOTE: Correlation fits are speculative numerology without theoretical derivation.")
    print("      More sophisticated fitting or theoretical work is needed.")

    # --- FINAL MODEL STATUS Printout (incorporating results) --- (Phase 5)
    print("\n" + "="*20 + " PHASE 5: FINAL MODEL STATUS " + "="*20)
    # Use a version number consistent with changes
    print("(LSC Comprehensive Particle Model - v17 w/ Refactored Fits & Multi-Var Test):")
    # Final assignment based on braid search complexity ordering
    electron_knot = final_assigned_knots.get('electron','N/A')
    muon_knot = final_assigned_knots.get('muon','N/A')
    tau_knot = final_assigned_knots.get('tau','N/A')
    # Ensure final_consistent is defined (might not be if braid search failed early)
    final_consistent = final_consistent if 'final_consistent' in locals() else False

    print(f"  - Lepton Hierarchy (Braid Search): Assigning simplest chiral knots by Nc ")
    print(f"     ({electron_knot}(e), {muon_knot}(mu), {tau_knot}(tau))")
    print(f"     Consistent Complexity Ordering Found: {final_consistent}")
    # Print required c1 values for this assignment
    c1_e = f"c1({electron_knot})≈{C1_REQUIRED_TARGETS.get('electron',0):.3f}" if electron_knot != 'N/A' else "c1(?)≈N/A"
    c1_mu = f"c1({muon_knot})≈{C1_REQUIRED_TARGETS.get('muon',0):.3f}" if muon_knot != 'N/A' else "c1(?)≈N/A"
    c1_tau = f"c1({tau_knot})≈{C1_REQUIRED_TARGETS.get('tau',0):.3f}" if tau_knot != 'N/A' else "c1(?)≈N/A"
    print(f"     Required c1 values: {c1_e}, {c1_mu}, {c1_tau}")

    # Report on Correlation Search based on FIXED hypothesis (3_1, 5_1, 5_2)
    print(f"  - Correlation Search (Fixed Hypothesis: 3_1(e), 5_1(mu), 5_2(tau)):")
    if fixed_fit_attempted: # Use the flag
        print(f"     Individual linear fits (c1 vs Invariant) performed for various invariants (see Phase 4.1 output).")
        # Specific previous fits removed, summary now refers to the loop output.
    else:
        print("     Fixed Fit attempt skipped due to missing invariant data or c1 targets.")
    # print("     Status: Simple correlations inconclusive.") # Moved status to new section

    # Report on Correlation Search based on DYNAMIC assignment vs -log|J(q)|
    print(f"  - Correlation Search (Dynamic Assignment vs -log|J(q)|):")
    if log_fit_data_valid:
         fit_desc = "Fits performed:"
         if r_squared_w3_dyn is not None: fit_desc += f" vs -log|J(ω3)| (R²={r_squared_w3_dyn:.3f})"
         else: fit_desc += " vs -log|J(ω3)| (Skipped/Failed)"
         if r_squared_w5_dyn is not None: fit_desc += f"; vs -log|J(ω5)| (R²={r_squared_w5_dyn:.3f})"
         else: fit_desc += "; vs -log|J(ω5)| (Skipped/Failed)"
         print(f"     {fit_desc}")
    else:
         print("     Log|J(q)| Fit attempt skipped due to missing assigned knots or invariant data.")
    print(f"     Overall Status: Correlations remain speculative.")

    # Report on Refined Hypothesis Test (Ansatz 3e)
    print(f"  - Correlation Search (Refined Hypothesis c1 ~ C_F*[-alpha*log|J(ω5)|]):")
    if refined_fit_r_squared is not None and refined_fit_slope is not None and refined_fit_offset is not None:
        print(f"     Fit: c1 ≈ {refined_fit_slope:.3f} * X + {refined_fit_offset:.3f} (R²={refined_fit_r_squared:.3f})")
        print(f"     Intercept near zero: {abs(refined_fit_offset) < offset_threshold}")
    else:
        print("     Refined fit skipped or failed.")

    # Placeholder for neutrino C_g value if not defined elsewhere
    neutrino_Cg_fit = 0.0 # Needs to be calculated/fitted elsewhere if desired
    print(f"  - Neutral Neutrino: Requires k ~ C_g*(E0/Ep)^2 (Needs C_g≈{neutrino_Cg_fit:.3f}).")

    print("  - CORE THEORETICAL CHALLENGES:")
    print("      1. Derive c1(Knot Topology) function yielding correct values for assigned chiral knots (from braid search).")
    print("         -> Test relationship with invariants like |Signature|, Det, or log|J(q)|.")
    print("         -> Test refined hypothesis c1 ≈ C_F * [-alpha*log|J(ω5)|]. If promising, derive C_F.")
    print("      2. Test FIXED c1 correlation hypothesis (3_1, 5_1, 5_2) with more sophisticated invariants/fits.") # Kept for reference
    print("      3. Confirm Achiral Knots (e.g., 4_1, 6_2) Map to Neutral Particles (Neutrinos?).")
    print("      4. Derive Gravitational Coupling C_g(Knot Topology) for Achiral Knots.")
    # Add other core challenges back if relevant
    print("      5. Formalize Spinor Phase, Vertex/Propagators, g-2, SU(3)/Confinement, Photon model...")

    # Report on Multi-Variable Fit Test
    print(f"  - Correlation Search (Multi-Variable Example Fit):")
    if SKLEARN_AVAILABLE and 'mv_model' in locals() and 'r_squared_multi' in locals(): # Check if fit was attempted and successful
        coeff_str_summary = " + ".join([f"{c:.3f}*{fn}" for c, fn in zip(coeffs_multi, feature_names)])
        print(f"     Fit: c1 ≈ {coeff_str_summary} + {intercept_multi:.3f} (R²={r_squared_multi:.4f})")
        print(f"     Intercept near zero: {is_mv_offset_zero}")
    elif not SKLEARN_AVAILABLE:
        print("     Multi-variable fit skipped: sklearn not available.")
    else:
        print("     Multi-variable fit skipped or failed (check Phase 4.4 output).")

    # Placeholder for neutrino C_g value if not defined elsewhere
    neutrino_Cg_fit = 0.0 # Needs to be calculated/fitted elsewhere if desired
    print(f"  - Neutral Neutrino: Requires k ~ C_g*(E0/Ep)^2 (Needs C_g≈{neutrino_Cg_fit:.3f}).")

    print("  - CORE THEORETICAL CHALLENGES:")
    print("      1. Derive c1(Knot Topology) function yielding correct values for assigned chiral knots (from braid search).")
    print("         -> Test relationship with invariants like |Signature|, Det, or log|J(q)|.")
    print("         -> Test refined hypothesis c1 ≈ C_F * [-alpha*log|J(ω5)|]. If promising, derive C_F.")
    print("      2. Test FIXED c1 correlation hypothesis (3_1, 5_1, 5_2) with more sophisticated invariants/fits.") # Kept for reference
    print("      3. Confirm Achiral Knots (e.g., 4_1, 6_2) Map to Neutral Particles (Neutrinos?).")
    print("      4. Derive Gravitational Coupling C_g(Knot Topology) for Achiral Knots.")
    # Add other core challenges back if relevant
    print("      5. Formalize Spinor Phase, Vertex/Propagators, g-2, SU(3)/Confinement, Photon model...")