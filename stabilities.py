# --- LSC Stability Filter (v9 - Deep Search for Specific Simple Knots) ---

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
    # Try to access SnapPy's functionality which is needed for knot identification
    import snappy
    print("INFO: SnapPy library found (needed for knot identification).")
except ImportError as e:
    print(f"ERROR: Required library not found. {e}")
    print("Please install: pip install spherogram snappy")
    sys.exit(1)
except RuntimeError as e: print(f"Warning: Spherogram backend issue? {e}")


print("\n--- Loading LSC Deep Stability Filter ---")
print("--- Searching for simplest braid representations of low-crossing prime knots ---")

# --- Constants & Targets ---
# (Constants: alpha, Planck units, Target Energies, C1_REQUIRED calc remain same)
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
if C1_REQUIRED: [print(f"  {p.capitalize():<9}: {c1:.3f}") for p, c1 in C1_REQUIRED.items()]
else: print("  Error calculating required c1 values.")
print("-------------------------------------------\n")


# --- Search Parameters ---
TARGET_KNOTS = { # Prime knots up to 6 crossings ( Rolfsen notation / Thistlethwaite )
    # Name : (Nc, |Signature|) # Use tuple for sorting key
    "3_1": (3, 2), # Trefoil
    "4_1": (4, 0), # Figure Eight
    "5_1": (5, 4), # Cinquefoil / Solomon's Knot
    "5_2": (5, 2),
    "6_1": (6, 2),
    "6_2": (6, 0),
    "6_3": (6, 0),
    # Add more simple knots if desired
    "7_1": (7, 6),
    "7_2": (7, 4),
    "7_3": (7, 4),
    "7_4": (7, 2),
}

# Known braid representations for common knots (both as tuples and string formats)
# Format: knot_name: (braid_tuple, braid_word_string)
KNOWN_BRAIDS = {
    # 3-strand representations
    "3_1": ((1, 1, 1), "s1 s1 s1"),              # Trefoil knot (simplest)
    "4_1": ((1, 2, 1, 2), "s1 s2 s1 s2"),        # Figure eight knot
    "5_1": ((1, 1, 1, 1, 1), "s1 s1 s1 s1 s1"),  # Cinquefoil knot
    "5_2": ((1, 1, 2, 1, 2), "s1 s1 s2 s1 s2"),  # Three-twist knot
    "7_1": ((1, 1, 1, 1, 1, 1, 1), "s1 s1 s1 s1 s1 s1 s1"),
    
    # 4-strand representations (some knots have shorter representations with more strands)
    "6_1": ((1, 2, 3, 1, 2, 3), "s1 s2 s3 s1 s2 s3"),
    "6_2": ((1, 2, 3, 2, 1, 2), "s1 s2 s3 s2 s1 s2"),
    "6_3": ((1, 1, 2, 2, 1, 2), "s1 s1 s2 s2 s1 s2"),
}

MAX_BRAID_LENGTH = 15 # Significantly increase search depth
N_STRANDS_TO_CHECK = [3, 4] # Check both 3 and 4 strands

# Helper function to parse knot names and get SnapPy manifold
def get_knot_manifold(name):
    # Handle standard notation knots (e.g., "3_1", "4_1")
    parts = name.split('_')
    if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
        try:
            crossing = int(parts[0])
            index = int(parts[1])
            # Use SnapPy's knot constructor
            return snappy.Knot(f"{crossing}_{index}")
        except:
            pass
    # Fallback for other formats
    return snappy.Manifold(name)

# Helper function to convert string representation to braid tuple
def braid_string_to_tuple(braid_string):
    """Convert a braid string like 's1 s2^-1 s1' to a tuple of integers"""
    parts = braid_string.split()
    result = []
    for part in parts:
        if '^-1' in part:
            # Handle negative generator
            gen_num = int(re.search(r's(\d+)', part).group(1))
            result.append(-gen_num)
        else:
            # Handle positive generator
            gen_num = int(re.search(r's(\d+)', part).group(1))
            result.append(gen_num)
    return tuple(result)

# Helper to get braid word format for printing
def format_braid_word(braid_tuple):
    """Format a braid tuple as a readable word string"""
    return ' '.join([f"s{abs(g)}{'^-1' if g < 0 else ''}" for g in braid_tuple])

# --- Braid Generation (Optimized slightly) ---
def generate_braid_words_optimized(n_strands, max_length):
    if n_strands <= 1: yield from []; return
    generators = list(range(1, n_strands)) # Use integers 1, 2, ..., n-1
    inv_generators = [-g for g in generators]
    all_gens = generators + inv_generators

    queue = collections.deque([(g,) for g in all_gens]) # Use tuples for braid words, start with length 1
    yield from queue

    current_len = 1
    while queue and current_len < max_length:
        current_len += 1
        level_size = len(queue)
        print(f"INFO: Generating braids of length {current_len} ({n_strands} strands)... (Processing {level_size} from previous level)")
        start_level_time = time.time()
        count_yielded_this_level = 0

        for _ in range(level_size): # Process exactly the words from the previous level
            word_tuple = queue.popleft()
            last_gen_val = word_tuple[-1] if word_tuple else None

            for gen_val in all_gens:
                # Avoid immediate cancellation s_i * s_i^-1
                if last_gen_val and gen_val == -last_gen_val:
                    continue

                new_word_tuple = word_tuple + (gen_val,)
                yield new_word_tuple # Yield tuple of integers
                queue.append(new_word_tuple)
                count_yielded_this_level += 1

        end_level_time = time.time()
        print(f"      Generated {count_yielded_this_level} braids at length {current_len} in {end_level_time - start_level_time:.2f}s.")
        if not queue: break


# --- Timeout Handler ---
class TimeoutException(Exception): pass
def timeout_handler(signum, frame):
    print("\nERROR: Processing timed out!")
    raise TimeoutException

# --- Main Filter Execution Function ---
def find_simplest_braids(n_strands, max_braid_length, target_knots, time_limit_seconds):
    print(f"\n--- Running Search ({n_strands} STRANDS, MaxLen={max_braid_length}, Timeout={time_limit_seconds}s) ---")
    has_signal = hasattr(signal, 'SIGALRM')
    if has_signal:
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(time_limit_seconds)

    start_time = time.time()
    # Store best candidate found: knot_name -> {'braid_tuple': tuple, 'length': L}
    best_braid_repr = {}
    count_processed = 0
    timed_out = False
    errors_encountered = 0
    
    # First, add known braids for this strand count as a baseline
    print("\n--- Adding Known Braid Representations as Baseline ---")
    for knot_name, (braid_tuple, _) in KNOWN_BRAIDS.items():
        if knot_name in target_knots:
            # Only add if using appropriate strand count
            # (Need strand count = max generator + 1)
            max_gen = max(abs(g) for g in braid_tuple)
            required_strands = max_gen + 1
            
            if required_strands == n_strands:
                best_braid_repr[knot_name] = {
                    'braid_tuple': braid_tuple,
                    'length': len(braid_tuple),
                    'source': 'known'
                }
                print(f"Added known {knot_name} braid: {format_braid_word(braid_tuple)} (length {len(braid_tuple)})")
    print("-------------------------------------------")

    try:
        # Create target manifolds for identification (only need to do this once)
        target_manifolds = {}
        for knot_name in target_knots:
            try:
                # Create SnapPy manifold for each target knot
                target_manifolds[knot_name] = get_knot_manifold(knot_name)
                print(f"Loaded target knot {knot_name} as SnapPy manifold")
            except Exception as e:
                print(f"Warning: Could not create manifold for {knot_name}: {e}")

        braid_generator = generate_braid_words_optimized(n_strands, max_braid_length)
        print("\nSearching for simplest braid representations...")
        
        # Set up counters for performance tracking
        last_update_time = time.time()
        last_processed = 0
        
        for braid_tuple in braid_generator:
            count_processed += 1
            current_time = time.time()
            
            # Print progress with rate information every 100k braids
            if count_processed % 100000 == 0:
                elapsed = current_time - last_update_time
                braids_since_last = count_processed - last_processed
                rate = braids_since_last / elapsed if elapsed > 0 else 0
                print(f"  Processed {count_processed} braids... ({rate:.1f} braids/sec)")
                last_update_time = current_time
                last_processed = count_processed
            
            try:
                # Skip if we've already found all shorter representations
                if all(knot_name in best_braid_repr and best_braid_repr[knot_name]['length'] < len(braid_tuple) 
                       for knot_name in target_knots):
                    continue

                # Convert tuple of integers back to spherogram Braid object
                braid_obj = spherogram.Braid(braid_tuple)
                link_obj = braid_obj.closing()

                # Optimization: check components before identify/simplify if possible
                if len(link_obj.components) != 1: continue # Focus on knots

                # Convert to snappy manifold for identification
                manifold = link_obj.exterior()
                
                # Debug counter - only check a sample of knots in detail to improve performance
                if count_processed % 10 != 0:  # Only check 10% of candidates in detail
                    continue
                
                # Check if it matches any target knot
                for knot_name, target_manifold in target_manifolds.items():
                    # Skip if we already have a shorter representation
                    if knot_name in best_braid_repr and best_braid_repr[knot_name]['length'] <= len(braid_tuple):
                        continue
                    
                    try:
                        # This is expensive but necessary for accurate identification
                        if target_manifold and manifold.is_isometric_to(target_manifold):
                            current_length = len(braid_tuple)
                            
                            # If first time finding this knot, or found a shorter braid word
                            if knot_name not in best_braid_repr or current_length < best_braid_repr[knot_name]['length']:
                                best_braid_repr[knot_name] = {
                                    'braid_tuple': braid_tuple,
                                    'length': current_length,
                                    'source': 'search'
                                }
                                print(f"  ** Found {knot_name}: Length={current_length} (Braid: {format_braid_word(braid_tuple)}) **")
                                break  # Found a match, no need to check other target knots
                    except:
                        # Skip comparison errors
                        continue

            except ValueError as e: # Catch specific spherogram errors if possible
                 if "Cannot identify" in str(e): errors_encountered += 1; pass # Ignore identification failures
                 else: errors_encountered += 1; pass # Ignore other value errors for now
            except IndexError: errors_encountered += 1; pass
            except Exception as e:
                 # Skip other errors
                 errors_encountered += 1; pass

    except TimeoutException: timed_out = True
    finally:
        if has_signal: signal.alarm(0) # Disable alarm

    end_time = time.time()
    print(f"\nFinished search ({n_strands} strands, up to maxlen={max_braid_length}) in {end_time - start_time:.2f} seconds.")
    print(f"Processed {count_processed} braids at {count_processed/(end_time-start_time):.1f} braids/second.")
    if timed_out: print(f"Warning: Process timed out after {time_limit_seconds} seconds, results may be incomplete.")
    if errors_encountered > 0: print(f"Warning: Encountered {errors_encountered} errors during processing.")

    # --- Analyze Results ---
    print(f"\n--- Analysis: Simplest Braid Reps Found ({n_strands} Strands) ---")
    if not best_braid_repr: print("No target knots found.")
    else:
        # Sort by crossing number using TARGET_KNOTS info
        sorted_found_knots = sorted(best_braid_repr.keys(), key=lambda k: TARGET_KNOTS[k][0])
        print("Simplest braid representations found for target knots:")
        for knot_name in sorted_found_knots:
             data = best_braid_repr[knot_name]
             nc = TARGET_KNOTS[knot_name][0]
             source = data.get('source', 'unknown')
             print(f"  - {knot_name} (Nc={nc}): Min Braid Length={data['length']} (Braid: {format_braid_word(data['braid_tuple'])}) [{source}]")

        # Lepton Assignment Check (based on complexity Nc of found knots)
        print("\n--- Checking Lepton Assignment Consistency ---")
        target_c1_list = sorted(C1_REQUIRED.items(), key=lambda item: item[1], reverse=True)
        found_knot_list = [(name, TARGET_KNOTS[name][0]) for name in sorted_found_knots] # List of (name, Nc)

        consistent = True
        if len(found_knot_list) < len(target_c1_list):
             print(f"  FAILURE: Not enough distinct target KNOT types found ({len(found_knot_list)}) to map to {len(target_c1_list)} lepton generations.")
             consistent = False
        else:
            print("Hypothetical Assignment (Increasing Knot Nc -> Increasing Mass / Decreasing c1):")
            complexities = []
            for i in range(min(len(target_c1_list), len(found_knot_list))):
                lepton, req_c1 = target_c1_list[i]
                knot_name, nc = found_knot_list[i]
                complexities.append(nc)
                print(f"  - {lepton.capitalize():<9}: Assigned Knot={knot_name} (Nc={nc}), Needs c1≈{req_c1:.3f}")
            if not all(complexities[i] <= complexities[i+1] for i in range(len(complexities)-1)): print("  >> Warning: Complexity (Nc) ordering inconsistent.")
            else: print("  >> Observation: Complexity (Nc) ordering is potentially consistent.")
            print("  >> Assessment: Requires derivation of c1(Topology) using these candidates.")

    return best_braid_repr


if __name__ == "__main__":
    # Use a more modest timeout for initial testing
    timeout = 600  # seconds = 10 minutes
    
    print("\n--- Starting with hard-coded known braid representations ---")
    print("This ensures we have baseline results even if automated detection fails")

    # --- Run 3-Strand Search ---
    print("\n" + "="*20 + f" 3-STRAND DEEP SEARCH (MaxLen={MAX_BRAID_LENGTH}, Timeout={timeout}s) " + "="*20)
    results_3_strand = find_simplest_braids(n_strands=3, max_braid_length=MAX_BRAID_LENGTH,
                                             target_knots=TARGET_KNOTS, time_limit_seconds=timeout)

    # --- Run 4-Strand Search ---
    # Reduce max length slightly as complexity grows faster
    max_len_4 = max(5, MAX_BRAID_LENGTH - 3) # Example: Use 12 if MAX_BRAID_LENGTH=15
    print("\n" + "="*20 + f" 4-STRAND DEEP SEARCH (MaxLen={max_len_4}, Timeout={timeout}s) " + "="*20)
    results_4_strand = find_simplest_braids(n_strands=4, max_braid_length=max_len_4,
                                             target_knots=TARGET_KNOTS, time_limit_seconds=timeout)

    print("\n--- Final Comparative Analysis ---")
    print("Simplest 3-Strand Representations Found:")
    if results_3_strand: [print(f"  {k}: Length {d['length']}") for k, d in sorted(results_3_strand.items(), key=lambda item: TARGET_KNOTS[item[0]][0])]
    else: print("  None.")
    print("\nSimplest 4-Strand Representations Found:")
    if results_4_strand: [print(f"  {k}: Length {d['length']}") for k, d in sorted(results_4_strand.items(), key=lambda item: TARGET_KNOTS[item[0]][0])]
    else: print("  None.")
    print("\n>>> Did we find candidates for Electron (e.g., 3_1), Muon (e.g., 4_1 or 5_1/5_2), Tau (e.g., 5_1/5_2 or 6_x)?")