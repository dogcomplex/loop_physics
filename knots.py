# --- knot_properties.py ---
# Module for defining the Knot object and calculating its invariants

import numpy as np
import math
import sys
import collections
import warnings # Import warnings module

# Suppress specific DeprecationWarning from spherogram if desired
# warnings.filterwarnings("ignore", category=DeprecationWarning, module="spherogram.links.invariants")

try:
    import spherogram
    import snappy # Often needed by spherogram for advanced features
    print("INFO: Spherogram & SnapPy libraries found.")
    # Check for Sage environment more robustly
    try:
        import sage.all
        # Test a simple Sage operation
        _ = sage.all.Integer(1) + sage.all.Integer(1)
        print("INFO: Running inside a functional SageMath environment.")
        SAGE_AVAILABLE = True
        # Define polynomial variables for Sage
        R = sage.all.LaurentPolynomialRing(sage.all.QQ, 'q_var')
        q_var = R.gen()
        T = sage.all.PolynomialRing(sage.all.ZZ, 't')
        t = T.gen()
        # Define roots of unity for evaluation
        ROOTS_OF_UNITY = {
            'w3': sage.all.exp(2 * sage.all.pi * sage.all.I / 3),
            'w4': sage.all.I, # exp(2*pi*i/4) = i
            'w5': sage.all.exp(2 * sage.all.pi * sage.all.I / 5),
            'w6': sage.all.exp(2 * sage.all.pi * sage.all.I / 6),
        }
    except (ImportError, NameError, TypeError) as e:
        print(f"WARNING: SageMath environment not fully functional or detected ({type(e).__name__}). Polynomial calculations will be limited.")
        SAGE_AVAILABLE = False
        ROOTS_OF_UNITY = {} # Ensure it's defined
        t = None # Ensure t is defined
except ImportError:
    print("CRITICAL WARNING: Spherogram library not found. Cannot perform knot calculations.")
    SAGE_AVAILABLE = False
    ROOTS_OF_UNITY = {}
    t = None
    # Optionally exit if spherogram is essential
    # sys.exit("Spherogram is required.")
except Exception as e:
    print(f"Warning: Issue during library import or Sage setup: {e}")
    SAGE_AVAILABLE = False
    ROOTS_OF_UNITY = {}
    t = None

class KnotPropertyCalculator:
    """
    Encapsulates a knot or link and calculates its topological invariants.
    Uses spherogram and potentially snappy. Requires Sage for some polynomials.
    """
    def __init__(self, identifier, source_description="Unknown"):
        """
        Initializes Knot object from an identifier (e.g., '3_1', braid tuple).
        identifier: Standard knot name (str), braid word (str), or braid tuple (tuple).
        source_description: String describing how this knot was generated (e.g., 'Known', 'Braid(1,1,1)').
        """
        self.identifier = identifier
        self.source_description = source_description
        self.link_obj = None         # Raw spherogram Link object
        self.knot_obj = None         # Spherogram Knot object (if single component)
        self.manifold = None
        self.properties = {
            # Core Identifiers
            'input_identifier': identifier, # Keep track of original input
            'source_description': source_description,
            'knot_atlas_name': None, # Standard name like 3_1, L10a140 (best guess)
            # Basic Topology
            'components': None,
            'crossing_number_min': None, # Minimal crossing number from identification (unavailable)
            'signature': None,
            'determinant': None,
            'log_abs_determinant': None, # Log of absolute determinant
            'volume': None, # Hyperbolic volume
            # Diagrammatic Properties
            'writhe': None, # Diagrammatic writhe
            'morse_number': None, # Diagrammatic Morse number (optional calc)
            'is_alternating': None,
            # Chirality / Symmetry
            'is_chiral': None, # True if chiral, False if amphichiral (unavailable)
            # Structural Properties
            'is_torus': None, # Is it a torus knot? (Often heuristic based on volume/name)
            'is_fibered': None, # Is it fibered? (From basic checks, may be inaccurate)
            'is_fibered_hfk': None, # Is it fibered? (From HFK calculation, more robust, optional)
            'seifert_genus': None, # Seifert genus (From HFK calculation, optional)
            # Polynomials & Evaluations
            'alexander_poly': None, # Store polynomial object if possible
            'jones_poly': None,     # Store polynomial object if possible
            'alex_at_neg1': None,   # |Alexander(-1)| = Determinant (usually)
            'alex_at_neg2': None,   # Value of Alexander Polynomial at t=-2
            'jones_at_roots': {}, # Dict: {'w3': val, 'w4': val, ...}
            'log_abs_jones_at_roots': {}, # Dict: {'w3': val, 'w4': val, ...}
            # Representations
            'braid_word_repr': None, # Braid word representation (tuple/string)
            # Properties needing more theoretical definition / postulates
            'framing_number': 'UNKNOWN', # How is framing defined/calculated? Placeholder.
            'topological_charge_Q': 'UNKNOWN', # Mechanism needed. Placeholder.
            'dynamic_stiffness_k': 'UNKNOWN', # Requires physics model. Placeholder.
            'effective_inertia_m_eff': 'UNKNOWN', # Requires physics model. Placeholder.
            'stability_heuristic': 'Unknown', # e.g., 'Prime', 'Composite', 'Unknot', 'Link', 'Error'
        }
        self._calculate_properties()

    def _create_link_object(self):
        """Attempts to create a spherogram Link object from the identifier."""
        if self.link_obj: return True
        try:
            # Handle different identifier types
            if isinstance(self.identifier, tuple): # Braid Tuple
                print(f"Info: Creating Link from Braid tuple {self.identifier}")
                braid_list = [int(g) for g in self.identifier] # Ensure integers
                # Create Braid, then get its Link object
                braid_obj = spherogram.ClosedBraid(braid_list)
                # ClosedBraid is a subclass of Link, so assign it directly
                self.link_obj = braid_obj
                # Store the Knot object too if it's a knot
                try:
                    if self.link_obj.num_components() == 1:
                        self.knot_obj = braid_obj # Assign the braid obj directly as the knot obj
                    else:
                        self.knot_obj = None
                except Exception as e:
                    print(f"Warn: Error checking components/assigning knot_obj for braid {self.identifier}: {e}")
                    self.knot_obj = None
                self.properties['knot_atlas_name'] = f"Braid{self.identifier}" # Temporary name
                self.properties['source_description'] = f"Braid{self.identifier}"

            elif isinstance(self.identifier, str): # String Identifier (Name)
                # Special case for unknot
                if self.identifier == '0_1':
                    print("Info: Creating Unknot via Link([])")
                    self.link_obj = spherogram.Link([])
                    self.properties['knot_atlas_name'] = '0_1'
                    # Explicitly set knot_obj to None for unknot? Or create Knot([])?
                    # Let's assume Link([]) is sufficient for now.
                else:
                    # spherogram.Knot() is unavailable in this environment.
                    # Start with spherogram.Link()
                    try:
                        # Try spherogram.Link()
                        print(f"Info: Attempting spherogram.Link('{self.identifier}')")
                        self.link_obj = spherogram.Link(self.identifier)
                        print(f"Info: Successfully created Link object for '{self.identifier}' via Link().")
                        self.properties['knot_atlas_name'] = self.identifier
                        # Since we used Link(), knot_obj is likely None unless Link() returns a Knot subclass?
                        # We'll try the conversion later anyway.
                        self.knot_obj = None
                    except (ValueError, KeyError, Exception) as e_link:
                        print(f"Info: spherogram.Link('{self.identifier}') failed ({type(e_link).__name__}: {e_link}). Falling back to knot_db...")
                        try:
                            # Fallback: Try knot_db
                            print(f"Info: Attempting spherogram.knot_db['{self.identifier}']")
                            k_obj_db = spherogram.knot_db[self.identifier]
                            # Ensure we have a Link object for consistency
                            if isinstance(k_obj_db, spherogram.Link):
                                self.link_obj = k_obj_db
                                self.knot_obj = None # Might not be a knot
                            elif hasattr(k_obj_db, 'link'): # Check if it behaves like a Knot (has link method)
                                self.link_obj = k_obj_db.link()
                                # Can we assign it to knot_obj? Check if it has expected methods.
                                if hasattr(k_obj_db, 'identify'): # Heuristic check
                                    self.knot_obj = k_obj_db # Assume it's a knot object
                                else:
                                    self.knot_obj = None
                            else:
                                raise TypeError("Object from knot_db is not Link or Knot-like.")
                            self.properties['knot_atlas_name'] = self.identifier # Name is confirmed by db lookup
                            print(f"Info: Successfully loaded '{self.identifier}' from knot_db.")
                        except (KeyError, TypeError, AttributeError, Exception) as e_db:
                            print(f"ERROR: Identifier '{self.identifier}' not found via Link() or knot_db ({type(e_db).__name__}: {e_db}). Cannot create object.")
                            self.properties['stability_heuristic'] = 'Error'
                            return False
            else:
                print(f"ERROR: Invalid identifier type: {type(self.identifier)}")
                self.properties['stability_heuristic'] = 'Error'
                return False

            # Basic check after creation
            if self.link_obj is None:
                 print(f"ERROR: Link object creation resulted in None for '{self.identifier}'.")
                 self.properties['stability_heuristic'] = 'Error'
                 return False

            return True
        except ImportError:
             print("ERROR: Spherogram library not available for _create_link_object.")
             self.properties['stability_heuristic'] = 'Error'
             return False
        except Exception as e:
            print(f"ERROR: Failed during Link/Braid object creation for '{self.identifier}': {type(e).__name__} - {e}")
            self.properties['stability_heuristic'] = 'Error'
            return False

    def _calculate_properties(self):
        """Calculate standard topological invariants."""
        # --- Special case: Unknot --- #
        if self.identifier == '0_1':
            if not self.link_obj:
                # Attempt to create if not already done (should be done by __init__)
                if not self._create_link_object():
                     print("ERROR: Failed to create link object for unknot 0_1.")
                     self.properties['stability_heuristic'] = 'Error'
                     return

            props = self.properties
            print("Info: Setting properties for Unknot (0_1) directly.")
            props['knot_atlas_name'] = "0_1"
            props['stability_heuristic'] = 'Unknot'
            props['components'] = 1
            props['crossing_number_min'] = 0
            props['signature'] = 0
            props['determinant'] = 1
            props['log_abs_determinant'] = 0.0 # log(1)
            props['volume'] = 0.0
            props['is_alternating'] = True
            props['is_chiral'] = False
            props['is_torus'] = True
            props['is_fibered'] = True
            if SAGE_AVAILABLE:
                try: props['alexander_poly'] = str(T(1)) # Alexander poly is 1
                except Exception: props['alexander_poly'] = '1' # Fallback string
                try: props['jones_poly'] = str(R(1)) # Jones poly is 1
                except Exception: props['jones_poly'] = '1' # Fallback string
            else:
                props['alexander_poly'] = '1'
                props['jones_poly'] = '1'
            props['alex_at_neg1'] = 1
            props['alex_at_neg2'] = 1
            # Jones at roots of unity for J=1
            for name in ROOTS_OF_UNITY:
                props['jones_at_roots'][name] = 1.0
                props['log_abs_jones_at_roots'][name] = 0.0
            return # Finished with unknot

        # --- Proceed for non-Unknot cases --- #
        if not self._create_link_object():
            print(f"Skipping invariant calculation for '{self.identifier}' due to creation failure.")
            # Ensure stability heuristic reflects failure
            if props['stability_heuristic'] != 'Error': props['stability_heuristic'] = 'LinkCreationFailed'
            return # Stability heuristic already set in _create_link_object on error

        l = self.link_obj
        props = self.properties

        try:
            props['components'] = l.num_components() # Use num_components method
        except AttributeError:
             print(f"Warn: Cannot get num_components for {self.identifier} (potentially old spherogram?). Assuming 1.")
             props['components'] = 1 # Fallback assumption
        except Exception as e:
            print(f"Warn: Failed to get number of components for {self.identifier}: {type(e).__name__} - {e}")
            props['components'] = None # Mark as unknown

        # --- Handle Links vs Knots ---
        if props['components'] is None:
            props['stability_heuristic'] = 'ErrorCalculatingComponents'
            return # Cannot proceed
        elif props['components'] > 1:
            props['stability_heuristic'] = 'Link'
            # Try to get a link name (often complex, e.g., L10a140)
            try:
                # Identify may work differently for links, might return list of tuples
                ident_result = l.identify()
                if ident_result and isinstance(ident_result, list):
                     # Simple representation for now
                     props['knot_atlas_name'] = "Link:" + "_".join([spherogram.knot_db.name_from_tuple(comp) for comp in ident_result])
                else:
                     props['knot_atlas_name'] = f"Link ({props['components']} components)"
            except AttributeError:
                print(f"Warn: Link identification ('identify') not available for {self.identifier}.")
                props['knot_atlas_name'] = f"Link ({props['components']} components, unidentified)"
            except Exception as e:
                print(f"Warn: Link identification failed for {self.identifier}: {type(e).__name__} - {e}")
                props['knot_atlas_name'] = f"Link ({props['components']} components, identify error)"
            # Calculate some link properties if needed - focusing on knots for now
            return # Stop further knot-specific calculations
        elif props['components'] == 0:
             props['stability_heuristic'] = 'Empty'
             props['knot_atlas_name'] = "Empty"
             return # Stop

        # --- Single Component: Assume Knot ---
        # Try to convert to a Knot object for potentially more methods
        try:
            # This might fail if the link isn't actually a knot or for other reasons
            self.knot_obj = l.Knot() # Try converting to Knot object
            k = self.knot_obj
            print(f"Info: Successfully created Knot object for {self.identifier}")
        except AttributeError:
             print(f"Warn: Could not convert Link to Knot object for {self.identifier}. Using Link object directly.")
             k = l # Fallback to using the Link object 'l'
        except Exception as e:
             print(f"Warn: Error converting Link to Knot object for {self.identifier}: {type(e).__name__} - {e}. Using Link object.")
             k = l # Fallback

        # --- Knot Identification and Basic Properties ---
        try:
            # Use identify() - should be reliable on Knot object, might work on Link too
            ident_tuple_list = k.identify() # Returns a list of tuples, e.g., [(3,1)] for trefoil
            if not ident_tuple_list: # Empty list
                 print(f"Warn: identify() returned empty list for {self.identifier}. Treating as UnknownPrime.")
                 props['knot_atlas_name'] = "UnknownPrime"
                 props['stability_heuristic'] = 'UnknownPrime'
            elif ident_tuple_list == [(0, 1)]:
                 props['knot_atlas_name'] = "0_1"
                 props['stability_heuristic'] = 'Unknot'
                 # Set invariants for unknot explicitly
                 props['crossing_number_min'] = 0
                 props['signature'] = 0
                 props['determinant'] = 1
                 props['log_abs_determinant'] = -np.inf # log(0) technically, use -inf
                 props['volume'] = 0.0
                 props['is_alternating'] = True # Conventionally
                 props['is_chiral'] = False
                 props['is_torus'] = True # Often considered a (0,1) or (1,1) torus knot
                 props['is_fibered'] = True # Unknot is fibered
                 if SAGE_AVAILABLE:
                     try: props['alexander_poly'] = str(t(1)) # Alexander poly is 1
                     except Exception: pass
                     try: props['jones_poly'] = str(q_var(1)) # Jones poly is 1
                     except Exception: pass
                 props['alex_at_neg1'] = 1
                 props['alex_at_neg2'] = 1
                 # Don't return yet, allow other calculations if needed, but most are set
            else:
                # Assume it's a single prime knot if identify gives one tuple != (0,1)
                # If list has multiple tuples, it's a composite knot identified as components
                if len(ident_tuple_list) == 1:
                    props['knot_atlas_name'] = spherogram.knot_db.name_from_tuple(ident_tuple_list[0])
                    props['stability_heuristic'] = 'Prime'
                else:
                    props['knot_atlas_name'] = " # ".join([spherogram.knot_db.name_from_tuple(comp) for comp in ident_tuple_list])
                    props['stability_heuristic'] = 'CompositeByIdentify' # Identified as composite

                # Try getting crossing number from identified name if possible
                try:
                    if props['stability_heuristic'] == 'Prime':
                         cn_tuple = ident_tuple_list[0]
                         props['crossing_number_min'] = cn_tuple[0]
                except Exception: pass # Ignore errors here

        except AttributeError:
            print(f"Warn: Knot identification ('identify') method not found or failed for {self.identifier}. Status Unknown.")
            props['stability_heuristic'] = 'IdentifyFailed'
        except Exception as e:
            print(f"Warn: Knot identification failed for {self.identifier}: {type(e).__name__} - {e}")
            props['stability_heuristic'] = 'IdentifyError'


        # --- Calculate other invariants (use object 'k' if available, else 'l') ---
        calc_obj = k if k else l # Use Knot object if conversion succeeded, else Link object

        # Check if we already determined it's an unknot
        if props['stability_heuristic'] != 'Unknot':
            # Basic Invariants
            try: props['signature'] = calc_obj.signature()
            except AttributeError: print(f"Warn: 'signature' method not found for {self.identifier}.")
            except ImportError: print(f"Warn: Signature calc failed (likely missing SnapPy/Sage) for {self.identifier}.")
            except Exception as e: print(f"Warn: Signature calc failed for {self.identifier}: {type(e).__name__} - {e}")

            try:
                det = calc_obj.determinant()
                # Spherogram determinant() often returns rational, ensure it's int for knots
                if SAGE_AVAILABLE and isinstance(det, sage.all.Rational):
                    if det.denominator() == 1:
                        det = det.numerator()
                    else: # Should not happen for knot determinant, but handle anyway
                         print(f"Warn: Non-integer determinant {det} for {self.identifier}. Using numerator.")
                         det = det.numerator()

                props['determinant'] = abs(int(det)) # Determinant is usually positive integer
                if props['determinant'] is not None and props['determinant'] > 0:
                     props['log_abs_determinant'] = np.log(props['determinant'])
                elif props['determinant'] == 0 : # Should not happen for knots other than unknot?
                    props['log_abs_determinant'] = -np.inf
            except AttributeError: print(f"Warn: 'determinant' method not found for {self.identifier}.")
            except ImportError: print(f"Warn: Determinant calc failed (likely missing SnapPy/Sage) for {self.identifier}.")
            except Exception as e: print(f"Warn: Determinant calc failed for {self.identifier}: {type(e).__name__} - {e}")

            # Diagrammatic Properties
            try: props['is_alternating'] = calc_obj.is_alternating()
            except AttributeError: print(f"Warn: 'is_alternating' method not found for {self.identifier}.")
            except Exception as e: print(f"Warn: Alternating check failed: {type(e).__name__} - {e}")

            try: props['writhe'] = calc_obj.writhe()
            except AttributeError: print(f"Warn: 'writhe' method not found for {self.identifier}.")
            except Exception as e: print(f"Warn: Writhe calculation failed: {type(e).__name__} - {e}")

            # Chirality (method known to be missing)
            try:
                # This will likely fail based on previous runs, but keep for completeness
                props['is_chiral'] = not calc_obj.is_amphichiral()
            except AttributeError: pass # Expected failure: print(f"Warn: 'is_amphichiral' method not found for {self.identifier}.")
            except Exception as e: print(f"Warn: Chirality check failed: {type(e).__name__} - {e}")

            # Fibering (basic check)
            try: props['is_fibered'] = calc_obj.is_fibered()
            except AttributeError: pass # Method might not exist depending on version/backend
            except NotImplementedError: pass # Some knots might not have this implemented
            except Exception as e: print(f"Warn: Fibered check (basic) failed: {type(e).__name__} - {e}")

            # --- Manifold Properties (using SnapPy via Spherogram) ---
            try:
                # Use exterior() on the original link object 'l' might be safer
                self.manifold = l.exterior()
                vol = self.manifold.volume() # CALL the method!
                # SnapPy volume can return complex with zero imaginary part, take real
                try:
                    # Directly convert Sage RealNumber (or other numeric types) to float
                    props['volume'] = float(vol)
                except TypeError:
                    # Fallback if direct conversion fails (shouldn't happen for RealNumber, but be safe)
                    print(f"Warn: Direct float(vol) failed for {self.identifier}. Type was {type(vol)}. Trying vol.real.")
                    props['volume'] = float(vol.real) if hasattr(vol, 'real') else None
                except Exception as e_vol_conv:
                     print(f"Warn: Exception during volume float conversion for {self.identifier}: {e_vol_conv}")
                     props['volume'] = None

                # Check if Torus based on volume (heuristic) & name
                # Only apply heuristic if volume is near zero AND it wasn't identified as unknot
                is_zero_vol = np.isclose(props['volume'], 0.0, atol=1e-6)
                is_unknot = props['stability_heuristic'] == 'Unknot'
                props['is_torus'] = is_zero_vol and not is_unknot

            except ImportError: pass # SnapPy likely missing
            except AttributeError: pass # exterior() might be missing
            except Exception as e:
                # Don't warn if volume is just zero (common for torus knots)
                current_vol = props.get('volume') # Get volume *after* potential assignment/error
                is_likely_zero = ('Manifold has zero volume' in str(e) or 
                                  (isinstance(current_vol, (float, int)) and np.isclose(current_vol, 0.0)))
                if not is_likely_zero:
                     print(f"Warn: Manifold/Volume calc failed for {self.identifier}: {type(e).__name__} - {e}")
                # Ensure volume is None if calculation truly failed, but keep 0.0 if that was the result
                if props['volume'] is None: # If it never got set correctly in try block
                    props['volume'] = None

        # --- Polynomial Calculations (Requires Sage) ---
        if SAGE_AVAILABLE:
            if props['stability_heuristic'] != 'Unknot': # Skip for unknot if already set
                # Jones Polynomial
                try:
                    jones_poly_q = calc_obj.jones_polynomial(variable=q_var)
                    props['jones_poly'] = str(jones_poly_q) # Store as string
                    # Evaluate at roots of unity
                    for name, root in ROOTS_OF_UNITY.items():
                        try:
                            jones_eval = jones_poly_q(root)
                            # Handle potential complex results: take absolute value (magnitude)
                            # Convert Sage complex number to float before abs if needed
                            if hasattr(jones_eval, 'abs'):
                                jones_abs = jones_eval.abs() # Use Sage's abs() for complex
                            else:
                                jones_abs = abs(complex(jones_eval)) # Convert to Python complex, then abs

                            props['jones_at_roots'][name] = float(jones_abs)
                            # Calculate log only if abs value is significantly non-zero
                            if jones_abs > 1e-9:
                                props['log_abs_jones_at_roots'][name] = float(np.log(jones_abs))
                            else:
                                props['log_abs_jones_at_roots'][name] = -np.inf # log(0) -> -inf

                        except TypeError as te:
                             # Can happen if polynomial eval fails for a root
                             print(f"Warn: Jones eval failed for root '{name}' on {self.identifier}: {te}")
                             props['jones_at_roots'][name] = None
                             props['log_abs_jones_at_roots'][name] = None
                        except Exception as e_eval:
                             print(f"Warn: Jones eval failed unexpectedly for root '{name}' on {self.identifier}: {type(e_eval).__name__} - {e_eval}")
                             props['jones_at_roots'][name] = None
                             props['log_abs_jones_at_roots'][name] = None

                except AttributeError: print(f"Warn: 'jones_polynomial' method not found for {self.identifier}.")
                except ImportError: print(f"Warn: Jones calc failed (likely missing SnapPy/Sage) for {self.identifier}.")
                except Exception as e_jones:
                    print(f"Warn: Jones Polynomial calculation failed for {self.identifier}: {type(e_jones).__name__} - {e_jones}")

                # Alexander Polynomial
                try:
                    # The variable= keyword caused TypeError in last run, remove it.
                    alex_poly_raw = l.alexander_polynomial()

                    # Convert to Sage Polynomial in variable t if possible and needed
                    if SAGE_AVAILABLE and alex_poly_raw is not None and t is not None:
                        try:
                            # Attempt conversion to the polynomial ring T(variable=t)
                            alex_poly_t = T(alex_poly_raw)
                        except Exception as e_conv:
                            print(f"Warn: Could not convert Alexander result to Sage polynomial T(t): {e_conv}")
                            alex_poly_t = None # Failed conversion
                    else:
                        alex_poly_t = alex_poly_raw # Use raw value if not using Sage or it's None

                    if alex_poly_t is not None: # Can be None for multi-component links if unnormalized
                         if alex_poly_t.is_zero():
                             # Alexander poly is 0 for split links, but shouldn't be for knots (unless unknot?)
                             print(f"Info: Alexander polynomial is zero for {self.identifier}.")
                             props['alexander_poly'] = "0"
                             props['alex_at_neg1'] = 0
                             props['alex_at_neg2'] = 0
                         else:
                             # Normalize (make constant term positive, typically)
                             # Spherogram might do this already, but double check
                             # alex_poly_t = alex_poly_t.monic() # Or another normalization? Check convention
                             props['alexander_poly'] = str(alex_poly_t)
                             try:
                                 val_neg1 = alex_poly_t(-1)
                                 # Ensure integer result if possible
                                 props['alex_at_neg1'] = abs(int(val_neg1)) if val_neg1.is_integer() else abs(float(val_neg1))
                             except Exception as e_alex1: print(f"Warn: Alex(-1) eval failed: {e_alex1}")

                             try:
                                 val_neg2 = alex_poly_t(-2)
                                 props['alex_at_neg2'] = abs(int(val_neg2)) if val_neg2.is_integer() else abs(float(val_neg2))
                             except Exception as e_alex2: print(f"Warn: Alex(-2) eval failed: {e_alex2}")
                    else:
                        # If alex_poly_t is None, but raw was not (e.g. Sage conversion failed)
                        if alex_poly_raw is not None:
                            props['alexander_poly'] = str(alex_poly_raw) # Store raw string representation
                        else:
                            print(f"Info: Alexander polynomial returned None for {self.identifier}.")
                            props['alexander_poly'] = None # Explicitly None

                except AttributeError: print(f"Warn: 'alexander_polynomial' method not found for {self.identifier}.")
                except ImportError: print(f"Warn: Alexander calc failed (likely missing SnapPy/Sage) for {self.identifier}.")
                except Exception as e_alex:
                    print(f"Warn: Alexander Polynomial calculation failed for {self.identifier}: {type(e_alex).__name__} - {e_alex}")

        # Final check: if determinant was calculated and alex_at_neg1 was calculated, they should match for a knot.
        det_prop = props.get('determinant')
        alex1_prop = props.get('alex_at_neg1')
        if det_prop is not None and alex1_prop is not None and det_prop != alex1_prop:
             print(f"Warn: Determinant ({det_prop}) and |Alexander(-1)| ({alex1_prop}) mismatch for {self.identifier}.")

        # --- Other Representations and Optional Calculations ---
        if props['stability_heuristic'] != 'Unknot':
            # Braid Word
            try:
                braid_word = calc_obj.braid_word()
                # Store as tuple of integers for potential later use
                props['braid_word_repr'] = tuple(braid_word)
            except AttributeError: print(f"Warn: 'braid_word' method not found for {self.identifier}.")
            except Exception as e:
                print(f"Warn: Braid word calculation failed: {type(e).__name__} - {e}")

            # Morse Number (Optional)
            try:
                props['morse_number'] = calc_obj.morse_number()
            except AttributeError: pass # Morse number method might not exist
            except ImportError: pass # May depend on external solver like GLPK
            except KeyboardInterrupt:
                print(f"Warn: Morse number calculation interrupted for {self.identifier}.")
                props['morse_number'] = None # Set to None if interrupted
            except Exception as e:
                # Avoid verbose warnings for optional calcs unless clearly an error
                if 'GLPK' not in str(e) and 'CBC' not in str(e):
                     print(f"Warn: Morse number calculation failed: {type(e).__name__} - {e}")

            # Knot Floer Homology (Optional - requires external HFKCalculator)
            try:
                hfk_results = calc_obj.knot_floer_homology()
                if isinstance(hfk_results, dict):
                    props['seifert_genus'] = hfk_results.get('seifert_genus')
                    props['is_fibered_hfk'] = hfk_results.get('fibered')
                    # Could store total_rank or rank details if needed
                    # props['hfk_total_rank'] = hfk_results.get('total_rank')
            except ImportError:
                # Usually means HFKCalculator command is not found
                print("Info: Knot Floer Homology calculation skipped (requires external HFKCalculator)." )
                # Only print this info once maybe?
                pass
            except FileNotFoundError:
                 print("Info: Knot Floer Homology calculation skipped (HFKCalculator command not found in PATH)." )
                 pass
            except Exception as e:
                # Avoid printing error for every knot if tool is just missing
                if 'HFKCalculator' not in str(e):
                    print(f"Warn: Knot Floer Homology calculation failed: {type(e).__name__} - {e}")

    def get_property(self, prop_name):
        return self.properties.get(prop_name, None)

    def report(self):
        """Prints a formatted report of the calculated properties."""
        # Use the most definitive name available
        display_name = self.properties.get('knot_atlas_name') or self.properties.get('input_identifier', 'N/A')
        print(f"\n--- Properties for Knot/Link: {display_name} ---")
        if self.properties.get('input_identifier') != display_name and self.properties.get('knot_atlas_name'):
             print(f"  (Input Identifier: {self.properties.get('input_identifier')})")
        if self.properties.get('source_description') != 'Unknown':
            print(f"  (Source: {self.properties.get('source_description')})")

        # Define order or skip certain keys if needed
        skip_keys = {'input_identifier', 'source_description', 'knot_atlas_name'}
        key_order = [
            # Core ID & Topology
            'stability_heuristic', 'components', 'crossing_number_min',
            'signature', 'determinant', 'log_abs_determinant', 'volume',
            # Diagrammatic
            'writhe', 'morse_number', 'is_alternating',
            # Structure
            'is_chiral', 'is_torus', 'is_fibered', 'is_fibered_hfk', 'seifert_genus',
            # Polynomials & Evaluations
            'alexander_poly', 'jones_poly',
            'alex_at_neg1', 'alex_at_neg2', 'jones_at_roots', 'log_abs_jones_at_roots',
            # Representations
            'braid_word_repr',
            # Properties needing more theoretical definition / postulates
            'framing_number', 'topological_charge_Q', 'dynamic_stiffness_k', 'effective_inertia_m_eff'
        ]

        reported_keys = set()

        for key in key_order:
             if key in self.properties and key not in skip_keys:
                 value = self.properties[key]
                 if value == 'UNKNOWN': continue # Skip placeholders explicitly marked
                 if value is None:
                      print(f"  {key:<22}: None")
                 elif isinstance(value, float) or isinstance(value, np.floating):
                      # Handle potential infinities from logs
                      if np.isinf(value):
                          print(f"  {key:<22}: {value}") # Print inf or -inf directly
                      else:
                          print(f"  {key:<22}: {value:.5f}") # More precision for floats
                 elif isinstance(value, dict):
                     print(f"  {key:<22}:")
                     if not value: print("    (empty)")
                     for k,v in sorted(value.items()): # Sort dict items for consistent output
                         if v is None:
                              print(f"    {k}: None")
                         elif isinstance(v, float) or isinstance(v, np.floating):
                             if np.isinf(v): print(f"    {k}: {v}")
                             else: print(f"    {k}: {v:.5f}")
                         else: print(f"    {k}: {v}")
                 elif isinstance(value, str) and len(value) > 60: # Truncate long strings (like polynomials)
                      print(f"  {key:<22}: {value[:57]}...")
                 else:
                      print(f"  {key:<22}: {value}")
                 reported_keys.add(key)

        # Report any remaining properties not in the defined order
        print("  --- Other Properties ---")
        other_keys = set(self.properties.keys()) - reported_keys - skip_keys
        if not other_keys:
            print("    (None)")
        else:
            for key in sorted(list(other_keys)):
                value = self.properties[key]
                if value == 'UNKNOWN': continue
                # Basic print for other keys
                print(f"  {key:<22}: {value}")


        print("-" * 35)

# --- knot_interactions.py ---
# Module for simulating knot interactions (Connected Sum)

def connected_sum(knot_calc1, knot_calc2):
    """
    Performs connected sum conceptually and calculates invariants of result.
    Uses known ADDITIVE/MULTIPLICATIVE properties where possible.
    """
    k1_name = knot_calc1.get_property('knot_atlas_name') or knot_calc1.identifier
    k2_name = knot_calc2.get_property('knot_atlas_name') or knot_calc2.identifier
    print(f"\n--- Simulating Connected Sum: {k1_name} # {k2_name} ---")

    # Basic check: Ensure inputs are single-component knots
    if knot_calc1.get_property('components') != 1 or knot_calc2.get_property('components') != 1:
        print("ERROR: Connected sum is typically defined for knots (1 component).")
        return None
    # Ensure inputs aren't unknots (sum with unknot is the original knot)
    if knot_calc1.get_property('stability_heuristic') == 'Unknot':
        print("Info: Connected sum with Unknot results in the other knot.")
        return knot_calc2.properties # Return properties of the non-unknot summand
    if knot_calc2.get_property('stability_heuristic') == 'Unknot':
        print("Info: Connected sum with Unknot results in the other knot.")
        return knot_calc1.properties # Return properties of the non-unknot summand


    result_name = f"{k1_name} # {k2_name}"
    # Initialize result props based on k1, then update
    result_props = {
        'input_identifier': result_name,
        'source_description': f"ConnectedSum({k1_name}, {k2_name})",
        'knot_atlas_name': result_name,
        'components': 1,
        'stability_heuristic': 'Composite',
        'is_alternating': None, # Alternating property not guaranteed
        'is_torus': False, # Connected sum of non-trivial knots is not torus
        'is_fibered': None, # Fibered property is complex for sums
        'volume': None, # Volume is NOT additive for connected sum (it's 0 for sum)
         'log_abs_determinant': None, # Recalculate below
        # Properties needing interaction model
         'framing_number': 'UNKNOWN',
         'topological_charge_Q': 'UNKNOWN',
         'dynamic_stiffness_k': 'UNKNOWN',
         'effective_inertia_m_eff': 'UNKNOWN',
    }


    # --- Properties that add ---
    sig1 = knot_calc1.get_property('signature')
    sig2 = knot_calc2.get_property('signature')
    if sig1 is not None and sig2 is not None:
        result_props['signature'] = sig1 + sig2
    else: result_props['signature'] = None

    # Crossing Number: Min crossing number is additive for connected sums
    cn1 = knot_calc1.get_property('crossing_number_min')
    cn2 = knot_calc2.get_property('crossing_number_min')
    if cn1 is not None and cn2 is not None:
         result_props['crossing_number_min'] = cn1 + cn2
    else: result_props['crossing_number_min'] = None # If either is unknown


    # --- Properties that multiply ---
    det1 = knot_calc1.get_property('determinant')
    det2 = knot_calc2.get_property('determinant')
    if det1 is not None and det2 is not None:
        result_props['determinant'] = det1 * det2
        if result_props['determinant'] > 0:
             result_props['log_abs_determinant'] = np.log(result_props['determinant'])
        else: # Should only be 0 if one summand was unknot, handled above?
             result_props['log_abs_determinant'] = -np.inf
    else: result_props['determinant'] = None

    alex1_str = knot_calc1.get_property('alexander_poly')
    alex2_str = knot_calc2.get_property('alexander_poly')
    if alex1_str and alex2_str and SAGE_AVAILABLE and t:
        try:
            # Need to parse string back to Sage polynomial using the correct ring T
            p1 = T(alex1_str)
            p2 = T(alex2_str)
            p_prod = p1 * p2
            result_props['alexander_poly'] = str(p_prod)
            # Recalculate evaluations
            try:
                val_neg1 = p_prod(-1)
                result_props['alex_at_neg1'] = abs(int(val_neg1)) if val_neg1.is_integer() else abs(float(val_neg1))
            except Exception: result_props['alex_at_neg1'] = None
            try:
                val_neg2 = p_prod(-2)
                result_props['alex_at_neg2'] = abs(int(val_neg2)) if val_neg2.is_integer() else abs(float(val_neg2))
            except Exception: result_props['alex_at_neg2'] = None
        except Exception as e_alex_sum:
            print(f"Warn: Failed to compute Alexander polynomial for connected sum: {e_alex_sum}")
            result_props['alexander_poly'] = None
            result_props['alex_at_neg1'] = None
            result_props['alex_at_neg2'] = None
    else:
        result_props['alexander_poly'] = None # Cannot calculate if Sage unavailable or polys missing
        result_props['alex_at_neg1'] = None
        result_props['alex_at_neg2'] = None

    # --- Properties that don't combine simply ---
    result_props['jones_poly'] = 'UNKNOWN (Non-multiplicative)'
    result_props['jones_at_roots'] = {} # Cannot easily combine
    result_props['log_abs_jones_at_roots'] = {}

    # Chirality: Sum is chiral if either component is chiral
    # (Sum is amphichiral only if BOTH components are amphichiral)
    chiral1 = knot_calc1.get_property('is_chiral')
    chiral2 = knot_calc2.get_property('is_chiral')
    if chiral1 is None or chiral2 is None:
        result_props['is_chiral'] = None # Unknown if either is unknown
    else:
        result_props['is_chiral'] = chiral1 or chiral2


    print("  Resulting Properties (estimated):")
    # Use the report formatting logic for consistency (simplified here)
    for key, value in sorted(result_props.items()):
         if key in ['framing_number', 'topological_charge_Q', 'dynamic_stiffness_k', 'effective_inertia_m_eff']: continue # Skip placeholders
         if isinstance(value, float): print(f"    {key:<22}: {value:.5f}")
         elif isinstance(value, dict) and not value: print(f"    {key:<22}: {{}}") # Empty dict
         else: print(f"    {key:<22}: {value}")
    print("-" * 35)

    # Return the dictionary of properties. We don't create a full KnotPropertyCalculator
    # as we lack the underlying geometric object (Link/Knot) for the sum.
    return result_props

# --- Braid Multiplication --- #
def multiply_braids(knot_calc1, knot_calc2):
    """
    Performs braid multiplication conceptually and calculates invariants of the result.
    Takes two KnotPropertyCalculator objects, assumes they have braid representations,
    concatenates their braid words, and creates a new KnotPropertyCalculator for the closure.
    """
    k1_name = knot_calc1.get_property('knot_atlas_name') or knot_calc1.identifier
    k2_name = knot_calc2.get_property('knot_atlas_name') or knot_calc2.identifier
    print(f"\n--- Calculating Braid Product: ({k1_name}) * ({k2_name}) ---")

    # Get braid words
    braid1 = knot_calc1.get_property('braid_word_repr')
    braid2 = knot_calc2.get_property('braid_word_repr')

    if braid1 is None or braid2 is None:
        print(f"ERROR: Cannot multiply braids. Missing braid representation for {k1_name} or {k2_name}.")
        return None

    # Concatenate braid words (tuples)
    combined_braid_word = braid1 + braid2

    # Check if the combined braid is trivial (e.g., empty tuple)
    if not combined_braid_word:
        print("Info: Combined braid word is empty. Result is Unknot.")
        # Return a pre-calculated unknot object if available, or create one
        # For simplicity, create one on the fly (might duplicate calculation if 0_1 is already done)
        return KnotPropertyCalculator('0_1', source_description=f"BraidProduct({k1_name}, {k2_name}) -> TrivialBraid")


    print(f"  Combined Braid Word: {combined_braid_word}")
    # Create a new KnotPropertyCalculator for the resulting closed braid
    # The identifier IS the braid word tuple
    try:
        result_calculator = KnotPropertyCalculator(combined_braid_word,
                                                   source_description=f"BraidProduct({k1_name}, {k2_name})")
        return result_calculator
    except Exception as e:
        print(f"ERROR: Failed to create/analyze knot from combined braid {combined_braid_word}: {e}")
        return None

# --- Combination Analysis Function --- #
def analyze_combinations(knot_calculators_dict):
    """
    Iterates through pairs of knots in the dictionary and analyzes their
    connected sum (conceptual) and braid product (generating a new object).
    """
    print("\n" + "="*20 + " Analyzing Knot Combinations " + "="*20)
    names = sorted([name for name in knot_calculators_dict if name != '0_1']) # Exclude unknot
    analyzed_pairs = set()

    for i in range(len(names)):
        for j in range(i, len(names)): # Include self-combinations like 3_1 # 3_1
            name1 = names[i]
            name2 = names[j]

            # Avoid duplicate pairs in reporting if order doesn't matter conceptually
            # (though braid product IS ordered)
            # pair = tuple(sorted((name1, name2)))
            # if pair in analyzed_pairs: continue
            # analyzed_pairs.add(pair)

            print(f"\n--- Analyzing Combination: {name1} and {name2} ---")
            calc1 = knot_calculators_dict[name1]
            calc2 = knot_calculators_dict[name2]

            # 1. Conceptual Connected Sum
            connected_sum_props = connected_sum(calc1, calc2)
            # (Reporting is handled within connected_sum)

            # 2. Braid Multiplication (Order matters: calc1 * calc2)
            braid_product_calc = multiply_braids(calc1, calc2)
            if braid_product_calc:
                braid_product_calc.report()
                # Store results if needed
                # collision_results[(name1, name2)] = braid_product_calc

            # Optional: Calculate braid product in reverse order (calc2 * calc1)
            # if name1 != name2:
            #     print(f"\n--- Analyzing Combination: {name2} and {name1} (Reversed Braid Product) ---")
            #     braid_product_calc_rev = multiply_braids(calc2, calc1)
            #     if braid_product_calc_rev:
            #         braid_product_calc_rev.report()

            print("-"*40) # Separator between pairs

# --- main_analysis.py ---
# Main script to calculate properties and explore

# Define placeholder targets for C1 fit (replace with actual values)
C1_REQUIRED_TARGETS = {
    "electron": 1.0,  # Example value
    "muon": 1.5,      # Example value
    "tau": 2.0,       # Example value
}


if __name__ == "__main__":
    print("\n" + "="*20 + " LSC Knot Property Analysis Engine " + "="*20)

    # --- Phase 1: Calculate Properties of Base Knots ---
    print("\n--- PHASE 1: Calculating Invariants for Base Knots ---")
    # Analyze first ~10 knots (up to 6 crossings) + unknot
    base_knots = [
        "0_1", "3_1", "4_1", "5_1", "5_2",
        "6_1", "6_2", "6_3", # Knots with 6 crossings
        "7_1", "7_2", "7_3", "7_4" # Knots with 7 crossings (first few)
        # Add more standard knot names here as needed
        ]
    knot_calculators = {}
    for name in base_knots:
        print(f"\nCalculating properties for: {name}")
        knot_calculators[name] = KnotPropertyCalculator(name, source_description="Known Atlas Knot")
        knot_calculators[name].report()

    # --- Phase 2: Analyze Knot Combinations --- 
    analyze_combinations(knot_calculators)

    # --- Phase 3: Data Analysis for Physics Mapping (Example) ---
    print("\n--- PHASE 3: Data Ready for Physics Correlation ---")
    print("Invariant data calculated and stored in 'knot_calculators' dictionary.")
    # Example: Extract data needed for c1 fit (as done in previous script)
    leptons_assigned = {"electron": "3_1", "muon": "5_1", "tau": "5_2"}
    print("\nData for Lepton Candidates (e=3_1, mu=5_1, tau=5_2):")
    fit_data = {}
    for lepton, knot_name in leptons_assigned.items():
        if knot_name in knot_calculators:
             calc = knot_calculators[knot_name]
             target_c1 = C1_REQUIRED_TARGETS.get(lepton) # Get target value

             # Extract specific invariants needed for c1 fit
             # Use .get() with default None for safety if a property wasn't calculated
             log_det_val = calc.get_property('log_abs_determinant')
             sig_val = calc.get_property('signature')
             jones_w5_logabs_dict = calc.get_property('log_abs_jones_at_roots')
             jones_w5_val = jones_w5_logabs_dict.get('w5') if isinstance(jones_w5_logabs_dict, dict) else None
             alex_neg1_val = calc.get_property('alex_at_neg1')

             data_pt = {
                 'req_c1': target_c1,
                 'log_det': log_det_val,
                 'sig': sig_val,
                 'jones_w5_logabs': jones_w5_val,
                 'alex_neg1': alex_neg1_val
             }
             fit_data[knot_name] = data_pt

             # Format output nicely, handling None values
             log_det_str = f"{log_det_val:.3f}" if log_det_val is not None else 'N/A'
             sig_str = f"{sig_val}" if sig_val is not None else 'N/A'
             jones_str = f"{-jones_w5_val:.3f}" if jones_w5_val is not None and np.isfinite(jones_w5_val) else ('0.000' if jones_w5_val == -np.inf else 'N/A') # Treat log(0) as 0? Check requirement
             alex_str = f"{alex_neg1_val}" if alex_neg1_val is not None else 'N/A'
             req_c1_str = f"{target_c1:.3f}" if target_c1 is not None else 'N/A'

             print(f"  {knot_name}: Req_c1={req_c1_str}, log|Det|={log_det_str}, Sig={sig_str}, "
                   f"-log|J(ω5)|={jones_str}, |Δ(-1)|={alex_str}")
        else:
             print(f"  {knot_name}: Data not calculated.")

    # Here, one would import the fitting functions and run the correlation
    # analysis on the 'fit_data' dictionary as performed in the previous script.
    # Example: result = perform_linear_fit(fit_data)
    # print(f"\nFit Result: {result}")
    print("\n>>> Analysis can now proceed using the calculated invariants <<<")

    print("\n--- Analysis Engine Finished ---")