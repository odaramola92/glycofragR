"""
High-level unified API for glycopeptide analysis and structure prediction.

Single-line initialization that handles everything internally:
- Glycopeptide creation (Peptide + Glycan)
- Fragment generation
- Theoretical table generation (raw + clean)
- Summary sheet creation
- Optional visualization

Example:
    >>> from glycofrag import GlycoPeptideAnalysis
    >>> 
    >>> # Everything in one call
    >>> analysis = GlycoPeptideAnalysis(
    ...     peptide_sequence='EEQYNSTYR',
    ...     glycan_code='4501',
    ...     glycosylation_site=5,
    ...     glycan_type='N',
    ...     peptide_fragment_types=['by', 'cz'],
    ...     glycan_fragment_types=['BY'],
    ...     visualize_structure='best'  # or 'all', [1,2,3], None
    ... )
    >>> 
    >>> # Access generated data
    >>> print(analysis.theoretical_df)  # raw theoretical fragments
    >>> print(analysis.clean_theoretical_df)  # deduplicated
    >>> print(analysis.summary_df)  # mass information
    >>> 
    >>> # For MS/MS matching
    >>> matched = analysis.match_fragments(spectrum, ppm_tolerance=20)
    >>> prediction = analysis.predict_structures(matched)
    >>> analysis.visualize_best_structure('structure.png')
"""

import os
from typing import Optional, Dict, List, Tuple, Union
import numpy as np
import pandas as pd
from glycofrag.glycopeptide import Glycopeptide
from glycofrag.io.tables import (
    build_clean_theoretical,
    score_structures_by_uniqueness,
    create_structure_prediction_sheet,
    create_theoretical_summary_sheet,
    get_best_structure,
    generate_fragment_table,
    match_fragments_from_table,
)
from glycofrag.io.visualizer import GlycanVisualizer
from glycofrag.core.mass_calculator import GlycanMassCalculator


class GlycoPeptideAnalysis:
    """
    All-in-one interface for complete glycopeptide MS/MS analysis pipeline.
    
    Initialization handles ALL setup:
    - Creates Glycopeptide (which includes Peptide + Glycan)
    - Generates fragment ions
    - Builds theoretical fragment tables (raw + deduplicated)
    - Creates summary sheet with masses/modifications
    - Optionally visualizes structures
    
    Then provides simple methods for:
    - Fragment matching (experimental)
    - Structure prediction/ranking
    - Visualization and export
    
    Example:
        >>> analysis = GlycoPeptideAnalysis(
        ...     peptide_sequence='EEQYNSTYR',
        ...     glycan_code='4501',
        ...     glycosylation_site=5,
        ...     glycan_type='N',
        ...     peptide_fragment_types=['by', 'cz'],
        ...     glycan_fragment_types=['BY'],
        ...     visualize_structure='best'
        ... )
        >>>
        >>> # Access all generated data
        >>> print(analysis.theoretical_df)  # 45 fragments
        >>> print(analysis.clean_theoretical_df)  # 43 deduplicated
        >>> print(analysis.mass)  # glycopeptide neutral mass
        >>>
        >>> # Match experimental peaks
        >>> matched_df = analysis.match_fragments(spectrum, ppm_tolerance=20)
        >>>
        >>> # Rank structures
        >>> prediction = analysis.predict_structures(matched_df)
        >>> analysis.visualize_best_structure('best_structure.png')
    """
    
    def __init__(
        self,
        peptide_sequence: str,
        glycan_code: str,
        glycosylation_site: int,
        glycan_type: str = 'N',
        use_cam: bool = True,
        max_structures: int = 100,
        mod_string: Optional[str] = None,
        peptide_fragment_types: Optional[List[str]] = None,
        glycan_fragment_types: Optional[List[str]] = None,
        isomer_sensitive: bool = False,
        charges: Optional[List[int]] = None,
        custom_mono_masses: Optional[Dict[str, float]] = None,
        visualize_structure: Optional[Union[str, int, List[int]]] = None,
        preferred_core: Optional[int] = None,
    ):
        """
        Initialize and run complete glycopeptide analysis pipeline.
        
        Args:
            peptide_sequence: Amino acid sequence (e.g., 'EEQYNSTYR')
            glycan_code: Glycan composition code (e.g., '4501')
            glycosylation_site: 1-indexed N position in N-X-S/T motif
            glycan_type: 'N' or 'O' for N/O-glycan (default 'N')
            use_cam: Apply carbamidomethylation on cysteines (default True)
            max_structures: Max structures to predict (default 100)
            mod_string: Modification string, e.g., 'M:Ox; K:Ac' (optional)
            peptide_fragment_types: Peptide fragment types to generate
                                   e.g., ['by', 'cz'] or None to skip
            glycan_fragment_types: Glycan fragment types to generate
                                  e.g., ['BY'] or None to skip
            isomer_sensitive: If True, treats mirror images as distinct structures
                            (exhaustive enumeration). If False, deduplicates
                            topologically identical mirror images (default False)
            charges: Charge states for m/z calculation (default [1, 2, 3])
                    Example: [1, 2, 3, 4] for extended charge states
            custom_mono_masses: Dict of mass (Da) to ADD to specific monosaccharides.
                              Keys: 'HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc'
                              Values: Mass to add (can be positive or negative)
                              Examples:
                                - {'NeuAc': 42.0106} for O-acetylated NeuAc
                                - {'NeuAc': 50.0657} for custom modification
                                - {'Hex': -18.0106} for mass loss (rare)
            visualize_structure: Which structures to visualize
                               - 'best': Top-ranked structure
                               - 'all': All predicted structures
                               - int (1, 2, 3...): Specific structure
                               - list [1, 2, 3]: Multiple structures
                               - None: No visualization
            preferred_core: For O-glycans with ambiguous cores (e.g., Core 1 vs 8).
                          Specify which core to use: 1, 3, 5, 6, 7, or 8.
                          If None, defaults to lowest core number.
                          For N-glycans, this parameter is ignored.
        
        Generates and stores internally:
            - theoretical_df: All theoretical fragments (raw)
            - clean_theoretical_df: Deduplicated by mass+composition
            - summary_df: Masses, modifications, QC info
            - structures: Glycan structure objects
        """
        # Store parameters
        self.peptide_sequence = peptide_sequence
        self.glycan_code = glycan_code
        self.glycosylation_site = glycosylation_site
        self.glycan_type = glycan_type
        self.mod_string = mod_string
        self.charges = charges if charges is not None else [1, 2, 3]
        self.glycan = None
        
        # === STEP 1: Create Glycopeptide (internal) ===
        self.gp = Glycopeptide(
            peptide_sequence,
            glycan_code,
            glycosylation_site,
            glycan_type,
            use_cam=use_cam,
            max_structures=max_structures,
            mod_string=mod_string or "",
            isomer_sensitive=isomer_sensitive,
            custom_mono_masses=custom_mono_masses,
            preferred_core=preferred_core,
        )
        self.glycan = self.gp.glycan
        
        # === STEP 2: Generate fragments (if requested) ===
        self.fragments = None
        if peptide_fragment_types or glycan_fragment_types:
            self.fragments = self.gp.generate_fragments(
                peptide_fragment_types=peptide_fragment_types or [],
                glycan_fragment_types=glycan_fragment_types or [],
                charges=self.charges,
            )
        
        # === STEP 3: Build theoretical tables (internal) ===
        if self.fragments:
            self.theoretical_df = generate_fragment_table(
                fragments=self.fragments,
                glycan_code=glycan_code,
                modification_type='glycopeptide',
                peptide=peptide_sequence,
                charge_states=self.charges,
            )
        else:
            self.theoretical_df = pd.DataFrame()
        
        # Deduplicate by mass + composition
        self.clean_theoretical_df = build_clean_theoretical(
            self.theoretical_df, 
            mass_round=6
        ) if not self.theoretical_df.empty else pd.DataFrame()
        
        # === STEP 4: Create summary sheet (internal) ===
        self.summary_df = create_theoretical_summary_sheet(
            peptide_sequence=peptide_sequence,
            glycan_code=glycan_code,
            peptide_mass=self.gp.peptide.residue_mass,
            glycan_mass=self.gp.mass_calculator.calculate_glycan_mass(glycan_code),
            glycopeptide_mass=self.gp.mass,
            mod_string=mod_string,
            glycosylation_site=glycosylation_site,
            glycan_type=glycan_type,
        )
        
        # Cache for analysis results
        self.prediction_df: Optional[pd.DataFrame] = None
        self.scores: Optional[Dict] = None
        
        # === STEP 5: Optional visualization ===
        if visualize_structure:
            self._visualize_init(visualize_structure)
    
    def _visualize_init(self, visualize_structure: Union[str, int, List[int]]):
        """Visualize structures during initialization."""
        if visualize_structure == 'all':
            # Visualize all structures
            for i in range(1, min(len(self.structures) + 1, 10)):  # Max 9 to prevent too many
                try:
                    GlycanVisualizer.visualize_with_matplotlib(
                        self.structures[i - 1],
                        title=f"Structure {i}",
                        show=True
                    )
                except Exception as e:
                    print(f"[WARNING] Could not visualize structure {i}: {e}")
        
        elif visualize_structure == 'best' or visualize_structure == 1:
            # Visualize best (first) structure
            try:
                GlycanVisualizer.visualize_with_matplotlib(
                    self.structures[0],
                    title=f"Structure 1 (Best)",
                    show=True
                )
            except Exception as e:
                print(f"[WARNING] Could not visualize best structure: {e}")
        
        elif isinstance(visualize_structure, int):
            # Visualize specific structure
            idx = visualize_structure - 1
            if 0 <= idx < len(self.structures):
                try:
                    GlycanVisualizer.visualize_with_matplotlib(
                        self.structures[idx],
                        title=f"Structure {visualize_structure}",
                        show=True
                    )
                except Exception as e:
                    print(f"[WARNING] Could not visualize structure: {e}")
        
        elif isinstance(visualize_structure, list):
            # Visualize multiple structures
            for struct_id in visualize_structure:
                idx = struct_id - 1
                if 0 <= idx < len(self.structures):
                    try:
                        GlycanVisualizer.visualize_with_matplotlib(
                            self.structures[idx],
                            title=f"Structure {struct_id}",
                            show=True
                        )
                    except Exception as e:
                        print(f"[WARNING] Could not visualize structure {struct_id}: {e}")
    
    def match_fragments(
        self,
        experimental_spectrum: Tuple[np.ndarray, np.ndarray],
        exp_rt_map: Optional[Dict[float, float]] = None,
        exp_area_map: Optional[Dict[float, float]] = None,
        ppm_tolerance: float = 20.0,
        charge_states: List[int] = None,
        intensity_threshold: float = 0.01,
    ) -> pd.DataFrame:
        """
        Match experimental peaks with theoretical fragments.
        
        Args:
            experimental_spectrum: (mz_array, intensity_array) tuple
            exp_rt_map: Dict mapping m/z to retention time (optional)
            exp_area_map: Dict mapping m/z to peak area (optional)
            ppm_tolerance: Fragment matching tolerance (default 20 ppm)
            charge_states: Charge states to try (default [1, 2, 3])
            intensity_threshold: Min intensity as fraction of max (default 0.01)
        
        Returns:
            DataFrame with matched fragments:
                Columns: Fragment, Composition, Structure, Charge,
                        Experimental m/z, Theoretical m/z, PPM Error,
                        Intensity (%), RT, Area
        
        Example:
            >>> matched = analysis.match_fragments(spectrum, ppm_tolerance=10)
            >>> print(len(matched), 'fragments matched')
        """
        if charge_states is None:
            charge_states = [1, 2, 3]
        
        matched_df = match_fragments_from_table(
            self.clean_theoretical_df,
            experimental_spectrum,
            exp_rt_map,
            exp_area_map,
            ppm_tolerance=ppm_tolerance,
            charge_states=charge_states,
            intensity_threshold=intensity_threshold,
        )
        
        return matched_df
    
    def predict_structures(
        self,
        matched_df: pd.DataFrame,
        mass_round: int = 6,
    ) -> pd.DataFrame:
        """
        Rank structures based on fragment matching uniqueness.
        
        For each structure, calculates:
        - Match %: matched unique / total unique in structure
        - Match Share %: contribution to total matched fragments
        - Confidence Level: relative strength (0-10)
        
        Args:
            matched_df: DataFrame from match_fragments()
            mass_round: Decimal places for deduplication (default 6)
        
        Returns:
            DataFrame with structure ranking:
                Columns: Structure, Unique Theoretical, Matched Unique,
                        Match %, Match Share %, Rank, Confidence Level (0-10)
        
        Example:
            >>> prediction = analysis.predict_structures(matched_df)
            >>> print(prediction[['Structure', 'Match %', 'Confidence Level']])
               Structure  Match %  Confidence Level (0-10)
            0          1     88.9                       8.9
            1          2     53.2                       5.3
        """
        if matched_df.empty:
            raise ValueError("No matched fragments available for prediction")
        
        # Score each structure
        self.scores = score_structures_by_uniqueness(
            self.clean_theoretical_df,
            matched_df
        )
        
        # Create prediction report
        self.prediction_df = create_structure_prediction_sheet(
            self.clean_theoretical_df,
            matched_df,
            self.scores
        )
        
        return self.prediction_df
    
    def visualize_predicted_structures(
        self,
        matched_df: pd.DataFrame,
        which: Union[str, List[int]] = "best",
        output_dir: Optional[str] = None,
        output_excel: Optional[str] = None,
        figsize: tuple = (5, 4),
    ) -> Dict[int, Optional[str]]:
        """
        Visualize predicted structures with flexible selection options.
        Internally handles structure prediction and ranking based on fragment matching.
        
        This is the primary method for structure visualization - no need to call 
        predict_structures() separately.
        
        Args:
            matched_df: DataFrame from match_fragments() with experimental matches
            which: Which structures to visualize:
                - "best" (default): Visualize all structures tied for top-ranked.
                  When multiple structures have identical best score, visualizes all.
                - "all": Visualize all predicted structures
                - [1, 2, 5]: Visualize specific structure IDs by number
            output_dir: Directory to save PNG files. If None, displays interactively.
                       Filenames: structure_1.png, structure_2.png, etc.
            output_excel: Path to Excel file to append Structure_Prediction sheet
            figsize: Tuple of (width, height) in inches (default: (5, 4))
        
        Returns:
            Dict mapping structure ID to output path:
            {1: 'structure_1.png', 2: 'structure_2.png', ...}
            or {1: None} if displayed interactively
        
        Examples:
            >>> # Visualize best structure only
            >>> analysis.visualize_predicted_structures(
            ...     matched_df,
            ...     which="best",
            ...     output_dir="output/",
            ...     figsize=(8, 6)
            ... )
            {1: 'output/structure_1.png'}
            
            >>> # Visualize all structures
            >>> analysis.visualize_predicted_structures(
            ...     matched_df,
            ...     which="all",
            ...     output_dir="output/"
            ... )
            {1: 'output/structure_1.png', 2: 'output/structure_2.png', ...}
            
            >>> # Visualize specific structures
            >>> analysis.visualize_predicted_structures(
            ...     matched_df,
            ...     which=[1, 3],
            ...     output_dir="output/"
            ... )
            {1: 'output/structure_1.png', 3: 'output/structure_3.png'}
        """
        if matched_df.empty:
            raise ValueError("No matched fragments provided")
        
        # Auto-predict structures if not already done
        if self.prediction_df is None or self.prediction_df.empty:
            self.predict_structures(matched_df)
        
        # Check if prediction was successful
        if self.prediction_df is None or self.prediction_df.empty:
            print(f"\n[ERROR] Structure prediction failed!")
            print(f"\nDiagnostics:")
            print(f"  Glycan code: {self.glycan_code}")
            print(f"  Glycan type: {self.glycan_type}")
            print(f"  Structures generated: {len(self.structures)}")
            print(f"  Theoretical fragments: {len(self.theoretical_df)}")
            print(f"  Matched fragments: {len(matched_df)}")
            
            if len(self.structures) == 0:
                print(f"\n[ROOT CAUSE] No glycan structures were generated for code '{self.glycan_code}'")
                print(f"  This glycan composition may not be valid for {self.glycan_type}-glycans.")
                print(f"  Glycan code format: HexNAc-Hex-NeuAc-Fuc (4 digits)")
                print(f"  Example valid codes: '4501' (HexNAc=4,Hex=5,NeuAc=0,Fuc=1)")
                print(f"                       '2210' (HexNAc=2,Hex=2,NeuAc=1,Fuc=0)")
                raise ValueError(f"Invalid glycan code '{self.glycan_code}' - no structures generated")
            elif len(matched_df) == 0:
                print(f"\n[ROOT CAUSE] No fragments matched experimental data")
                raise ValueError("No matched fragments available for structure prediction")
            else:
                print(f"\n[ROOT CAUSE] Structure scoring returned empty results")
                print(f"  Check that matched fragments have valid 'Structure' column values")
                raise ValueError("Structure prediction scoring failed")
        
        # Determine which structures to visualize
        if 'Structure' not in self.prediction_df.columns:
            print(f"[ERROR] 'Structure' column not found in prediction DataFrame")
            print(f"Available columns: {list(self.prediction_df.columns)}")
            print(f"Clean theoretical DataFrame columns: {list(self.clean_theoretical_df.columns)}")
            raise KeyError("'Structure' column missing from prediction results. Check theoretical fragment generation.")
        
        # Sort structures by Rank (NOT by structure ID number!)
        # prediction_df already has Rank column with proper ordering
        ranked_df = self.prediction_df.sort_values('Rank')
        all_structure_ids = [int(s) for s in ranked_df['Structure'].astype(str).str.strip()]
        
        if which == "best":
            # Select all structures with the best (highest) Weighted Score
            # This handles ties - when multiple structures have identical scores
            best_score = self.prediction_df['Weighted Score'].max()
            best_structures = self.prediction_df[self.prediction_df['Weighted Score'] == best_score]
            struct_ids = [int(s) for s in best_structures['Structure'].astype(str).str.strip()]
            
            if len(struct_ids) > 1:
                print(f"[INFO] {len(struct_ids)} structures tied for best (score={best_score}). Visualizing all.")
        elif which == "all":
            struct_ids = all_structure_ids
        elif isinstance(which, (list, tuple)):
            struct_ids = [int(s) for s in which]
            # Validate that requested structures exist
            invalid = set(struct_ids) - set(all_structure_ids)
            if invalid:
                raise ValueError(f"Requested structure ID(s) {invalid} not found. Available: {all_structure_ids}")
        else:
            raise ValueError(f"Invalid 'which' parameter: {which}. Use 'best', 'all', or list of IDs")
        
        # Append prediction sheet to Excel once
        if output_excel:
            try:
                with pd.ExcelWriter(
                    output_excel,
                    engine='openpyxl',
                    mode='a',
                    if_sheet_exists='replace'
                ) as writer:
                    self.prediction_df.to_excel(
                        writer,
                        sheet_name='Structure_Prediction',
                        index=False
                    )
                print(f"[OK] Structure_Prediction sheet appended: {output_excel}")
            except Exception as e:
                print(f"[WARNING] Could not update Excel: {e}")
        
        # Visualize structures using Glycopeptide's built-in tree visualization
        # This uses the proper hierarchical tree layout from the structure objects
        output_paths = {}
        
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            # Use Glycopeptide's visualization method which handles tree-based rendering
            saved_files = self.gp.visualize_structures(
                structures_to_draw=struct_ids,
                figsize=figsize,
                save_images=True,
                show_plots=False,
                output_dir=output_dir
            )
            
            # Map structure IDs to output paths
            for i, struct_id in enumerate(struct_ids, 1):
                output_path = f"{output_dir}/structure_{struct_id}.png"
                output_paths[struct_id] = output_path
                print(f"[OK] Structure {struct_id} saved: {output_path}")
        else:
            # Display interactively
            self.gp.visualize_structures(
                structures_to_draw=struct_ids,
                figsize=figsize,
                save_images=False,
                show_plots=True,
                output_dir=None
            )
            
            for struct_id in struct_ids:
                output_paths[struct_id] = None
        
        return output_paths
    
    def visualize_best_structure(
        self,
        output_image: Optional[str] = None,
        output_excel: Optional[str] = None,
        figsize: tuple = (5, 4),
    ) -> Optional[str]:
        """
        DEPRECATED: Use visualize_predicted_structures() instead.
        
        This method is kept for backwards compatibility.
        Visualizes the best structure from a previous predict_structures() call.
        
        Args:
            output_image: Path to save PNG
            output_excel: Path to Excel file
            figsize: Figure size in inches
        
        Returns:
            Path to saved image
        """
        import warnings
        warnings.warn(
            "visualize_best_structure() is deprecated. "
            "Use visualize_predicted_structures(matched_df, which='best') instead.",
            DeprecationWarning,
            stacklevel=2
        )
        
        if self.prediction_df is None or self.prediction_df.empty:
            raise ValueError("Run predict_structures() first")
        
        best_struct_id = get_best_structure(self.prediction_df)
        if best_struct_id is None:
            raise ValueError("No best structure found")
        
        struct_idx = int(best_struct_id)
        structure = self.structures[struct_idx - 1]
        
        # Update Excel if provided
        if output_excel:
            try:
                with pd.ExcelWriter(
                    output_excel,
                    engine='openpyxl',
                    mode='a',
                    if_sheet_exists='replace'
                ) as writer:
                    self.prediction_df.to_excel(
                        writer,
                        sheet_name='Structure_Prediction',
                        index=False
                    )
                print(f"[OK] Structure_Prediction appended: {output_excel}")
            except Exception as e:
                print(f"[WARNING] Could not update Excel: {e}")
        
        # Visualize
        show = output_image is None
        GlycanVisualizer.visualize_with_matplotlib(
            structure,
            title=f"Most Probable Structure (Structure {struct_idx})",
            figsize=figsize,
            show=show,
            output_path=output_image
        )
        
        if output_image:
            print(f"[OK] Structure saved: {output_image}")
            return output_image
        return None
    
    def export_to_excel(
        self,
        output_file: str,
        matched_df: Optional[pd.DataFrame] = None,
    ) -> str:
        """
        Export all generated data to Excel.
        
        Creates sheets:
        - Summary: Masses, modifications, QC info
        - Raw_Theoretical: All generated fragments (70+)
        - Clean_Theoretical: Deduplicated fragments (68)
        - Matched: Matched experimental fragments (optional)
        - Structure_Prediction: Ranking results (if predict_structures called)
        
        Args:
            output_file: Output Excel file path
            matched_df: Optional matched fragments DataFrame
        
        Returns:
            Confirmation message
        
        Example:
            >>> analysis.export_to_excel('analysis_results.xlsx', matched_df)
            [OK] Results exported to: analysis_results.xlsx
               Sheet 0: Summary (Masses & Modifications)
               Sheet 1: Raw_Theoretical (70 fragments)
               Sheet 2: Clean_Theoretical (68 fragments)
               Sheet 3: Matched (14 fragments)
        """
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Sheet 1: Summary
            self.summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Sheet 2: Raw Theoretical
            self.theoretical_df.to_excel(
                writer,
                sheet_name='Raw_Theoretical',
                index=False
            )
            
            # Sheet 3: Clean Theoretical
            self.clean_theoretical_df.to_excel(
                writer,
                sheet_name='Clean_Theoretical',
                index=False
            )
            
            # Sheet 4: Matched (if provided)
            if matched_df is not None and not matched_df.empty:
                matched_df.to_excel(
                    writer,
                    sheet_name='Matched',
                    index=False
                )
            
            # Sheet 5: Structure Prediction (if available)
            if self.prediction_df is not None and not self.prediction_df.empty:
                self.prediction_df.to_excel(
                    writer,
                    sheet_name='Structure_Prediction',
                    index=False
                )
        
        sheet_count = 3
        if matched_df is not None and not matched_df.empty:
            sheet_count += 1
        if self.prediction_df is not None:
            sheet_count += 1
        
        print(f"\n[OK] Results exported to: {output_file}")
        print(f"   Summary, Raw_Theoretical ({len(self.theoretical_df)} fragments),")
        print(f"   Clean_Theoretical ({len(self.clean_theoretical_df)} fragments)")
        if matched_df is not None and not matched_df.empty:
            print(f"   Matched ({len(matched_df)} fragments)")
        if self.prediction_df is not None:
            print(f"   Structure_Prediction (ranking results)")
        
        return f"Exported to {output_file}"
    
    # ========== R-COMPATIBLE CONVERSION METHODS ==========
    # These methods convert pandas DataFrames to R-friendly formats
    # (list of dictionaries) for seamless interchange with R via reticulate
    
    def get_theoretical_df_as_list(self) -> List[Dict]:
        """
        Convert theoretical_df to R-compatible list of dictionaries.
        
        Each row becomes a dictionary, making it easy for R to work with.
        Handles NaN/None values appropriately.
        
        Returns:
            List of dictionaries, one per row in theoretical_df.
            Empty list if no theoretical fragments.
        """
        if self.theoretical_df.empty:
            return []
        return self.theoretical_df.fillna("").to_dict('records')
    
    def get_clean_df_as_list(self) -> List[Dict]:
        """
        Convert clean_theoretical_df to R-compatible list of dictionaries.
        
        Each row becomes a dictionary, making it easy for R to work with.
        Handles NaN/None values appropriately.
        
        Returns:
            List of dictionaries, one per row in clean_theoretical_df.
            Empty list if no clean fragments.
        """
        if self.clean_theoretical_df.empty:
            return []
        return self.clean_theoretical_df.fillna("").to_dict('records')
    
    def get_summary_df_as_list(self) -> List[Dict]:
        """
        Convert summary_df to R-compatible list of dictionaries.
        
        Each row becomes a dictionary, making it easy for R to work with.
        Handles NaN/None values appropriately.
        
        Returns:
            List of dictionaries, one per row in summary_df.
            Empty list if no summary data.
        """
        if self.summary_df.empty:
            return []
        return self.summary_df.fillna("").to_dict('records')
    
    def get_structures_as_list(self) -> List[Dict]:
        """
        Convert predicted structures to R-compatible format.
        
        Returns:
            List of dictionaries representing each predicted structure.
            Each dict contains structure information and string representation.
        """
        structures_list = []
        for idx, struct in enumerate(self.structures, 1):
            structures_list.append({
                'structure_id': idx,
                'structure_str': str(struct),
                'structure_repr': repr(struct)
            })
        return structures_list
    
    # Properties for clean access
    @property
    def mass(self) -> float:
        """Glycopeptide neutral mass (Da)."""
        return self.gp.mass
    
    @property
    def peptide_mass(self) -> float:
        """Peptide residue mass (Da)."""
        return self.gp.peptide.residue_mass
    
    @property
    def glycan_mass(self) -> float:
        """Glycan neutral mass (Da)."""
        return self.gp.mass_calculator.calculate_glycan_mass(self.glycan_code)
    
    @property
    def structures(self) -> List:
        """Predicted glycan structures (list of structure objects)."""
        return self.gp.glycan_structures
    
    @property
    def num_structures(self) -> int:
        """Number of predicted structures."""
        return len(self.structures)
