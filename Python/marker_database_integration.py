#!/usr/bin/env python3

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union
import logging

class ScoreCalculator:
    """
    Calculate Evidence Consistency Score (ECs) and Specificity Score (SPs) 
    for marker gene databases, following the methodology from CellMarkerAccordion R package.
    """
    
    def __init__(self, log_level: str = 'INFO'):
        """Initialize the score calculator."""
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
        
    def calculate_ecs_global(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate global Evidence Consistency Score (ECs).
        
        ECs measures how frequently a marker gene is reported for a specific 
        cell type across different resources.
        
        Args:
            df: DataFrame with columns: species, CL_celltype, CL_ID, marker, marker_type, resource
            
        Returns:
            DataFrame with ECs_global scores
        """
        # Group by marker-celltype combination and count unique resources
        required_cols = ['species', 'CL_celltype', 'CL_ID', 'marker', 'marker_type', 'resource']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for ECs calculation: {missing_cols}")
        
        # Get unique entries to avoid double counting
        unique_entries = df[required_cols].drop_duplicates()
        
        # Count occurrences (equivalent to R's ddply nrow)
        ecs_counts = unique_entries.groupby(
            ['species', 'CL_celltype', 'CL_ID', 'marker', 'marker_type']
        ).size().reset_index(name='ECs_global')
        
        self.logger.debug(f"Calculated global ECs for {len(ecs_counts)} marker-celltype combinations")
        return ecs_counts
    
    def calculate_ecs_tissue_specific(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate tissue-specific Evidence Consistency Score (ECs).
        
        Args:
            df: DataFrame with columns: species, Uberon_tissue, Uberon_ID, CL_celltype, CL_ID, marker, marker_type, resource
            
        Returns:
            DataFrame with ECs_tissue_specific scores
        """
        required_cols = ['species', 'Uberon_tissue', 'Uberon_ID', 'CL_celltype', 'CL_ID', 'marker', 'marker_type', 'resource']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for tissue-specific ECs calculation: {missing_cols}")
        
        # Get unique entries to avoid double counting
        unique_entries = df[required_cols].drop_duplicates()
        
        # Count occurrences per tissue
        ecs_counts = unique_entries.groupby(
            ['species', 'Uberon_tissue', 'Uberon_ID', 'CL_celltype', 'CL_ID', 'marker', 'marker_type']
        ).size().reset_index(name='ECs_tissue_specific')
        
        self.logger.debug(f"Calculated tissue-specific ECs for {len(ecs_counts)} marker-celltype-tissue combinations")
        return ecs_counts
    
    def calculate_sps_global(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate global Specificity Score (SPs).
        
        SPs measures how specific a marker is to particular cell types.
        Higher specificity (lower number of cell types) = higher score.
        
        Args:
            df: DataFrame with columns: species, marker, marker_type, CL_celltype
            
        Returns:
            DataFrame with SPs_global scores
        """
        required_cols = ['species', 'marker', 'marker_type', 'CL_celltype']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for global SPs calculation: {missing_cols}")
        
        # Get unique marker-celltype combinations
        unique_entries = df[required_cols].drop_duplicates()
        
        # Count unique cell types per marker
        celltype_counts = unique_entries.groupby(
            ['species', 'marker', 'marker_type']
        )['CL_celltype'].nunique().reset_index(name='celltype_count')
        
        # Calculate specificity: SPs = 1 / number_of_cell_types
        celltype_counts['SPs_global'] = (1 / celltype_counts['celltype_count']).round(2)
        
        # Drop intermediate column
        sps_scores = celltype_counts.drop('celltype_count', axis=1)
        
        self.logger.debug(f"Calculated global SPs for {len(sps_scores)} marker combinations")
        return sps_scores
    
    def calculate_sps_tissue_specific(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate tissue-specific Specificity Score (SPs).
        
        Args:
            df: DataFrame with columns: species, Uberon_tissue, Uberon_ID, marker, marker_type, CL_celltype
            
        Returns:
            DataFrame with SPs_tissue_specific scores
        """
        required_cols = ['species', 'Uberon_tissue', 'Uberon_ID', 'marker', 'marker_type', 'CL_celltype']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for tissue-specific SPs calculation: {missing_cols}")
        
        # Get unique marker-celltype-tissue combinations
        unique_entries = df[required_cols].drop_duplicates()
        
        # Count unique cell types per marker within each tissue
        celltype_counts = unique_entries.groupby(
            ['species', 'Uberon_tissue', 'Uberon_ID', 'marker', 'marker_type']
        )['CL_celltype'].nunique().reset_index(name='celltype_count')
        
        # Calculate specificity: SPs = 1 / number_of_cell_types
        celltype_counts['SPs_tissue_specific'] = (1 / celltype_counts['celltype_count']).round(2)
        
        # Drop intermediate column
        sps_scores = celltype_counts.drop('celltype_count', axis=1)
        
        self.logger.debug(f"Calculated tissue-specific SPs for {len(sps_scores)} marker-tissue combinations")
        return sps_scores
    
    def merge_scores_with_data(self, df: pd.DataFrame, ecs_global: pd.DataFrame, 
                              ecs_tissue: pd.DataFrame, sps_global: pd.DataFrame, 
                              sps_tissue: pd.DataFrame) -> pd.DataFrame:
        """
        Merge calculated scores with original data.
        
        Args:
            df: Original DataFrame
            ecs_global: Global ECs scores
            ecs_tissue: Tissue-specific ECs scores  
            sps_global: Global SPs scores
            sps_tissue: Tissue-specific SPs scores
            
        Returns:
            DataFrame with scores merged
        """
        result = df.copy()
        
        # Merge global ECs
        merge_cols_ecs_global = ['species', 'CL_celltype', 'CL_ID', 'marker', 'marker_type']
        result = result.merge(ecs_global, on=merge_cols_ecs_global, how='left')
        
        # Merge tissue-specific ECs
        merge_cols_ecs_tissue = ['species', 'Uberon_tissue', 'Uberon_ID', 'CL_celltype', 'CL_ID', 'marker', 'marker_type']
        result = result.merge(ecs_tissue, on=merge_cols_ecs_tissue, how='left')
        
        # Merge global SPs
        merge_cols_sps_global = ['species', 'marker', 'marker_type']
        result = result.merge(sps_global, on=merge_cols_sps_global, how='left')
        
        # Merge tissue-specific SPs
        merge_cols_sps_tissue = ['species', 'Uberon_tissue', 'Uberon_ID', 'marker', 'marker_type']
        result = result.merge(sps_tissue, on=merge_cols_sps_tissue, how='left')
        
        # Fill missing scores with appropriate defaults
        score_columns = ['ECs_global', 'ECs_tissue_specific', 'SPs_global', 'SPs_tissue_specific']
        for col in score_columns:
            if col in result.columns:
                if 'ECs' in col:
                    result[col] = result[col].fillna(1)  # Default ECs to 1
                else:  # SPs columns
                    result[col] = result[col].fillna(1.00)  # Default SPs to 1.00
                    # Ensure SPs columns are formatted to 2 decimal places
                    result[col] = result[col].round(2)
        
        self.logger.debug(f"Merged scores with {len(result)} records")
        return result
    
    def add_scores_to_dataframe(self, df: pd.DataFrame, species_col: str = 'species') -> pd.DataFrame:
        """
        Main function to add all scores to a DataFrame.
        
        Args:
            df: Input DataFrame with marker data
            species_col: Name of species column (default: 'species')
            
        Returns:
            DataFrame with calculated ECs and SPs scores added
        """
        if species_col not in df.columns:
            # Add species column if missing, default to 'Human'
            self.logger.warning(f"Species column '{species_col}' not found, defaulting to 'Human'")
            df = df.copy()
            df['species'] = 'Human'
        
        self.logger.debug(f"Starting score calculation for {len(df)} records")
        
        # Calculate all scores
        ecs_global = self.calculate_ecs_global(df)
        ecs_tissue = self.calculate_ecs_tissue_specific(df)
        sps_global = self.calculate_sps_global(df)
        sps_tissue = self.calculate_sps_tissue_specific(df)
        
        # Merge scores with original data
        result = self.merge_scores_with_data(df, ecs_global, ecs_tissue, sps_global, sps_tissue)
        
        self.logger.debug("Score calculation completed successfully")
        return result
    
    def validate_scores_against_reference(self, calculated_df: pd.DataFrame, 
                                        reference_df: pd.DataFrame,
                                        tolerance: float = 0.01) -> Dict[str, bool]:
        """
        Validate calculated scores against reference implementation (e.g., R output).
        
        Args:
            calculated_df: DataFrame with calculated scores
            reference_df: DataFrame with reference scores
            tolerance: Tolerance for numerical comparison
            
        Returns:
            Dictionary with validation results for each score type
        """
        score_columns = ['ECs_global', 'ECs_tissue_specific', 'SPs_global', 'SPs_tissue_specific']
        results = {}
        
        # Merge dataframes for comparison
        merge_cols = ['species', 'Uberon_tissue', 'Uberon_ID', 'CL_celltype', 'CL_ID', 'marker', 'marker_type']
        common_cols = [col for col in merge_cols if col in calculated_df.columns and col in reference_df.columns]
        
        merged = calculated_df.merge(reference_df, on=common_cols, suffixes=('_calc', '_ref'))
        
        for score_col in score_columns:
            calc_col = f"{score_col}_calc"
            ref_col = f"{score_col}_ref"
            
            if calc_col in merged.columns and ref_col in merged.columns:
                # Handle missing values
                valid_rows = merged[[calc_col, ref_col]].notna().all(axis=1)
                if valid_rows.sum() > 0:
                    diff = np.abs(merged.loc[valid_rows, calc_col] - merged.loc[valid_rows, ref_col])
                    results[score_col] = (diff <= tolerance).all()
                    
                    if not results[score_col]:
                        max_diff = diff.max()
                        self.logger.warning(f"Validation failed for {score_col}: max difference = {max_diff:.4f}")
                else:
                    results[score_col] = True  # No valid comparisons
                    
        self.logger.debug(f"Validation results: {results}")
        return results


def main():
    """Example usage of ScoreCalculator."""
    # This would typically be called from convert_original_to_light.py
    import sys
    if len(sys.argv) < 2:
        print("Usage: python score_calculator.py <input_csv> [output_csv]")
        return
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Load data
    df = pd.read_csv(input_file)
    
    # Calculate scores
    calculator = ScoreCalculator()
    result = calculator.add_scores_to_dataframe(df)
    
    # Save result
    if output_file:
        result.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("Results:")
        print(result.head())


if __name__ == "__main__":
    main()