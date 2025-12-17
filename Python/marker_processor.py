#!/usr/bin/env python3

import pandas as pd
import os
from typing import List, Dict, Optional
import logging
from marker_database_integration import ScoreCalculator
from ensembl_mapper import EnsemblMapper


class MarkerProcessor:
    """Handle marker file processing operations."""
    
    def __init__(self, ensembl_mapper: Optional[EnsemblMapper] = None, 
                 calculate_scores: bool = False, log_level: str = 'INFO'):
        """Initialize MarkerProcessor."""
        self.ensembl_mapper = ensembl_mapper or EnsemblMapper(log_level=log_level)
        self.calculate_scores = calculate_scores
        self.score_calculator = ScoreCalculator(log_level) if calculate_scores else None
        
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
        
        # Define columns for light version
        self.light_columns = [
            'Uberon_tissue',
            'Uberon_ID', 
            'CL_celltype',
            'CL_ID',
            'marker',
            'marker_type',
            'resource'
        ]
        
        # Add score columns if calculating scores
        if self.calculate_scores:
            self.light_columns.extend([
                'ECs_global',
                'ECs_tissue_specific', 
                'SPs_global',
                'SPs_tissue_specific'
            ])
    
    def load_marker_data(self, file_path: str) -> pd.DataFrame:
        """Load marker data from CSV file."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Marker file not found: {file_path}")
        
        try:
            data = pd.read_csv(file_path)
            self.logger.debug(f"Loaded {len(data)} records from {file_path}")
            return data
        except Exception as e:
            raise RuntimeError(f"Failed to load marker data from {file_path}: {str(e)}")
    
    def get_unique_markers(self, data: pd.DataFrame) -> List[str]:
        """Extract unique marker symbols from data."""
        if 'marker' not in data.columns:
            self.logger.warning("No 'marker' column found in data")
            return []
        
        unique_markers = data['marker'].dropna().unique().tolist()
        self.logger.debug(f"Found {len(unique_markers)} unique markers")
        return unique_markers
    
    def add_ensembl_ids_to_source_data(self, data: pd.DataFrame, job_config: Dict) -> pd.DataFrame:
        """Add Ensembl IDs to source species data."""
        # Get unique markers
        unique_markers = self.get_unique_markers(data)
        
        if not unique_markers:
            self.logger.warning(f"No markers found for {job_config['species']}")
            data_copy = data.copy()
            data_copy['ensembl_gene_id'] = ''
            return data_copy
        
        # Map symbols to Ensembl IDs
        mapping = self.ensembl_mapper.map_symbols_to_ensembl(unique_markers, job_config)
        
        # Join data with mapping
        result = self.ensembl_mapper.join_data_with_ensembl_mapping(data, mapping)
        
        self.logger.debug(f"Added Ensembl IDs to {len(result)} records for {job_config['species']}")
        return result
    
    def add_ensembl_ids_to_derived_data(self, source_data: pd.DataFrame, 
                                       source_config: Dict, target_config: Dict) -> pd.DataFrame:
        """Add Ensembl IDs to derived species data using ortholog mapping."""
        # Get unique source Ensembl IDs
        source_ensembl_ids = source_data['ensembl_gene_id'].dropna()
        source_ensembl_ids = source_ensembl_ids[source_ensembl_ids != ''].unique().tolist()
        
        if not source_ensembl_ids:
            self.logger.warning(f"No source Ensembl IDs found for {target_config['species']} ortholog mapping")
            result = source_data.copy()
            result['ensembl_gene_id'] = ''
            return result
        
        # Map Ensembl IDs to orthologs
        ortholog_mapping = self.ensembl_mapper.map_ensembl_to_orthologs(
            source_ensembl_ids, source_config, target_config)
        
        # Join data with ortholog mapping
        result = self.ensembl_mapper.join_data_with_ortholog_mapping(source_data, ortholog_mapping)
        
        self.logger.debug(f"Added ortholog Ensembl IDs to {len(result)} records for {target_config['species']}")
        return result
    
    def add_species_column_if_missing(self, data: pd.DataFrame, species: str = None) -> pd.DataFrame:
        """Add species column if missing."""
        if 'species' not in data.columns:
            if species is None:
                # Try to infer species from marker symbols
                unique_markers = self.get_unique_markers(data)
                species = self.ensembl_mapper.infer_species_from_symbols(unique_markers)
                self.logger.debug(f"Inferred species: {species}")
            
            data_copy = data.copy()
            data_copy['species'] = species
            return data_copy
        
        return data
    
    def add_cl_ids_if_missing(self, data: pd.DataFrame) -> pd.DataFrame:
        """Add CL_ID column if missing."""
        if 'CL_ID' not in data.columns and 'CL_celltype' in data.columns:
            data_copy = data.copy()
            unique_celltypes = data_copy['CL_celltype'].unique()
            cl_id_map = {celltype: f"CL:{str(i).zfill(7)}" 
                        for i, celltype in enumerate(unique_celltypes, start=1)}
            data_copy['CL_ID'] = data_copy['CL_celltype'].map(cl_id_map)
            self.logger.debug("Added CL_ID mappings for cell types")
            return data_copy
        
        return data
    
    def calculate_and_add_scores(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate and add ECs/SPs scores to data."""
        if not self.calculate_scores or self.score_calculator is None:
            return data
        
        self.logger.debug("Calculating ECs and SPs scores...")
        
        # Check required columns
        required_cols = ['CL_celltype', 'marker', 'marker_type', 'resource']
        missing_cols = [col for col in required_cols if col not in data.columns]
        
        if missing_cols:
            self.logger.warning(f"Missing required columns for score calculation: {missing_cols}")
            return data
        
        # Add missing columns if needed
        data_with_species = self.add_species_column_if_missing(data)
        data_with_cl_ids = self.add_cl_ids_if_missing(data_with_species)
        
        try:
            result = self.score_calculator.add_scores_to_dataframe(data_with_cl_ids)
            self.logger.debug(f"Successfully calculated scores for {len(result)} records")
            
            # Clean up score column names from pandas merge operations
            score_mapping = {}
            for col in result.columns:
                if col.endswith('_y') and any(score in col for score in 
                    ['ECs_global', 'ECs_tissue_specific', 'SPs_global', 'SPs_tissue_specific']):
                    clean_name = col.replace('_y', '')
                    score_mapping[col] = clean_name
            
            if score_mapping:
                result = result.rename(columns=score_mapping)
                # Drop the _x versions if they exist
                x_columns = [col for col in result.columns if col.endswith('_x') and any(score in col for score in 
                    ['ECs_global', 'ECs_tissue_specific', 'SPs_global', 'SPs_tissue_specific'])]
                if x_columns:
                    result = result.drop(columns=x_columns)
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error calculating scores: {str(e)}")
            return data
    
    def filter_to_light_columns(self, data: pd.DataFrame) -> pd.DataFrame:
        """Filter DataFrame to keep only light version columns."""
        # Add ensembl_gene_id to columns to keep
        columns_to_keep = self.light_columns + ['ensembl_gene_id']
        
        # Only keep columns that exist in the data
        existing_columns = [col for col in columns_to_keep if col in data.columns]
        
        result = data[existing_columns]
        self.logger.debug(f"Filtered to {len(existing_columns)} columns: {existing_columns}")
        return result
    
    def save_processed_data(self, data: pd.DataFrame, output_path: str):
        """Save processed data to CSV file."""
        try:
            # Create output directory if it doesn't exist
            output_dir = os.path.dirname(output_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            
            data.to_csv(output_path, index=False)
            self.logger.debug(f"Saved {len(data)} records to {output_path}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to save data to {output_path}: {str(e)}")
    
    def process_source_file(self, job_config: Dict) -> pd.DataFrame:
        """Process a source species file end-to-end."""
        file_path = os.path.join("v1.0.0", job_config['file'])
        
        # Load data
        data = self.load_marker_data(file_path)
        
        # Add Ensembl IDs
        data_with_ensembl = self.add_ensembl_ids_to_source_data(data, job_config)
        
        # Calculate scores if requested
        data_with_scores = self.calculate_and_add_scores(data_with_ensembl)
        
        # Filter to light columns
        result = self.filter_to_light_columns(data_with_scores)
        
        self.logger.debug(f"Processed source file for {job_config['species']}: {len(result)} records")
        return result
    
    def process_derived_file(self, source_data: pd.DataFrame, source_config: Dict, 
                           target_config: Dict) -> pd.DataFrame:
        """Process a derived species file using source data."""
        # Add ortholog Ensembl IDs
        data_with_orthologs = self.add_ensembl_ids_to_derived_data(
            source_data, source_config, target_config)
        
        # Calculate scores if requested (scores are calculated on the ortholog data)
        data_with_scores = self.calculate_and_add_scores(data_with_orthologs)
        
        # Filter to light columns
        result = self.filter_to_light_columns(data_with_scores)
        
        self.logger.debug(f"Processed derived file for {target_config['species']}: {len(result)} records")
        return result