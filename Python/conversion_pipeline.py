#!/usr/bin/env python3

import os
from typing import Dict, Optional
import logging
from .config_manager import ConfigManager, ConfigurationError
from .marker_processor import MarkerProcessor
from .ensembl_mapper import EnsemblMapper
from .biomart_client import BioMartClient, BioMartAPIException


class ConversionPipeline:
    """Orchestrate the marker file conversion pipeline."""
    
    def __init__(self, config_file: str = "biomart_config.json", 
                 calculate_scores: bool = False, log_level: str = 'INFO'):
        """Initialize ConversionPipeline."""
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
        
        # Initialize components
        self.config_manager = ConfigManager(config_file, log_level)
        self.biomart_client = BioMartClient(log_level)
        self.ensembl_mapper = EnsemblMapper(self.biomart_client, log_level)
        self.marker_processor = MarkerProcessor(self.ensembl_mapper, calculate_scores, log_level)
        
        # Cache for processed source species data
        self.source_data_cache = {}
        self.calculate_scores = calculate_scores
    
    def run_conversion(self, config_id: int = 1):
        """Run the complete conversion pipeline for a configuration ID."""
        try:
            self.logger.info(f"Converting markers for config ID {config_id}")
            
            # Get configuration
            config_set = self.config_manager.get_config_by_id(config_id)
            
            # Process each version block
            for block in config_set['list']:
                self.logger.debug(f"Processing version {block['version']}")
                self._process_version_block(block)
            
            self.logger.info("Conversion completed!")
            
        except (ConfigurationError, BioMartAPIException) as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in pipeline: {str(e)}")
            raise RuntimeError(f"Pipeline execution failed: {str(e)}")
    
    def _process_version_block(self, version_block: Dict):
        """Process a single version block."""
        biomart_configs = version_block['biomart_configs']
        
        # First, process all source species
        for job_config in biomart_configs:
            if job_config['role'] == 'source':
                self._process_source_species(job_config)
        
        # Then, process all derived species
        for job_config in biomart_configs:
            if job_config['role'] == 'derived':
                self._process_derived_species(job_config)
    
    def _process_source_species(self, job_config: Dict):
        """Process a source species job."""
        species = job_config['species']
        self.logger.info(f"Processing {species}...")
        
        try:
            # Process the source file
            result = self.marker_processor.process_source_file(job_config)
            
            # Cache the result for derived species
            self.source_data_cache[species] = result
            
            # Save the light file
            output_file = self._get_output_filename(job_config)
            output_path = os.path.join("v1.0.0", output_file)
            self.marker_processor.save_processed_data(result, output_path)
            
            self.logger.info(f"Saved {species}: {len(result)} records")
            
        except Exception as e:
            self.logger.error(f"Failed to process source species {species}: {str(e)}")
            raise
    
    def _process_derived_species(self, job_config: Dict):
        """Process a derived species job."""
        species = job_config['species']
        source_species = job_config['source_species']
        
        self.logger.info(f"Processing {species} (from {source_species})...")
        
        # Check if source data is available
        if source_species not in self.source_data_cache:
            raise RuntimeError(f"Source species {source_species} not found in cache for {species}")
        
        try:
            # Get source data and config
            source_data = self.source_data_cache[source_species]
            source_config = self._get_source_config_for_species(source_species)
            
            # Process the derived file
            result = self.marker_processor.process_derived_file(
                source_data, source_config, job_config)
            
            # Cache the result
            self.source_data_cache[species] = result
            
            # Save the light file
            output_file = self._get_derived_output_filename(job_config)
            output_path = os.path.join("v1.0.0", output_file)
            self.marker_processor.save_processed_data(result, output_path)
            
            self.logger.info(f"Saved {species}: {len(result)} records")
            
        except Exception as e:
            self.logger.error(f"Failed to process derived species {species}: {str(e)}")
            raise
    
    def _get_source_config_for_species(self, species: str) -> Dict:
        """Get source configuration for a species."""
        # This is a simplified approach - in a more complete implementation,
        # we'd need to track the source configs used
        if species == "Human":
            return {
                'ensembl_dataset': 'hsapiens_gene_ensembl',
                'symbol_attribute': 'hgnc_symbol',
                'release': 110
            }
        elif species == "Mouse":
            return {
                'ensembl_dataset': 'mmusculus_gene_ensembl',
                'symbol_attribute': 'mgi_symbol',
                'release': 110
            }
        else:
            raise ValueError(f"Unknown source species: {species}")
    
    def _get_output_filename(self, job_config: Dict) -> str:
        """Generate output filename for source species."""
        original_file = job_config['file']
        return original_file.replace('original_', 'light_')
    
    def _get_derived_output_filename(self, job_config: Dict) -> str:
        """Generate output filename for derived species."""
        return job_config['outfile']
    
    def process_single_species(self, config_id: int, species: str):
        """Process a single species from the configuration."""
        try:
            # Get job config for the species
            job_config = self.config_manager.get_job_by_species(config_id, species)
            if not job_config:
                raise ValueError(f"Species {species} not found in configuration {config_id}")
            
            self.logger.info(f"Processing {species}...")
            
            if job_config['role'] == 'source':
                self._process_source_species(job_config)
            elif job_config['role'] == 'derived':
                # Need to process source species first
                source_species = job_config['source_species']
                if source_species not in self.source_data_cache:
                    source_config = self.config_manager.get_job_by_species(config_id, source_species)
                    if source_config:
                        self._process_source_species(source_config)
                    else:
                        raise ValueError(f"Source species {source_species} not found in configuration")
                
                self._process_derived_species(job_config)
            
            self.logger.info(f"Completed {species}")
            
        except Exception as e:
            self.logger.error(f"Failed to process species {species}: {str(e)}")
            raise
    
    def list_available_species(self, config_id: int):
        """List available species in configuration."""
        try:
            species_list = self.config_manager.list_available_species(config_id)
            self.logger.debug(f"Available species in configuration {config_id}: {species_list}")
            return species_list
        except Exception as e:
            self.logger.error(f"Failed to list species: {str(e)}")
            raise
    
    def clear_cache(self):
        """Clear all cached data."""
        self.source_data_cache.clear()
        self.ensembl_mapper.clear_cache()
        self.logger.debug("Pipeline cache cleared")
    
    def get_pipeline_status(self) -> Dict:
        """Get current pipeline status."""
        return {
            'cached_species': list(self.source_data_cache.keys()),
            'score_calculation_enabled': self.calculate_scores,
            'config_file': self.config_manager.config_file
        }