#!/usr/bin/env python3

import json
import os
from typing import Dict, List, Optional
import logging


class ConfigurationError(Exception):
    """Exception raised for configuration-related errors."""
    pass


class ConfigManager:
    """Manage BioMart configuration loading and validation."""
    
    def __init__(self, config_file: str, log_level: str = 'INFO'):
        """Initialize ConfigManager with configuration file."""
        self.config_file = config_file
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
        self._config = None
    
    def load_config(self) -> Dict:
        """Load and validate configuration from JSON file."""
        if self._config is not None:
            return self._config
        
        if not os.path.exists(self.config_file):
            raise ConfigurationError(f"Configuration file not found: {self.config_file}")
        
        try:
            with open(self.config_file, 'r') as f:
                self._config = json.load(f)
            
            self._validate_config()
            self.logger.debug(f"Successfully loaded configuration from {self.config_file}")
            return self._config
            
        except json.JSONDecodeError as e:
            raise ConfigurationError(f"Invalid JSON in configuration file: {str(e)}")
        except Exception as e:
            raise ConfigurationError(f"Failed to load configuration: {str(e)}")
    
    def _validate_config(self):
        """Validate configuration structure."""
        if not isinstance(self._config, dict):
            raise ConfigurationError("Configuration must be a JSON object")
        
        if 'marker_list_set' not in self._config:
            raise ConfigurationError("Configuration missing 'marker_list_set' key")
        
        if not isinstance(self._config['marker_list_set'], list):
            raise ConfigurationError("'marker_list_set' must be an array")
        
        for i, marker_set in enumerate(self._config['marker_list_set']):
            self._validate_marker_set(marker_set, i)
    
    def _validate_marker_set(self, marker_set: Dict, index: int):
        """Validate individual marker set configuration."""
        required_keys = ['id', 'list']
        for key in required_keys:
            if key not in marker_set:
                raise ConfigurationError(f"Marker set {index} missing required key: {key}")
        
        if not isinstance(marker_set['list'], list):
            raise ConfigurationError(f"Marker set {index} 'list' must be an array")
        
        for j, version_block in enumerate(marker_set['list']):
            self._validate_version_block(version_block, index, j)
    
    def _validate_version_block(self, version_block: Dict, set_index: int, block_index: int):
        """Validate version block configuration."""
        required_keys = ['version', 'biomart_configs']
        for key in required_keys:
            if key not in version_block:
                raise ConfigurationError(f"Version block {set_index}.{block_index} missing required key: {key}")
        
        if not isinstance(version_block['biomart_configs'], list):
            raise ConfigurationError(f"Version block {set_index}.{block_index} 'biomart_configs' must be an array")
        
        for k, job_config in enumerate(version_block['biomart_configs']):
            self._validate_job_config(job_config, set_index, block_index, k)
    
    def _validate_job_config(self, job_config: Dict, set_index: int, block_index: int, job_index: int):
        """Validate individual job configuration."""
        location = f"Job {set_index}.{block_index}.{job_index}"
        
        required_keys = ['role', 'species', 'tax_id', 'ensembl_dataset', 'symbol_attribute', 'backend', 'release']
        for key in required_keys:
            if key not in job_config:
                raise ConfigurationError(f"{location} missing required key: {key}")
        
        # Validate role
        if job_config['role'] not in ['source', 'derived']:
            raise ConfigurationError(f"{location} invalid role: {job_config['role']}")
        
        # Source jobs need file attribute
        if job_config['role'] == 'source' and 'file' not in job_config:
            raise ConfigurationError(f"{location} source jobs must have 'file' attribute")
        
        # Derived jobs need source_species attribute
        if job_config['role'] == 'derived' and 'source_species' not in job_config:
            raise ConfigurationError(f"{location} derived jobs must have 'source_species' attribute")
    
    def get_config_by_id(self, config_id: int) -> Dict:
        """Get configuration set by ID."""
        config = self.load_config()
        
        for marker_set in config['marker_list_set']:
            if marker_set['id'] == config_id:
                return marker_set
        
        raise ConfigurationError(f"Configuration ID {config_id} not found")
    
    def get_source_jobs(self, config_id: int) -> List[Dict]:
        """Get all source jobs for a configuration ID."""
        config_set = self.get_config_by_id(config_id)
        source_jobs = []
        
        for version_block in config_set['list']:
            for job_config in version_block['biomart_configs']:
                if job_config['role'] == 'source':
                    source_jobs.append(job_config)
        
        return source_jobs
    
    def get_derived_jobs(self, config_id: int) -> List[Dict]:
        """Get all derived jobs for a configuration ID."""
        config_set = self.get_config_by_id(config_id)
        derived_jobs = []
        
        for version_block in config_set['list']:
            for job_config in version_block['biomart_configs']:
                if job_config['role'] == 'derived':
                    derived_jobs.append(job_config)
        
        return derived_jobs
    
    def get_job_by_species(self, config_id: int, species: str) -> Optional[Dict]:
        """Get job configuration for specific species."""
        config_set = self.get_config_by_id(config_id)
        
        for version_block in config_set['list']:
            for job_config in version_block['biomart_configs']:
                if job_config['species'] == species:
                    return job_config
        
        return None
    
    def list_available_species(self, config_id: int) -> List[str]:
        """List all available species in configuration."""
        config_set = self.get_config_by_id(config_id)
        species = []
        
        for version_block in config_set['list']:
            for job_config in version_block['biomart_configs']:
                if job_config['species'] not in species:
                    species.append(job_config['species'])
        
        return species
    
    def get_version_blocks(self, config_id: int) -> List[Dict]:
        """Get all version blocks for a configuration ID."""
        config_set = self.get_config_by_id(config_id)
        return config_set['list']