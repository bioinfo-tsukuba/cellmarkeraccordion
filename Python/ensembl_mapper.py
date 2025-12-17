#!/usr/bin/env python3

import pandas as pd
from typing import List, Dict, Optional
import logging
from biomart_client import BioMartClient, BioMartAPIException


class EnsemblMapper:
    """Handle Ensembl ID mapping operations using BioMart."""
    
    def __init__(self, biomart_client: Optional[BioMartClient] = None, log_level: str = 'INFO'):
        """Initialize EnsemblMapper."""
        self.biomart_client = biomart_client or BioMartClient(log_level)
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
        self.mapping_cache = {}
    
    def map_symbols_to_ensembl(self, symbols: List[str], job_config: Dict) -> pd.DataFrame:
        """Map gene symbols to Ensembl IDs using job configuration."""
        if not symbols:
            return pd.DataFrame(columns=['symbol', 'ensembl_gene_id'])
        
        # Create cache key
        cache_key = (
            tuple(sorted(symbols)),
            job_config['release'],
            job_config['ensembl_dataset'],
            job_config['symbol_attribute']
        )
        
        # Check cache
        if cache_key in self.mapping_cache:
            self.logger.debug(f"Using cached mapping for {len(symbols)} symbols")
            return self.mapping_cache[cache_key]
        
        try:
            mapping = self.biomart_client.query_symbol_to_ensembl(
                release=job_config['release'],
                dataset=job_config['ensembl_dataset'],
                symbol_attribute=job_config['symbol_attribute'],
                symbols=symbols
            )
            
            # Cache the result
            self.mapping_cache[cache_key] = mapping
            return mapping
            
        except BioMartAPIException as e:
            self.logger.error(f"Failed to map symbols to Ensembl IDs: {str(e)}")
            return pd.DataFrame(columns=['symbol', 'ensembl_gene_id'])
    
    def map_ensembl_to_orthologs(self, ensembl_ids: List[str], source_config: Dict, 
                                target_config: Dict) -> pd.DataFrame:
        """Map Ensembl IDs to orthologs in target species."""
        if not ensembl_ids:
            return pd.DataFrame(columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
        
        # Create cache key
        cache_key = (
            tuple(sorted(ensembl_ids)),
            source_config['ensembl_dataset'],
            target_config['ensembl_dataset'],
            target_config['release']
        )
        
        # Check cache
        if cache_key in self.mapping_cache:
            self.logger.debug(f"Using cached ortholog mapping for {len(ensembl_ids)} IDs")
            return self.mapping_cache[cache_key]
        
        try:
            mapping = self.biomart_client.query_ortholog_mapping(
                release=target_config['release'],
                source_dataset=source_config['ensembl_dataset'],
                source_ensembl_ids=ensembl_ids,
                target_dataset=target_config['ensembl_dataset']
            )
            
            # Cache the result
            self.mapping_cache[cache_key] = mapping
            return mapping
            
        except BioMartAPIException as e:
            self.logger.error(f"Failed to map orthologs: {str(e)}")
            return pd.DataFrame(columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
    
    def join_data_with_ensembl_mapping(self, data: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
        """Join data with Ensembl mapping and handle multiple mappings."""
        if mapping.empty:
            data_copy = data.copy()
            data_copy['ensembl_gene_id'] = ''
            return data_copy
        
        # Join data with mapping on marker symbol
        result = data.merge(
            mapping,
            left_on='marker',
            right_on='symbol',
            how='left'
        )
        
        # Fill missing ensembl_gene_id with empty string
        result['ensembl_gene_id'] = result['ensembl_gene_id'].fillna('')
        
        # Remove the duplicate symbol column from mapping
        if 'symbol' in result.columns:
            result = result.drop('symbol', axis=1)
        
        return result
    
    def join_data_with_ortholog_mapping(self, data: pd.DataFrame, ortholog_mapping: pd.DataFrame) -> pd.DataFrame:
        """Join data with ortholog mapping."""
        if ortholog_mapping.empty:
            data_copy = data.copy()
            data_copy['ensembl_gene_id'] = ''
            return data_copy
        
        # Join source data with ortholog mapping
        result = data.merge(
            ortholog_mapping,
            left_on='ensembl_gene_id',
            right_on='source_ensembl_id',
            how='left'
        )
        
        # Use target Ensembl ID as the main ensembl_gene_id
        result['ensembl_gene_id'] = result.get('target_ensembl_id', '').fillna('')
        
        # Replace marker symbols with target species symbols where available
        if 'target_symbol' in result.columns:
            # Use target symbol when available, keep original marker when not
            result['marker'] = result['target_symbol'].fillna(result['marker'])
        
        # Clean up temporary columns
        columns_to_drop = ['source_ensembl_id', 'target_ensembl_id', 'target_symbol']
        result = result.drop([col for col in columns_to_drop if col in result.columns], axis=1)
        
        return result
    
    def infer_species_from_symbols(self, symbols: List[str]) -> str:
        """Infer species from gene symbol format."""
        if not symbols:
            return 'Human'
        
        # Simple heuristic: Human genes are typically all uppercase, mouse genes are title case
        sample_symbols = symbols[:10]  # Sample first 10 symbols
        if all(symbol.isupper() for symbol in sample_symbols if symbol):
            return 'Human'
        else:
            return 'Mouse'
    
    def clear_cache(self):
        """Clear mapping cache."""
        self.mapping_cache.clear()
        self.logger.debug("Mapping cache cleared")