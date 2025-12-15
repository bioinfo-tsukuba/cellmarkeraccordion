#!/usr/bin/env python3

import json
import pandas as pd
import requests
import time
import os
from typing import Dict, List, Optional, Tuple
import argparse

class BioMartConverter:
    """Convert original marker files to light format with Ensembl IDs using BioMart."""
    
    def __init__(self, config_file: str = "biomart_config.json"):
        """Initialize converter with configuration."""
        self.config_file = config_file
        self.config = self._load_config()
        self.cache = {}  # Cache for derived species processing
        
        # Columns to keep in light version
        self.light_columns = [
            'Uberon_tissue',
            'Uberon_ID', 
            'CL_celltype',
            'CL_ID',
            'marker',
            'marker_type',
            'resource'
        ]
    
    def _load_config(self) -> Dict:
        """Load biomart configuration from JSON file."""
        with open(self.config_file, 'r') as f:
            return json.load(f)
    
    def _get_biomart_url(self, release: int) -> str:
        """Get BioMart URL for specified Ensembl release."""
        if release == 110:
            return "https://www.ensembl.org/biomart/martservice"
        else:
            return f"https://{release}.archive.ensembl.org/biomart/martservice"
    
    def _query_biomart(self, xml_query: str, release: int, query_type: str = 'symbol_to_ensembl') -> pd.DataFrame:
        """Query BioMart with XML and return results as DataFrame."""
        url = self._get_biomart_url(release)
        
        response = requests.post(
            url,
            data={'query': xml_query},
            headers={'Content-Type': 'application/x-www-form-urlencoded'}
        )
        
        if response.status_code != 200:
            print(f"BioMart query failed with status {response.status_code}")
            return pd.DataFrame()
        
        # Parse tab-separated response
        lines = response.text.strip().split('\n')
        if not lines or lines[0] == '':
            return pd.DataFrame()
        
        # Create DataFrame from response
        data = []
        for line in lines:
            if line.strip():
                data.append(line.split('\t'))
        
        if not data:
            return pd.DataFrame()
        
        # Set column names based on query type
        if query_type == 'symbol_to_ensembl':
            df = pd.DataFrame(data, columns=['symbol', 'ensembl_gene_id'])
            # Remove empty results
            df = df[df['ensembl_gene_id'].str.strip() != '']
        elif query_type == 'ortholog':
            # For ortholog queries, we expect 3 columns
            df = pd.DataFrame(data)
            if len(df.columns) >= 3:
                df.columns = ['source_ensembl_id', 'target_ensembl_id', 'target_symbol']
                # Remove empty target results
                df = df[df['target_ensembl_id'].str.strip() != '']
            else:
                # Handle case where we get fewer columns than expected
                df = pd.DataFrame(columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
        
        return df
    
    def biomart_symbol_to_ensembl(self, release: int, dataset: str, symbol_attribute: str, symbols: List[str]) -> pd.DataFrame:
        """Convert gene symbols to Ensembl Gene IDs using BioMart."""
        
        # Create XML query for symbol to Ensembl ID mapping
        symbols_filter = ','.join(symbols)
        
        xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        <Filter name="{symbol_attribute}" value="{symbols_filter}"/>
        <Attribute name="{symbol_attribute}"/>
        <Attribute name="ensembl_gene_id"/>
    </Dataset>
</Query>"""
        
        return self._query_biomart(xml_query, release, 'symbol_to_ensembl')
    
    def biomart_ortholog_mapping(self, release: int, source_ensembl_ids: List[str], target_dataset: str, target_symbol_attribute: str) -> pd.DataFrame:
        """Convert source Ensembl IDs to target species orthologs using BioMart."""
        
        # For ortholog mapping, we need to determine the source dataset from the first ID
        # This is a simplified approach - in reality you'd need to track the source dataset
        source_dataset = "mmusculus_gene_ensembl"  # Assuming mouse as source for rat
        
        ids_filter = ','.join(source_ensembl_ids)
        
        xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{source_dataset}" interface="default">
        <Filter name="ensembl_gene_id" value="{ids_filter}"/>
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="{target_dataset.replace('_gene_ensembl', '')}_homolog_ensembl_gene"/>
        <Attribute name="{target_dataset.replace('_gene_ensembl', '')}_homolog_associated_gene_name"/>
    </Dataset>
</Query>"""
        
        return self._query_biomart(xml_query, release, 'ortholog')
    
    def join_and_expand(self, data: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
        """Join data with mapping and expand rows for multiple Ensembl IDs."""
        if mapping.empty:
            # If no mapping found, add empty ensembl_gene_id column
            data['ensembl_gene_id'] = ''
            return data
        
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
    
    def filter_columns(self, data: pd.DataFrame) -> pd.DataFrame:
        """Filter DataFrame to keep only light version columns."""
        # Add ensembl_gene_id to the columns to keep
        columns_to_keep = self.light_columns + ['ensembl_gene_id']
        
        # Only keep columns that exist in the data
        existing_columns = [col for col in columns_to_keep if col in data.columns]
        
        return data[existing_columns]
    
    def process_source_species(self, job_config: Dict) -> pd.DataFrame:
        """Process source species (Human/Mouse) files."""
        file_path = os.path.join("v1.0.0", job_config['file'])
        
        print(f"Loading {job_config['species']} data from {file_path}")
        data = pd.read_csv(file_path)
        
        # Filter to light columns first
        data = self.filter_columns(data)
        
        # Get unique markers for BioMart query
        unique_markers = data['marker'].dropna().unique().tolist()
        
        if not unique_markers:
            print(f"No markers found in {job_config['species']} data")
            data['ensembl_gene_id'] = ''
            return data
        
        print(f"Querying BioMart for {len(unique_markers)} unique {job_config['species']} markers")
        
        # Query BioMart for symbol to Ensembl mapping
        mapping = self.biomart_symbol_to_ensembl(
            release=job_config['release'],
            dataset=job_config['ensembl_dataset'],
            symbol_attribute=job_config['symbol_attribute'],
            symbols=unique_markers
        )
        
        print(f"Found {len(mapping)} symbol-to-Ensembl mappings for {job_config['species']}")
        
        # Join and expand data with mapping
        result = self.join_and_expand(data, mapping)
        
        return result
    
    def process_derived_species(self, job_config: Dict) -> pd.DataFrame:
        """Process derived species (Rat) using ortholog mapping from source species."""
        source_species = job_config['source_species']
        
        if source_species not in self.cache:
            raise ValueError(f"Source species {source_species} not found in cache")
        
        source_data = self.cache[source_species]
        
        # Get unique source Ensembl IDs
        source_ensembl_ids = source_data['ensembl_gene_id'].dropna()
        source_ensembl_ids = source_ensembl_ids[source_ensembl_ids != ''].unique().tolist()
        
        if not source_ensembl_ids:
            print(f"No source Ensembl IDs found for {job_config['species']} ortholog mapping")
            # Return empty DataFrame with correct columns
            result = source_data.copy()
            result['ensembl_gene_id'] = ''
            return result
        
        print(f"Querying BioMart for {len(source_ensembl_ids)} ortholog mappings to {job_config['species']}")
        
        # Query BioMart for ortholog mapping
        ortholog_mapping = self.biomart_ortholog_mapping(
            release=job_config['release'],
            source_ensembl_ids=source_ensembl_ids,
            target_dataset=job_config['ensembl_dataset'],
            target_symbol_attribute=job_config['symbol_attribute']
        )
        
        print(f"Found {len(ortholog_mapping)} ortholog mappings for {job_config['species']}")
        
        # Join source data with ortholog mapping
        if not ortholog_mapping.empty:
            result = source_data.merge(
                ortholog_mapping,
                left_on='ensembl_gene_id',
                right_on='source_ensembl_id',
                how='left'
            )
            
            # Use target Ensembl ID as the main ensembl_gene_id
            result['ensembl_gene_id'] = result.get('target_ensembl_id', '').fillna('')
            
            # Clean up temporary columns
            columns_to_drop = ['source_ensembl_id', 'target_ensembl_id', 'target_symbol']
            result = result.drop([col for col in columns_to_drop if col in result.columns], axis=1)
        else:
            result = source_data.copy()
            result['ensembl_gene_id'] = ''
        
        return result
    
    def add_ensembl_id(self, config_id: int = 1) -> None:
        """Main pipeline to add Ensembl IDs to marker files."""
        
        # Find config by ID
        config_set = None
        for marker_set in self.config['marker_list_set']:
            if marker_set['id'] == config_id:
                config_set = marker_set
                break
        
        if not config_set:
            raise ValueError(f"Config ID {config_id} not found")
        
        # Process each version block
        for block in config_set['list']:
            print(f"Processing version {block['version']}")
            
            # Process biomart jobs
            for job in block['biomart_configs']:
                
                if job['role'] == 'source':
                    print(f"\nProcessing source species: {job['species']}")
                    
                    # Process source file
                    result = self.process_source_species(job)
                    
                    # Cache result for derived species
                    self.cache[job['species']] = result
                    
                    # Save light file
                    output_file = job['file'].replace('original_', 'light_')
                    output_path = os.path.join("v1.0.0", output_file)
                    result.to_csv(output_path, index=False)
                    print(f"Saved {output_path}")
                
                elif job['role'] == 'derived':
                    print(f"\nProcessing derived species: {job['species']}")
                    
                    # Process derived file using ortholog mapping
                    result = self.process_derived_species(job)
                    
                    # Cache result
                    self.cache[job['species']] = result
                    
                    # Save light file
                    species_name = job['species'].lower()
                    output_file = f"light_{species_name}_healthy.csv"
                    output_path = os.path.join("v1.0.0", output_file)
                    result.to_csv(output_path, index=False)
                    print(f"Saved {output_path}")
        
        print("\nConversion completed!")

def main():
    """Main function to run the converter."""
    parser = argparse.ArgumentParser(description='Convert original marker files to light format with Ensembl IDs')
    parser.add_argument('--config-id', type=int, default=1, help='Configuration ID to use (default: 1)')
    parser.add_argument('--config-file', type=str, default='biomart_config.json', help='Configuration file (default: biomart_config.json)')
    
    args = parser.parse_args()
    
    # Create converter and run pipeline
    converter = BioMartConverter(config_file=args.config_file)
    converter.add_ensembl_id(config_id=args.config_id)

if __name__ == "__main__":
    main()