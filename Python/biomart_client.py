#!/usr/bin/env python3

import requests
import pandas as pd
from typing import List, Dict
import logging


class BioMartAPIException(Exception):
    """Exception raised for BioMart API errors."""
    pass


class BioMartClient:
    """Pure BioMart API client for querying Ensembl data."""
    
    def __init__(self, log_level: str = 'INFO'):
        """Initialize BioMart client."""
        logging.basicConfig(level=getattr(logging, log_level.upper()))
        self.logger = logging.getLogger(__name__)
    
    def get_biomart_url(self, release: int) -> str:
        """Get BioMart URL for specified Ensembl release."""
        if release == 110:
            return "https://www.ensembl.org/biomart/martservice"
        else:
            return f"https://{release}.archive.ensembl.org/biomart/martservice"
    
    def build_symbol_to_ensembl_query(self, dataset: str, symbol_attribute: str, symbols: List[str]) -> str:
        """Build XML query for symbol to Ensembl ID mapping."""
        symbols_filter = ','.join(symbols)
        
        return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        <Filter name="{symbol_attribute}" value="{symbols_filter}"/>
        <Attribute name="{symbol_attribute}"/>
        <Attribute name="ensembl_gene_id"/>
    </Dataset>
</Query>"""
    
    def build_ortholog_query(self, source_dataset: str, source_ensembl_ids: List[str], 
                            target_dataset: str) -> str:
        """Build XML query for ortholog mapping."""
        ids_filter = ','.join(source_ensembl_ids)
        target_prefix = target_dataset.replace('_gene_ensembl', '')
        
        return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="{source_dataset}" interface="default">
        <Filter name="ensembl_gene_id" value="{ids_filter}"/>
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="{target_prefix}_homolog_ensembl_gene"/>
        <Attribute name="{target_prefix}_homolog_associated_gene_name"/>
    </Dataset>
</Query>"""
    
    def query_biomart(self, xml_query: str, release: int) -> requests.Response:
        """Execute BioMart query and return response."""
        url = self.get_biomart_url(release)
        
        try:
            response = requests.post(
                url,
                data={'query': xml_query},
                headers={'Content-Type': 'application/x-www-form-urlencoded'},
                timeout=300  # 5 minutes timeout
            )
            response.raise_for_status()
            return response
            
        except requests.exceptions.RequestException as e:
            raise BioMartAPIException(f"BioMart API request failed: {str(e)}")
    
    def parse_symbol_to_ensembl_response(self, response: requests.Response) -> pd.DataFrame:
        """Parse BioMart response for symbol to Ensembl mapping."""
        lines = response.text.strip().split('\n')
        if not lines or lines[0] == '':
            return pd.DataFrame(columns=['symbol', 'ensembl_gene_id'])
        
        data = []
        for line in lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    data.append(parts[:2])
        
        if not data:
            return pd.DataFrame(columns=['symbol', 'ensembl_gene_id'])
        
        df = pd.DataFrame(data, columns=['symbol', 'ensembl_gene_id'])
        # Remove empty results
        df = df[df['ensembl_gene_id'].str.strip() != '']
        return df
    
    def parse_ortholog_response(self, response: requests.Response) -> pd.DataFrame:
        """Parse BioMart response for ortholog mapping."""
        lines = response.text.strip().split('\n')
        if not lines or lines[0] == '':
            return pd.DataFrame(columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
        
        data = []
        for line in lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 3:
                    data.append(parts[:3])
        
        if not data:
            return pd.DataFrame(columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
        
        df = pd.DataFrame(data, columns=['source_ensembl_id', 'target_ensembl_id', 'target_symbol'])
        # Remove empty target results
        df = df[df['target_ensembl_id'].str.strip() != '']
        return df
    
    def query_symbol_to_ensembl(self, release: int, dataset: str, symbol_attribute: str, 
                               symbols: List[str]) -> pd.DataFrame:
        """Query BioMart for symbol to Ensembl ID mapping."""
        self.logger.debug(f"Querying BioMart for {len(symbols)} symbols using dataset {dataset}")
        
        xml_query = self.build_symbol_to_ensembl_query(dataset, symbol_attribute, symbols)
        response = self.query_biomart(xml_query, release)
        result = self.parse_symbol_to_ensembl_response(response)
        
        self.logger.debug(f"Found {len(result)} symbol-to-Ensembl mappings")
        return result
    
    def query_ortholog_mapping(self, release: int, source_dataset: str, source_ensembl_ids: List[str], 
                              target_dataset: str) -> pd.DataFrame:
        """Query BioMart for ortholog mapping."""
        self.logger.debug(f"Querying BioMart for {len(source_ensembl_ids)} ortholog mappings from {source_dataset} to {target_dataset}")
        
        xml_query = self.build_ortholog_query(source_dataset, source_ensembl_ids, target_dataset)
        response = self.query_biomart(xml_query, release)
        result = self.parse_ortholog_response(response)
        
        self.logger.debug(f"Found {len(result)} ortholog mappings")
        return result