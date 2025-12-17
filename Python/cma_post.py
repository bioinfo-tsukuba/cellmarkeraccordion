#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import sys
from typing import Dict, Optional

# Import CellMarkerAccordion functionality
from marker_database_integration import ScoreCalculator
from marker_processor import MarkerProcessor


def detect_tissue_ontology_id(target_srx):
    """
    Detect tissue ontology ID for CellMarkerAccordion database.
    
    For hippocampus samples, return the hippocampal formation UBERON ID.
    This can be extended to support other tissue types based on SRX metadata.
    """
    # For hippocampus studies, use hippocampal formation
    # This corresponds to the tissue available in CellMarkerAccordion database
    return "UBERON:0002421"  # hippocampal formation


def load_cellmarkeraccordion_database(process_id, tax_id, tissue_ontology=None):
    """
    Load CellMarkerAccordion database based on process_id configuration and tax_id.
    
    Parameters:
        process_id: Process ID from configuration
        tax_id: Taxonomic ID for species selection (e.g., 9606 for human, 10090 for mouse)
        tissue_ontology: Tissue ontology ID (e.g., UBERON:0002421)
    
    Returns:
        dict: Dictionary with cell type markers organized by cell type
    """
    print(f"--- Loading CellMarkerAccordion database for process_id: {process_id}")
    
    # Determine database file path based on configuration
    with open("../post_processing_config.json", 'r') as f:
        config = pd.read_json(f)
    
    # Find configuration for this process_id
    process_config = None
    for proc in config['post_processing']:
        if proc['id'] == str(process_id):
            process_config = proc['configs'][0]  # Take first config
            break
    
    if not process_config:
        raise ValueError(f"No configuration found for process_id: {process_id}")
    
    annotation_config = process_config['Annotation']
    if annotation_config['tool'] != 'cellmarkeraccordion':
        raise ValueError(f"Expected cellmarkeraccordion, got: {annotation_config['tool']}")
    
    # Load biomart_config.json to get version and outfile for this annotation id
    biomart_config_path = f"../tools/cellmarkeraccordion/{annotation_config['config']}"
    with open(biomart_config_path, 'r') as f:
        biomart_config = pd.read_json(f)
    
    # Find the marker list set with the matching annotation id
    marker_set = None
    for marker_list in biomart_config['marker_list_set']:
        if marker_list['id'] == int(annotation_config['id']):
            marker_set = marker_list['list'][0]  # Take first list entry
            break
    
    if not marker_set:
        raise ValueError(f"No marker set found for annotation id: {annotation_config['id']}")
    
    # Get version and find species-specific outfile based on tax_id
    db_version = f"v{marker_set['version']}"
    species_config = None
    for config in marker_set['biomart_configs']:
        if config['tax_id'] == int(tax_id):
            species_config = config
            break
    
    if not species_config:
        raise ValueError(f"No configuration found for tax_id: {tax_id} in biomart_config")
    
    db_file = species_config['outfile']  # e.g., "light_human_healthy.csv", "light_mouse_healthy.csv"
    db_path = f"../tools/cellmarkeraccordion/{db_version}/{db_file}"
    
    print(f"--- Loading CellMarkerAccordion database from: {db_path}")
    
    try:
        # Load CellMarkerAccordion data
        markers_df = pd.read_csv(db_path)
        print(f"--- Loaded {len(markers_df)} marker entries from CellMarkerAccordion database")
        
        # Filter by tissue if specified
        if tissue_ontology:
            tissue_filtered = markers_df[markers_df['Uberon_ID'] == tissue_ontology]
            if len(tissue_filtered) > 0:
                markers_df = tissue_filtered
                print(f"--- Filtered to {len(markers_df)} markers for tissue: {tissue_ontology}")
            else:
                print(f"--- Warning: No markers found for tissue {tissue_ontology}, using all tissues")
        
        # Group markers by cell type for annotation
        cell_type_markers = {}
        ontology_mapping = {}
        
        for cell_type in markers_df['CL_celltype'].unique():
            cell_markers = markers_df[markers_df['CL_celltype'] == cell_type]
            
            # Get positive markers only (CellMarkerAccordion format)
            positive_markers = cell_markers[cell_markers['marker_type'] == 'positive']['marker'].tolist()
            
            if positive_markers:
                cell_type_markers[cell_type] = {
                    'positive': positive_markers,
                    'scores': cell_markers[['marker', 'ECs_global', 'SPs_global']].drop_duplicates(subset=['marker']).set_index('marker').to_dict('index')
                }
                
                # Store ontology mapping (cell_type -> CL_ID)
                cl_id = cell_markers['CL_ID'].iloc[0]  # All rows for same cell_type should have same CL_ID
                ontology_mapping[cell_type] = cl_id
        
        print(f"--- Organized markers for {len(cell_type_markers)} cell types")
        return cell_type_markers, ontology_mapping
        
    except FileNotFoundError:
        # Fallback to v1.0.0 if specific version not found
        db_path_fallback = f"../tools/cellmarkeraccordion/v1.0.0/{db_file}"
        print(f"--- Fallback: Loading from {db_path_fallback}")
        
        markers_df = pd.read_csv(db_path_fallback)
        print(f"--- Loaded {len(markers_df)} marker entries from fallback database")
        
        # Filter by tissue if specified
        if tissue_ontology:
            tissue_filtered = markers_df[markers_df['Uberon_ID'] == tissue_ontology]
            if len(tissue_filtered) > 0:
                markers_df = tissue_filtered
                print(f"--- Filtered to {len(markers_df)} markers for tissue: {tissue_ontology}")
        
        # Group markers by cell type
        cell_type_markers = {}
        ontology_mapping = {}
        
        for cell_type in markers_df['CL_celltype'].unique():
            cell_markers = markers_df[markers_df['CL_celltype'] == cell_type]
            positive_markers = cell_markers[cell_markers['marker_type'] == 'positive']['marker'].tolist()
            
            if positive_markers:
                cell_type_markers[cell_type] = {
                    'positive': positive_markers,
                    'scores': cell_markers[['marker', 'ECs_global', 'SPs_global']].drop_duplicates(subset=['marker']).set_index('marker').to_dict('index')
                }
                
                # Store ontology mapping (cell_type -> CL_ID)
                cl_id = cell_markers['CL_ID'].iloc[0]  # All rows for same cell_type should have same CL_ID
                ontology_mapping[cell_type] = cl_id
        
        return cell_type_markers, ontology_mapping


def perform_cellmarkeraccordion_annotation(adata, cluster_key="louvain", 
                                         process_id=1, tax_id=9606, tissue_ontology=None):
                                         
    print(f"--- Starting CellMarkerAccordion annotation with process_id: {process_id}")
    
    # Load CellMarkerAccordion database
    cell_type_markers, ontology_mapping = load_cellmarkeraccordion_database(process_id, tax_id, tissue_ontology)
    
    # Get all marker genes from CellMarkerAccordion database
    all_marker_genes = set()
    for cell_type_data in cell_type_markers.values():
        all_marker_genes.update(cell_type_data['positive'])
    
    print(f"--- Found {len(all_marker_genes)} unique marker genes in CellMarkerAccordion database")
    
    # Use full expression data if provided, otherwise use adata.raw or current adata
    full_expr = adata.to_df().T  # Fallback to processed data
    
    # Set gene names as index - use gene_symbols if available, otherwise use var_names
    if 'gene_symbols' in adata.var.columns:
        full_expr.index = adata.var['gene_symbols']
        print(f"--- Using gene_symbols for gene names")
    else:
        full_expr.index = adata.var_names
        print(f"--- Warning: gene_symbols not found, using var_names instead")
    
    # Filter to only marker genes present in dataset (intersect)
    available_marker_genes = [gene for gene in all_marker_genes if gene in full_expr.index]
    expr_data = full_expr.loc[available_marker_genes]
    
    print(f"--- Using {len(available_marker_genes)} marker genes present in dataset")
    print(f"--- Expression data shape for CellMarkerAccordion: {expr_data.shape} (marker_genes × cells)")
    
    # Calculate CellMarkerAccordion scores using ECs and SPs weighting
    cluster_results = []
    for cluster in adata.obs[cluster_key].unique():
        # Get cells in this cluster
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_cells = adata.obs.index[cluster_mask]
        
        if len(cluster_cells) == 0:
            continue
            
        # Calculate cluster expression profile (mean expression)
        cluster_expr = expr_data[cluster_cells].mean(axis=1)
        
        # Score each cell type
        cell_type_scores = {}
        for cell_type, markers_info in cell_type_markers.items():
            positive_markers = markers_info['positive']
            marker_scores = markers_info['scores']
            
            # Get intersection of positive markers with available genes
            available_positive = [m for m in positive_markers if m in cluster_expr.index]
            
            if not available_positive:
                cell_type_scores[cell_type] = 0.0
                continue
            
            # Calculate weighted score using ECs and SPs
            total_score = 0.0
            total_weight = 0.0
            
            for marker in available_positive:
                expr_level = cluster_expr[marker]
                
                # Get ECs and SPs scores for this marker (default to 1.0 if not found)
                ecs_score = marker_scores.get(marker, {}).get('ECs_global', 1.0)
                sps_score = marker_scores.get(marker, {}).get('SPs_global', 1.0)
                
                # Combined weight: ECs (evidence) * SPs (specificity)
                weight = ecs_score * sps_score
                weighted_expr = expr_level * weight
                
                total_score += weighted_expr
                total_weight += weight
            
            # Average weighted score
            if total_weight > 0:
                cell_type_scores[cell_type] = total_score / total_weight
            else:
                cell_type_scores[cell_type] = 0.0
        
        # Find best scoring cell type
        if cell_type_scores:
            best_cell_type = max(cell_type_scores.keys(), key=lambda ct: cell_type_scores[ct])
            best_score = cell_type_scores[best_cell_type]
            
            if best_score > 0:
                print(f"------ Cluster {cluster}: {best_cell_type} (score: {best_score:.3f}, n_cells: {len(cluster_cells)})")
            else:
                best_cell_type = "Unknown"
                print(f"------ Cluster {cluster}: {best_cell_type} (no marker genes found, n_cells: {len(cluster_cells)})")
        else:
            best_cell_type = "Unknown"
            best_score = 0.0
            print(f"------ Cluster {cluster}: {best_cell_type} (no cell types available, n_cells: {len(cluster_cells)})")
        
        # Get ontology information for the best cell type
        cl_id = ontology_mapping.get(best_cell_type, "Unknown")
        
        cluster_results.append({
            'cluster': cluster,
            'cell_type': best_cell_type,
            'CL_celltype': best_cell_type,  # Same as cell_type for consistency
            'CL_ID': cl_id,
            'score': best_score,
            'ncells': len(cluster_cells)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(cluster_results)
    
    # Add cell type annotations to adata
    cluster_to_celltype = dict(zip(results_df['cluster'], results_df['cell_type']))
    adata.obs['cell_type'] = adata.obs[cluster_key].map(cluster_to_celltype)
    
    # Print summary
    print("--- CellMarkerAccordion annotation summary:")
    for _, row in results_df.iterrows():
        print(f"    Cluster {row['cluster']}: {row['cell_type']} (score: {row['score']:.3f}, n_cells: {row['ncells']})")
    
    return adata, results_df