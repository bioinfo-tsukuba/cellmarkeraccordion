#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import sys
from typing import Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import scanpy as sc

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
    # print(f"    Loading CellMarkerAccordion database for process_id: {process_id}")
    
    # Determine database file path based on configuration
    with open("post_processing_config.json", 'r') as f:
        config = pd.read_json(f)
    
    # Find configuration for this process_id
    process_config = None
    for proc in config['post_processing']:
        if proc['process_id'] == str(process_id):
            process_config = proc['configs'][0]  # Take first config
            break
    
    if not process_config:
        raise ValueError(f"No configuration found for process_id: {process_id}")
    
    annotation_config = process_config['Annotation']
    if annotation_config['tool'] != 'cellmarkeraccordion':
        raise ValueError(f"Expected cellmarkeraccordion, got: {annotation_config['tool']}")
    
    # Load biomart_config.json to get version and outfile for this annotation id
    biomart_config_path = f"tools/cellmarkeraccordion/{annotation_config['config']}"
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
    db_path = f"tools/cellmarkeraccordion/{db_version}/{db_file}"
    
    # print(f"    Loading CellMarkerAccordion database from: {db_path}")
    
    try:
        # Load CellMarkerAccordion data
        markers_df = pd.read_csv(db_path)
        # print(f"    Loaded {len(markers_df)} marker entries from CellMarkerAccordion database")
        
        # Filter by tissue if specified
        if tissue_ontology:
            tissue_filtered = markers_df[markers_df['Uberon_ID'] == tissue_ontology]
            if len(tissue_filtered) > 0:
                markers_df = tissue_filtered
                # print(f"    Filtered to {len(markers_df)} markers for tissue: {tissue_ontology}")
            else:
                print(f"    Warning: No markers found for tissue {tissue_ontology}, using all tissues")
        
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
        
        # print(f"    Organized markers for {len(cell_type_markers)} cell types")
        return cell_type_markers, ontology_mapping
        
    except FileNotFoundError:
        # Fallback to v1.0.0 if specific version not found
        db_path_fallback = f"tools/cellmarkeraccordion/v1.0.0/{db_file}"
        print(f"    Fallback: Loading from {db_path_fallback}")
        
        markers_df = pd.read_csv(db_path_fallback)
        print(f"    Loaded {len(markers_df)} marker entries from fallback database")
        
        # Filter by tissue if specified
        if tissue_ontology:
            tissue_filtered = markers_df[markers_df['Uberon_ID'] == tissue_ontology]
            if len(tissue_filtered) > 0:
                markers_df = tissue_filtered
                # print(f"    Filtered to {len(markers_df)} markers for tissue: {tissue_ontology}")
        
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


def perform_cellmarkeraccordion_annotation(adata, cluster_key="louvain", process_id=1, tax_id=9606, tissue_ontology=None):
                                         
    # print(f"    Starting CellMarkerAccordion annotation with process_id: {process_id}")
    
    # Load CellMarkerAccordion database
    cell_type_markers, ontology_mapping = load_cellmarkeraccordion_database(process_id, tax_id, tissue_ontology)
    
    # Get all marker genes from CellMarkerAccordion database
    all_marker_genes = set()
    for cell_type_data in cell_type_markers.values():
        all_marker_genes.update(cell_type_data['positive'])

    print(f"    Found {len(all_marker_genes)} unique marker genes for {tissue_ontology} in CellMarkerAccordion database")
    
    # Use full expression data if provided, otherwise use adata.raw or current adata
    full_expr = adata.to_df().T  # Fallback to processed data
    
    # Set gene names as index - use gene_symbols if available, otherwise use var_names
    if 'gene_symbols' in adata.var.columns:
        full_expr.index = adata.var['gene_symbols']
        # print(f"    Using gene_symbols for gene names")
    else:
        full_expr.index = adata.var_names
        print(f"    Warning: gene_symbols not found, using var_names instead")
    
    # Filter to only marker genes present in dataset (intersect)
    available_marker_genes = [gene for gene in all_marker_genes if gene in full_expr.index]
    expr_data = full_expr.loc[available_marker_genes]
    
    print(f"    Using {len(available_marker_genes)} marker genes present in dataset")
    # print(f"    Expression data shape for CellMarkerAccordion: {expr_data.shape} (marker_genes × cells)")
    
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
                
                # Get ECs and SPs scores for this marker (default to 0 if not found)
                ecs_score = marker_scores.get(marker, {}).get('ECs_global', 0)
                sps_score = marker_scores.get(marker, {}).get('SPs_global', 0)
                
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
        
        # Ensure ALL cell types from database get scores (fill missing with 0.0)
        all_cell_types_scores = {}
        for cell_type in cell_type_markers.keys():
            all_cell_types_scores[cell_type] = cell_type_scores.get(cell_type, 0.0)
        
        # Find best scoring cell type from computed scores only (not 0.0 defaults)
        if cell_type_scores:
            best_cell_type = max(cell_type_scores.keys(), key=lambda ct: cell_type_scores[ct])
            best_score = cell_type_scores[best_cell_type]
        else:
            # If no cell types had any markers, pick first available as default
            best_cell_type = list(all_cell_types_scores.keys())[0] if all_cell_types_scores else "Unknown"
            best_score = 0.0
        
        # Add entries for ALL cell types from database, marking which one is determined (best)
        for cell_type, score in all_cell_types_scores.items():
            cl_id = ontology_mapping.get(cell_type, "Unknown")
            is_determined = (cell_type == best_cell_type)
            
            cluster_results.append({
                'cluster': cluster,
                'cell_type': cell_type,
                'CL_celltype': cell_type,  # Same as cell_type for consistency
                'CL_ID': cl_id,
                'score': score,
                'ncells': len(cluster_cells),
                'is_determined': is_determined
            })
        
        # Handle case where no cell types are available in database
        if not all_cell_types_scores:
            cluster_results.append({
                'cluster': cluster,
                'cell_type': "Unknown",
                'CL_celltype': "Unknown", 
                'CL_ID': "Unknown",
                'score': 0.0,
                'ncells': len(cluster_cells),
                'is_determined': True
            })
            # print(f"       Cluster {cluster}: Unknown (no cell types available, n_cells: {len(cluster_cells)})")
    
    # Create results DataFrame
    results_df = pd.DataFrame(cluster_results)
    
    # Add cell type annotations to adata (only use determined cell types)
    determined_results = results_df[results_df['is_determined'] == True]
    cluster_to_celltype = dict(zip(determined_results['cluster'], determined_results['cell_type']))
    adata.obs['cell_type'] = adata.obs[cluster_key].map(cluster_to_celltype)
    
    # Print summary
    # print("    CellMarkerAccordion annotation summary:")
    # for _, row in results_df.iterrows():
    #     print(f"    Cluster {row['cluster']}: {row['cell_type']} (score: {row['score']:.3f}, n_cells: {row['ncells']})")
    
    return adata, results_df


def create_score_matrix(results_df):
    """
    Create cell type score matrix (cell types × clusters).
    
    Parameters:
        results_df: DataFrame from perform_cellmarkeraccordion_annotation
        
    Returns:
        pd.DataFrame: Score matrix with cell types as rows, clusters as columns
    """
    # Pivot the results to get score matrix
    score_matrix = results_df.pivot(index='cell_type', columns='cluster', values='score')
    score_matrix = score_matrix.fillna(0)  # Fill any missing values with 0
    
    # Sort clusters by name/number for consistent ordering
    score_matrix = score_matrix.reindex(sorted(score_matrix.columns), axis=1)
    
    return score_matrix


def create_determined_score_matrix(results_df):
    """
    Create filtered score matrix showing only determined cell types for each cluster.
    
    Parameters:
        results_df: DataFrame from perform_cellmarkeraccordion_annotation
        
    Returns:
        pd.DataFrame: Filtered score matrix with determined cell types as rows, clusters as columns
    """
    # Get only determined results
    determined_results = results_df[results_df['is_determined'] == True].copy()
    
    if len(determined_results) == 0:
        print("No determined cell types found")
        return None
    
    # Create cluster to determined cell type mapping
    cluster_to_celltype = dict(zip(determined_results['cluster'], determined_results['cell_type']))
    
    # Get all relevant scores for the determined cell types
    determined_celltypes = set(determined_results['cell_type'])
    filtered_results = results_df[
        (results_df['cell_type'].isin(determined_celltypes)) &
        (results_df['cluster'].isin(cluster_to_celltype.keys()))
    ].copy()
    
    # Create score matrix
    score_matrix = filtered_results.pivot(index='cell_type', columns='cluster', values='score')
    score_matrix = score_matrix.fillna(0)
    
    # Sort clusters by name/number for consistent ordering
    score_matrix = score_matrix.reindex(sorted(score_matrix.columns), axis=1)
    
    # Sort cell types to put determined ones first
    determined_order = []
    other_celltypes = []
    
    for cluster in sorted(score_matrix.columns):
        determined_celltype = cluster_to_celltype.get(cluster)
        if determined_celltype and determined_celltype not in determined_order:
            determined_order.append(determined_celltype)
    
    # Add any remaining cell types
    for celltype in score_matrix.index:
        if celltype not in determined_order:
            other_celltypes.append(celltype)
    
    final_order = determined_order + sorted(other_celltypes)
    score_matrix = score_matrix.reindex(final_order, axis=0)
    
    return score_matrix


def create_ordered_heatmap(score_matrix, results_df, output_path):
    """
    Create heatmap ordered by determined cell types.
    
    Parameters:
        score_matrix: Score matrix from create_score_matrix
        results_df: Original results DataFrame
        output_path: Path to save the heatmap PNG
    """
    # Get determined cell types for each cluster
    determined_results = results_df[results_df['is_determined'] == True]
    cluster_to_celltype = dict(zip(determined_results['cluster'], determined_results['cell_type']))
    
    # Order clusters by their determined cell type
    ordered_clusters = []
    celltype_order = []
    
    for cluster in sorted(score_matrix.columns):
        determined_celltype = cluster_to_celltype.get(cluster, 'Unknown')
        ordered_clusters.append(cluster)
        celltype_order.append(determined_celltype)
    
    # Reorder score matrix columns
    ordered_score_matrix = score_matrix[ordered_clusters]
    
    # Create cluster labels with determined cell type
    cluster_labels = [f"Cluster {cluster}\n({cluster_to_celltype.get(cluster, 'Unknown')})" 
                     for cluster in ordered_clusters]
    
    # Create heatmap
    plt.figure(figsize=(max(12, len(ordered_clusters) * 0.8), max(8, len(score_matrix.index) * 0.4)))
    
    sns.heatmap(ordered_score_matrix, 
                annot=True, 
                fmt='.2f', 
                cmap='viridis',
                cbar_kws={'label': 'CellMarkerAccordion Score'},
                xticklabels=cluster_labels,
                yticklabels=ordered_score_matrix.index)
    
    plt.title('CellMarkerAccordion Scores Heatmap\nOrdered by Determined Cell Types')
    plt.xlabel('Clusters (Determined Cell Type)')
    plt.ylabel('Cell Types')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def calculate_marker_correlations(adata, results_df, cell_type_markers, cluster_key="louvain"):
    """
    Calculate correlation between marker gene expression and cell type scores.
    
    Parameters:
        adata: AnnData object with expression data
        results_df: Results from cellmarkeraccordion annotation
        cell_type_markers: Dictionary of cell type markers from database
        cluster_key: Key for cluster annotations
        
    Returns:
        pd.DataFrame: Best marker genes for each cluster-celltype assignment
    """
    # Get expression data
    if 'gene_symbols' in adata.var.columns:
        expr_data = adata.to_df().T
        expr_data.index = adata.var['gene_symbols']
    else:
        expr_data = adata.to_df().T
        expr_data.index = adata.var_names
    
    # Get determined cell types
    determined_results = results_df[results_df['is_determined'] == True]
    
    # Calculate cluster mean expression for all genes
    cluster_expr = {}
    for cluster in determined_results['cluster'].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_cells = adata.obs.index[cluster_mask]
        cluster_expr[cluster] = expr_data[cluster_cells].mean(axis=1)
    
    # Convert to DataFrame for easier handling
    cluster_expr_df = pd.DataFrame(cluster_expr)
    
    # Calculate correlations for each determined cluster-celltype pair
    best_markers = []
    
    for _, row in determined_results.iterrows():
        cluster = row['cluster']
        cell_type = row['cell_type']
        
        if cell_type not in cell_type_markers:
            continue
            
        positive_markers = cell_type_markers[cell_type]['positive']
        
        # Get cell type scores across all clusters for this cell type
        celltype_scores = results_df[results_df['cell_type'] == cell_type]
        score_by_cluster = dict(zip(celltype_scores['cluster'], celltype_scores['score']))
        
        # Calculate correlations for each marker gene
        marker_correlations = []
        
        for marker in positive_markers:
            if marker not in cluster_expr_df.index:
                continue
                
            # Get marker expression across clusters
            marker_expr_by_cluster = []
            scores_by_cluster = []
            
            for cluster_id in cluster_expr_df.columns:
                if cluster_id in score_by_cluster:
                    marker_expr_by_cluster.append(cluster_expr_df.loc[marker, cluster_id])
                    scores_by_cluster.append(score_by_cluster[cluster_id])
            
            if len(marker_expr_by_cluster) >= 3:  # Need at least 3 points for correlation
                try:
                    corr_coef, p_value = pearsonr(marker_expr_by_cluster, scores_by_cluster)
                    if not np.isnan(corr_coef):
                        marker_correlations.append({
                            'marker': marker,
                            'correlation': abs(corr_coef),  # Use absolute correlation
                            'correlation_coef': corr_coef,
                            'p_value': p_value
                        })
                except:
                    continue
        
        # Find best marker (highest absolute correlation)
        if marker_correlations:
            best_marker = max(marker_correlations, key=lambda x: x['correlation'])
            best_markers.append({
                'cluster': cluster,
                'cell_type': cell_type,
                'marker_gene': best_marker['marker'],
                'correlation': best_marker['correlation_coef'],
                'p_value': best_marker['p_value']
            })
    
    return pd.DataFrame(best_markers)


def create_violin_plots(adata, best_markers_df, cluster_key="louvain", output_path=None):
    """
    Create violin plots for best marker genes.
    
    Parameters:
        adata: AnnData object
        best_markers_df: DataFrame with best marker genes from calculate_marker_correlations
        cluster_key: Key for cluster annotations
        output_path: Path to save the violin plot PNG
    """
    # Get expression data
    if 'gene_symbols' in adata.var.columns:
        expr_data = adata.to_df()
        gene_to_var = dict(zip(adata.var['gene_symbols'], adata.var_names))
    else:
        expr_data = adata.to_df()
        gene_to_var = dict(zip(adata.var_names, adata.var_names))
    
    # Prepare data for plotting
    plot_data = []
    
    for _, row in best_markers_df.iterrows():
        cluster = row['cluster']
        cell_type = row['cell_type']
        marker_gene = row['marker_gene']
        
        if marker_gene not in gene_to_var:
            continue
            
        var_name = gene_to_var[marker_gene]
        
        # Get cells in this cluster
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_cells = adata.obs.index[cluster_mask]
        
        # Get expression values
        for cell in cluster_cells:
            plot_data.append({
                'cluster': f"Cluster {cluster}",
                'cell_type': cell_type,
                'marker_gene': marker_gene,
                'expression_level': expr_data.loc[cell, var_name],
                'cluster_celltype': f"Cluster {cluster}\n({cell_type})"
            })
    
    if not plot_data:
        print("No data available for violin plots")
        return
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create violin plots faceted by cell type
    unique_celltypes = plot_df['cell_type'].unique()
    n_celltypes = len(unique_celltypes)
    
    fig, axes = plt.subplots(1, n_celltypes, figsize=(4 * n_celltypes, 6), sharey=False)
    
    if n_celltypes == 1:
        axes = [axes]
    
    for i, cell_type in enumerate(unique_celltypes):
        celltype_data = plot_df[plot_df['cell_type'] == cell_type]
        
        # Order clusters for this cell type
        cluster_order = sorted(celltype_data['cluster_celltype'].unique())
        
        sns.violinplot(data=celltype_data, 
                      x='cluster_celltype', 
                      y='expression_level',
                      ax=axes[i])
        
        axes[i].set_title(f'{cell_type}')
        axes[i].set_xlabel('Cluster')
        axes[i].set_ylabel('Expression Level')
        axes[i].tick_params(axis='x', rotation=45)
        
        # Add marker gene name to subplot
        marker_gene = celltype_data['marker_gene'].iloc[0]
        axes[i].text(0.02, 0.98, f'Gene: {marker_gene}', 
                    transform=axes[i].transAxes, 
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle('Expression of Best Marker Genes by Cluster')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def create_validation_violin_plot_scanpy_style(adata, best_markers_df, cluster_key="leiden", output_path="validation_violin.png"):
    """
    Create validation violin plot using scanpy style matching create_sanity_check_violin_plot.
    
    Parameters:
        adata: AnnData object
        best_markers_df: DataFrame with columns ['cluster', 'cell_type', 'marker_gene', 'correlation']
        cluster_key: Key for cluster annotations  
        output_path: Path to save the violin plot PNG
    """
    if best_markers_df is None or len(best_markers_df) == 0:
        print("No best marker genes available for violin plot")
        return
    
    # Deduplicate marker genes - keep the one with highest absolute correlation
    unique_markers = {}
    for _, row in best_markers_df.iterrows():
        marker_gene = row['marker_gene']
        correlation = abs(row.get('correlation', 0))  # Use absolute correlation
        
        if marker_gene not in unique_markers or correlation > unique_markers[marker_gene]['correlation']:
            unique_markers[marker_gene] = {
                'cluster': row['cluster'],
                'cell_type': row['cell_type'],
                'marker_gene': marker_gene,
                'correlation': correlation
            }
    
    # Convert to DataFrame for further processing
    unique_markers_df = pd.DataFrame(list(unique_markers.values()))
    
    # Convert unique best_markers to tissue_markers format: {celltype: {"marker": gene}}
    tissue_markers = {}
    for _, row in unique_markers_df.iterrows():
        cell_type = row['cell_type']
        marker_gene = row['marker_gene']
        
        # Use first occurrence if multiple cell types have same marker
        if cell_type not in tissue_markers:
            tissue_markers[cell_type] = {"marker": marker_gene}
    
    if not tissue_markers:
        print("No tissue markers generated from best marker genes")
        return
    
    # Create gene to celltype mapping (same as original function)
    gene2celltype = {
        v["marker"]: k
        for k, v in tissue_markers.items()
    }
    
    # Get cluster to celltype mapping from adata.obs
    cluster_celltype_gene = (
        adata.obs[[cluster_key, "cell_type"]]
        .dropna()
        .groupby(cluster_key)["cell_type"]
        .agg(lambda x: x.value_counts().idxmax())
        .sort_index()
    )
    cluster_celltype_df = cluster_celltype_gene.rename("cell_type").to_frame()
    cluster_celltype_df["gene"] = cluster_celltype_df["cell_type"].map(
        {k: v["marker"] for k, v in tissue_markers.items()}
    )
    
    # Order genes based on cluster assignments  
    all_genes = list(gene2celltype.keys())
    ordered_genes = cluster_celltype_df["gene"].dropna().tolist()
    rest = sorted([g for g in all_genes if g not in set(ordered_genes)])
    ordered_genes.extend(rest)
    
    if not ordered_genes:
        print("No genes available for violin plot")
        return
    
    # Filter ordered_genes to only include genes present in adata
    if 'gene_symbols' in adata.var.columns:
        available_genes = adata.var['gene_symbols'].tolist()
    else:
        available_genes = adata.var_names.tolist()
    
    ordered_genes = [g for g in ordered_genes if g in available_genes]
    
    if not ordered_genes:
        print("None of the marker genes found in dataset")
        return
    
    # Create stacked violin plot (same style as original)
    try:
        ax = sc.pl.stacked_violin(
            adata,
            ordered_genes,
            groupby=cluster_key,
            cmap="viridis_r",
            swap_axes=True,
            show=False,
        )
        
        # Get axis objects
        main_ax = ax["mainplot_ax"]
        cbar_ax = ax["color_legend_ax"]
        
        # Position colorbar (same as original)
        cbar_ax.set_position([0.88, 0.88, 0.2, 0.025])
        fig = main_ax.figure
        
        # Create right twin axis for cell type labels
        ax_r = main_ax.twinx()
        ax_r.set_ylim(main_ax.get_ylim())
        ax_r.set_yticks(main_ax.get_yticks())
        gene_celltypes = [gene2celltype.get(g, "") for g in ordered_genes]
        ax_r.set_yticklabels(gene_celltypes)
        ax_r.tick_params(axis="y", length=0)  # hide tick marks
        ax_r.set_ylabel("cell type")
        
        # Create top twin axis for cluster cell type labels
        ax_t = main_ax.twiny()
        ax_t.set_xlim(main_ax.get_xlim())
        xticks = main_ax.get_xticks()
        ax_t.set_xticks(xticks)
        cluster_labels = [t.get_text() for t in main_ax.get_xticklabels()]
        cluster2ct = cluster_celltype_df["cell_type"].astype(str)
        cluster2ct.index = cluster2ct.index.astype(str)
        top_labels = [cluster2ct.get(c, "") for c in cluster_labels]
        ax_t.set_xticklabels(top_labels, rotation=90)
        ax_t.tick_params(axis="x", length=0)          # remove tick marks
        ax_t.spines["top"].set_visible(False)         # optional
        ax_t.set_xlabel("cluster cell type")
        
        # Save figure (same parameters as original)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        
    except Exception as e:
        print(f"Error creating scanpy violin plot: {str(e)}")
        # Fallback to basic matplotlib if scanpy fails
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, f'Violin plot generation failed:\n{str(e)}', 
                ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Validation Violin Plot - Error')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


def create_complete_violin_plot_scanpy_style(adata, cell_type_markers, cluster_key="leiden", output_path="complete_violin.png"):
    """
    Create complete violin plot showing ALL marker genes from ALL cell types in database using scanpy style.
    
    Parameters:
        adata: AnnData object
        cell_type_markers: Dictionary of ALL cell type markers from database
        cluster_key: Key for cluster annotations  
        output_path: Path to save the violin plot PNG
    """
    if not cell_type_markers:
        print("No cell type markers available for complete violin plot")
        return
    
    # Convert ALL cell_type_markers to tissue_markers format 
    # Use original CDP_fun format: {celltype: {"marker": gene}}
    tissue_markers = {}
    for cell_type, markers_info in cell_type_markers.items():
        positive_markers = markers_info.get('positive', [])
        if positive_markers:
            # For complete view, use first positive marker for each cell type
            tissue_markers[cell_type] = {"marker": positive_markers[0]}
    
    if not tissue_markers:
        print("No tissue markers available for complete violin plot")
        return
    
    # Create gene to celltype mapping (clean cell type names without cluster info)
    gene2celltype = {
        v["marker"]: k  # k is clean cell type name
        for k, v in tissue_markers.items()
    }
    
    # Get cluster to celltype mapping from adata.obs
    cluster_celltype_gene = (
        adata.obs[[cluster_key, "cell_type"]]
        .dropna()
        .groupby(cluster_key)["cell_type"]
        .agg(lambda x: x.value_counts().idxmax())
        .sort_index()
    )
    cluster_celltype_df = cluster_celltype_gene.rename("cell_type").to_frame()
    cluster_celltype_df["gene"] = cluster_celltype_df["cell_type"].map(
        {k: v["marker"] for k, v in tissue_markers.items()}
    )
    
    # Order genes - prioritize determined cluster assignments first
    all_genes = list(gene2celltype.keys())
    ordered_genes = cluster_celltype_df["gene"].dropna().tolist()
    
    # Add remaining genes not in cluster order
    rest = sorted([g for g in all_genes if g not in set(ordered_genes)])
    ordered_genes.extend(rest)
    
    if not ordered_genes:
        print("No genes available for complete violin plot")
        return
    
    # Filter ordered_genes to only include genes present in adata
    if 'gene_symbols' in adata.var.columns:
        available_genes = adata.var['gene_symbols'].tolist()
    else:
        available_genes = adata.var_names.tolist()
    
    ordered_genes = [g for g in ordered_genes if g in available_genes]
    
    if not ordered_genes:
        print("None of the marker genes found in dataset for complete violin plot")
        return
    
    # Create stacked violin plot (same style as original)
    try:
        ax = sc.pl.stacked_violin(
            adata,
            ordered_genes,
            groupby=cluster_key,
            cmap="viridis_r",
            swap_axes=True,
            show=False,
        )
        
        # Get axis objects
        main_ax = ax["mainplot_ax"]
        cbar_ax = ax["color_legend_ax"]
        
        # Position colorbar at top-right (original good position)
        cbar_ax.set_position([0.88, 0.88, 0.2, 0.025])
        fig = main_ax.figure
        
        # Create right twin axis for cell type labels (clean names only)
        ax_r = main_ax.twinx()
        ax_r.set_ylim(main_ax.get_ylim())
        ax_r.set_yticks(main_ax.get_yticks())
        gene_celltypes = [gene2celltype.get(g, "") for g in ordered_genes]
        ax_r.set_yticklabels(gene_celltypes)
        ax_r.tick_params(axis="y", length=0)  # hide tick marks
        ax_r.set_ylabel("cell type")
        
        # Create top twin axis for cluster cell type labels
        ax_t = main_ax.twiny()
        ax_t.set_xlim(main_ax.get_xlim())
        xticks = main_ax.get_xticks()
        ax_t.set_xticks(xticks)
        cluster_labels = [t.get_text() for t in main_ax.get_xticklabels()]
        cluster2ct = cluster_celltype_df["cell_type"].astype(str)
        cluster2ct.index = cluster2ct.index.astype(str)
        top_labels = [cluster2ct.get(c, "") for c in cluster_labels]
        ax_t.set_xticklabels(top_labels, rotation=90)
        ax_t.tick_params(axis="x", length=0)          # remove tick marks
        ax_t.spines["top"].set_visible(False)         # optional
        ax_t.set_xlabel("cluster cell type")
        
        # Save figure (same parameters as original)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        
    except Exception as e:
        print(f"Error creating complete scanpy violin plot: {str(e)}")
        # Fallback to basic matplotlib if scanpy fails
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, f'Complete violin plot generation failed:\n{str(e)}', 
                ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Complete Violin Plot - Error')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


def create_combined_heatmap_violin_plot(adata, results_df, best_markers_df, cluster_key="leiden", output_path="combined_validation.png"):
    """
    Create combined side-by-side visualization with filtered heatmap and unique marker violin plot.
    
    Parameters:
        adata: AnnData object
        results_df: Results DataFrame from cellmarkeraccordion annotation
        best_markers_df: Best marker genes DataFrame
        cluster_key: Key for cluster annotations
        output_path: Path to save the combined plot
    """
    # Create filtered score matrix (determined cell types only)
    determined_score_matrix = create_determined_score_matrix(results_df)
    
    if determined_score_matrix is None:
        print("Cannot create combined plot: no determined cell types found")
        return
    
    # Deduplicate marker genes for violin plot - prioritize unique markers per cell type
    if best_markers_df is None or len(best_markers_df) == 0:
        print("Cannot create combined plot: no best marker genes available")
        return
    
    # Get unique markers with highest correlation, ensuring one marker per cell type
    unique_markers_per_celltype = {}
    unique_markers_per_gene = {}
    
    # Sort by absolute correlation to prioritize highest correlations
    sorted_markers = best_markers_df.copy()
    sorted_markers['abs_correlation'] = sorted_markers['correlation'].abs()
    sorted_markers = sorted_markers.sort_values('abs_correlation', ascending=False)
    
    for _, row in sorted_markers.iterrows():
        marker_gene = row['marker_gene']
        cell_type = row['cell_type']
        correlation = row.get('correlation', 0)
        
        # Priority 1: One unique marker per cell type (no sharing)
        # Priority 2: Highest absolute correlation
        if (cell_type not in unique_markers_per_celltype and 
            marker_gene not in unique_markers_per_gene):
            
            unique_markers_per_celltype[cell_type] = {
                'cluster': row['cluster'],
                'cell_type': cell_type,
                'marker_gene': marker_gene,
                'correlation': correlation
            }
            unique_markers_per_gene[marker_gene] = cell_type
    
    unique_markers_df = pd.DataFrame(list(unique_markers_per_celltype.values()))
    
    # Convert to tissue_markers format for scanpy violin
    tissue_markers = {}
    for _, row in unique_markers_df.iterrows():
        cell_type = row['cell_type']
        marker_gene = row['marker_gene']
        if cell_type not in tissue_markers:
            tissue_markers[cell_type] = {"marker": marker_gene}
    
    # Get cluster to celltype mapping
    determined_results = results_df[results_df['is_determined'] == True]
    cluster_to_celltype = dict(zip(determined_results['cluster'], determined_results['cell_type']))
    
    # Prepare gene ordering for violin plot
    gene2celltype = {v["marker"]: k for k, v in tissue_markers.items()}
    cluster_celltype_gene = (
        adata.obs[[cluster_key, "cell_type"]]
        .dropna()
        .groupby(cluster_key)["cell_type"]
        .agg(lambda x: x.value_counts().idxmax())
        .sort_index()
    )
    cluster_celltype_df = cluster_celltype_gene.rename("cell_type").to_frame()
    cluster_celltype_df["gene"] = cluster_celltype_df["cell_type"].map(
        {k: v["marker"] for k, v in tissue_markers.items()}
    )
    
    # Order genes
    all_genes = list(gene2celltype.keys())
    ordered_genes = cluster_celltype_df["gene"].dropna().tolist()
    rest = sorted([g for g in all_genes if g not in set(ordered_genes)])
    ordered_genes.extend(rest)
    
    # Filter genes available in dataset
    if 'gene_symbols' in adata.var.columns:
        available_genes = adata.var['gene_symbols'].tolist()
    else:
        available_genes = adata.var_names.tolist()
    ordered_genes = [g for g in ordered_genes if g in available_genes]
    ordered_genes = list(dict.fromkeys(ordered_genes))

    if not ordered_genes:
        print("Cannot create combined plot: no marker genes found in dataset")
        return
    
    # Calculate dynamic figure dimensions based on data components
    n_clusters = len(determined_score_matrix.columns)
    n_celltypes = len(determined_score_matrix.index) 
    n_genes = len(ordered_genes)
    
    # Dynamic sizing factors (reduced for smaller plots relative to text)
    cluster_width_factor = 1.0  # width per cluster (was 1.5)
    celltype_height_factor = 0.25  # height per cell type (was 0.8)
    gene_height_factor = 0.2     # height per gene (was 0.4)
    base_width = 6               # base width for layout (was 8)
    base_height = 3              # base height for layout (was 4)
    
    # Calculate figure size (reduced minimums for more compact plots)
    fig_width = max(8, n_clusters * cluster_width_factor + base_width)
    fig_height = max(4, (n_celltypes * celltype_height_factor) + (n_genes * gene_height_factor) + base_height)
    
    # Create figure with dynamic sizing
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        height_ratios=[4, 9],      # ← plot_spacer() の高さ比（好みで調整）
        width_ratios=[4, 6],       # ← heatmap : violin の横幅比
        hspace=0.0, wspace=0.15
    )
    ax_spacer = fig.add_subplot(gs[0, 0])
    ax_spacer.axis("off")
    
    try:
        # Left subplot: Filtered heatmap (smaller)
        ax_heatmap = fig.add_subplot(gs[1, 0])
        
        # Create heatmap of determined cell types
        sns.heatmap(determined_score_matrix, 
                   annot=False, 
                   cmap='viridis',
                   cbar_kws={'label': ' '},
                   ax=ax_heatmap)
        
        # Create simple cluster labels for bottom axis
        cluster_labels = [cluster for cluster in determined_score_matrix.columns]
        ax_heatmap.set_xticklabels(cluster_labels, rotation=0)
        ax_heatmap.set_yticklabels(determined_score_matrix.index, rotation=0)
        ax_heatmap.set_xlabel('Clusters')
        ax_heatmap.set_ylabel('Cell Types')
        ax_heatmap.set_title('CellMarkerAccordion Scores')
        
        # Add top x-axis with cell type labels rotated 90 degrees
        ax_heatmap_top = ax_heatmap.twiny()
        ax_heatmap_top.set_xlim(ax_heatmap.get_xlim())
        ax_heatmap_top.set_xticks(ax_heatmap.get_xticks())
        celltype_labels = [cluster_to_celltype.get(cluster, 'Unknown') for cluster in determined_score_matrix.columns]
        ax_heatmap_top.set_xticklabels(celltype_labels, rotation=90)
        ax_heatmap_top.set_xlabel('Cell Types')
        ax_heatmap_top.tick_params(axis="x", length=0)
        
        # Right subplot: Violin plot using scanpy (larger)
        ax_violin = fig.add_subplot(gs[:, 1])
        
        # Create temporary scanpy plot
        violin_axes = sc.pl.stacked_violin(
            adata,
            ordered_genes,
            groupby=cluster_key,
            cmap="viridis_r",
            swap_axes=True,
            show=False,
            ax=ax_violin
        )
        
        # Extract axes exactly like original CDP_fun function
        main_ax = violin_axes["mainplot_ax"]
        cbar_ax = violin_axes["color_legend_ax"]
        
        # Use same direct manipulation strategy as original function
        fig = main_ax.figure
        
        # Add twin axes for labels (clean cell type names without cluster info)
        ax_r = main_ax.twinx()
        ax_r.set_ylim(main_ax.get_ylim())
        ax_r.set_yticks(main_ax.get_yticks())
        gene_celltypes = [gene2celltype.get(g, "") for g in ordered_genes]
        ax_r.set_yticklabels(gene_celltypes)
        ax_r.tick_params(axis="y", length=0)  # hide tick marks
        ax_r.set_ylabel("Cell type for marker gene")
        
        # Style violin plot x-axis to match heatmap
        cluster_labels_violin = [t.get_text() for t in main_ax.get_xticklabels()]
        
        # Create simple cluster labels for bottom axis (matching heatmap)
        simple_cluster_labels = [cluster_label for cluster_label in cluster_labels_violin]
        main_ax.set_xticklabels(simple_cluster_labels, rotation=0)
        main_ax.set_xlabel('Clusters')
        
        # Add top twin axis with cell type labels rotated 90 degrees
        ax_t = main_ax.twiny()
        ax_t.set_xlim(main_ax.get_xlim())
        xticks = main_ax.get_xticks()
        ax_t.set_xticks(xticks)
        cluster2ct = cluster_celltype_df["cell_type"].astype(str)
        cluster2ct.index = cluster2ct.index.astype(str)
        top_celltype_labels = [cluster2ct.get(c, 'Unknown') for c in cluster_labels_violin]
        ax_t.set_xticklabels(top_celltype_labels, rotation=90)
        ax_t.set_xlabel('Cell Types')
        ax_t.tick_params(axis="x", length=0)
        ax_t.spines["top"].set_visible(True)
        
        main_ax.set_title('Marker Gene Expression')
        
        # Adjust layout and save
        plt.tight_layout()
        cbar_ax.set_position([0.85, 0.75, 0.1, 0.025])
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error creating combined plot: {str(e)}")
        # Fallback to simple layout
        plt.close()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Simple heatmap
        sns.heatmap(determined_score_matrix, annot=True, fmt='.2f', 
                   cmap='viridis', ax=ax1)
        ax1.set_title('Determined Cell Type Scores')
        
        # Simple text plot for violin
        ax2.text(0.5, 0.5, f'Violin plot generation failed:\n{str(e)}', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Best Marker Genes')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


def validate_cellmarkeraccordion(adata, results_df, cell_type_markers, 
                                cluster_key="louvain", output_dir="validation_output",
                                srx_id="unknown"):
    """
    Main validation function that creates all 4 required outputs.
    
    Parameters:
        adata: AnnData object with expression data
        results_df: Results DataFrame from perform_cellmarkeraccordion_annotation
        cell_type_markers: Cell type markers dictionary from database
        cluster_key: Key for cluster annotations
        output_dir: Directory to save output files
        srx_id: SRX ID for file naming
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Creating CellMarkerAccordion validation outputs for {srx_id}...")
    
    # 1. Create score matrix and save as CSV
    print("  1. Creating score matrix...")
    score_matrix = create_score_matrix(results_df)
    score_matrix_path = os.path.join(output_dir, f"{srx_id}_score_matrix.csv")
    score_matrix.to_csv(score_matrix_path)
    print(f"     Saved score matrix to: {score_matrix_path}")
    
    # 2. Create ordered heatmap
    print("  2. Creating ordered heatmap...")
    heatmap_path = os.path.join(output_dir, f"{srx_id}_score_heatmap.png")
    create_ordered_heatmap(score_matrix, results_df, heatmap_path)
    print(f"     Saved heatmap to: {heatmap_path}")
    
    # 3. Calculate marker correlations and create gene list
    print("  3. Calculating marker gene correlations...")
    best_markers_df = calculate_marker_correlations(adata, results_df, cell_type_markers, cluster_key)
    gene_list_path = os.path.join(output_dir, f"{srx_id}_best_marker_genes.csv")
    best_markers_df.to_csv(gene_list_path, index=False)
    print(f"     Saved best marker genes to: {gene_list_path}")
    
    # 4. Create violin plots
    print("  4. Creating violin plots...")
    violin_path = os.path.join(output_dir, f"{srx_id}_violin_plots.png")
    create_violin_plots(adata, best_markers_df, cluster_key, violin_path)
    print(f"     Saved violin plots to: {violin_path}")
    
    print(f"CellMarkerAccordion validation completed for {srx_id}!")
    
    return {
        'score_matrix': score_matrix,
        'best_markers': best_markers_df,
        'output_files': {
            'score_matrix': score_matrix_path,
            'heatmap': heatmap_path,
            'gene_list': gene_list_path,
            'violin_plots': violin_path
        }
    }