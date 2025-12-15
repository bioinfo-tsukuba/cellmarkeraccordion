"""
CellMarkerAccordion Python Package

A modular package for converting original marker files to light format with Ensembl IDs.
"""

from .biomart_client import BioMartClient, BioMartAPIException
from .ensembl_mapper import EnsemblMapper
from .config_manager import ConfigManager, ConfigurationError
from .marker_processor import MarkerProcessor
from .conversion_pipeline import ConversionPipeline
from .marker_database_integration import ScoreCalculator

__version__ = "1.0.0"
__author__ = "CellMarkerAccordion Team"

__all__ = [
    'BioMartClient',
    'BioMartAPIException',
    'EnsemblMapper', 
    'ConfigManager',
    'ConfigurationError',
    'MarkerProcessor',
    'ConversionPipeline',
    'ScoreCalculator'
]