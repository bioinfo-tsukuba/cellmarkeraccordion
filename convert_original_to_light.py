#!/usr/bin/env python3

import os
import sys
import argparse
import logging

# Add Python package to path for importing
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Python.conversion_pipeline import ConversionPipeline
from Python.config_manager import ConfigurationError
from Python.biomart_client import BioMartAPIException

class BioMartConverter:
    """Legacy wrapper for ConversionPipeline to maintain backwards compatibility."""
    
    def __init__(self, config_file: str = "biomart_config.json", calculate_scores: bool = False):
        """Initialize converter with configuration."""
        self.pipeline = ConversionPipeline(config_file, calculate_scores)
        self.logger = logging.getLogger(__name__)
    
    
    def add_ensembl_id(self, config_id: int = 1) -> None:
        """Main pipeline to add Ensembl IDs to marker files."""
        try:
            self.pipeline.run_conversion(config_id)
        except (ConfigurationError, BioMartAPIException) as e:
            self.logger.error(f"Conversion failed: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error during conversion: {str(e)}")
            raise

def main():
    """Main function to run the converter."""
    parser = argparse.ArgumentParser(
        description='Convert original marker files to light format with Ensembl IDs and calculate ECs/SPs scores (default behavior)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_original_to_light.py
  python convert_original_to_light.py --config-id 1 --no-scores
  python convert_original_to_light.py --config-file custom_config.json --calculate-scores
        """
    )
    parser.add_argument('--config-id', type=int, default=1, help='Configuration ID to use (default: 1)')
    parser.add_argument('--config-file', type=str, default='biomart_config.json', help='Configuration file (default: biomart_config.json)')
    parser.add_argument('--calculate-scores', action='store_true', default=True, help='Calculate ECs and SPs scores (default: True)')
    parser.add_argument('--no-scores', action='store_true', help='Skip score calculation (generate original light format)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose logging')
    parser.add_argument('--species', type=str, help='Process only specific species (optional)')
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = 'DEBUG' if args.verbose else 'INFO'
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    logger = logging.getLogger(__name__)
    
    try:
        # Validate arguments
        if not os.path.exists(args.config_file):
            logger.error(f"Configuration file not found: {args.config_file}")
            return 1
        
        # Determine whether to calculate scores
        calculate_scores = args.calculate_scores and not args.no_scores
        
        # Create converter and run pipeline
        converter = BioMartConverter(config_file=args.config_file, calculate_scores=calculate_scores)
        
        if args.species:
            # Process single species
            logger.debug(f"Processing single species: {args.species}")
            converter.pipeline.process_single_species(args.config_id, args.species)
        else:
            # Process all species
            logger.debug(f"Processing all species in configuration ID: {args.config_id}")
            converter.add_ensembl_id(config_id=args.config_id)
        
        logger.info("Conversion completed!")
        return 0
        
    except ConfigurationError as e:
        logger.error(f"Configuration error: {str(e)}")
        return 1
    except BioMartAPIException as e:
        logger.error(f"BioMart API error: {str(e)}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"File not found: {str(e)}")
        return 1
    except KeyboardInterrupt:
        logger.warning("Conversion interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        if args.verbose:
            import traceback
            logger.debug(traceback.format_exc())
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
    