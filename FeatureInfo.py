#!/usr/bin/env python

# -*- coding: utf-8 -*-


# Module metadata variables
__author__ = "Rafael Barrero Rodriguez"
__credits__ = ["Rafael Barrero Rodriguez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "rbarreror@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"


# Import modules 
import os
import sys
import argparse
import configparser
import logging
from pathlib import Path
import numpy as np
import pandas as pd

import pdb


###################
# Local functions #
###################

def readInfile(infile, row):
    '''
    Input:
        - infile: Path where the file is located
        - row: Index (0-based) of column headers, where the table starts
    Output:
        - df: Pandas dataframe with the table
    '''

    log_str = f'Reading input file: {str(Path(infile))}'
    logging.info(log_str)

    try:
        df = pd.read_excel(infile, header=row)
    
    except:
        log_str = f'Error when reading {str(Path(infile))}'
        logging.info(log_str)
        
        # Log error class and message
        exctype, value = sys.exc_info()[:2]
        log_str = f'{exctype}: {value}'
        logging.info(log_str)

        sys.exit()
    
    log_str = f'{str(Path(infile))} was read'
    logging.info(log_str)

    return df


def mergeTable(feature_table, identification_table):
    '''
    Input:
        - feature_table: Pandas dataframe containing feature ID, Apex m/z and RT
        - identification_table: Pandas dataframe containing identifications
    
    Output:
        - merged_table: Pandas dataframe obtained from merging feature_table and identification_table
    '''

    # Assert that column names are correct
    assert ('Apex m/z' in feature_table.columns) or ('Experimental mass' in feature_table.columns),\
        "Name of the column with feature mass should be 'Apex m/z' or 'Experimental mass'"
    
    assert ('Apex m/z' in identification_table.columns) or ('Experimental mass' in identification_table.columns),\
        "Name of the column with feature mass should be 'Apex m/z' or 'Experimental mass'"
    
    
    # Get number of digits to round numbers
    n_digits = config_param.getint('Parameters', 'N_Digits')

    # Rename table column: Apex m/z --> Experimental mass
    feature_table.columns = [name if name != "Apex m/z" else "Experimental mass" for name in feature_table.columns]

    identification_table.columns = [name if name != "Apex m/z" else "Experimental mass" for name in identification_table.columns]

    # Round numbers in both dataframes
    feature_table = feature_table.round({'Experimental mass': n_digits})
    identification_table = identification_table.round({'Experimental mass': n_digits})

    # Merge both tables
    merged_table = pd.merge(feature_table, identification_table, how = 'outer', on = 'Experimental mass')

    return merged_table


def getOutputFilename():
    '''
    Output:
        - filename: String containing the name of the output file
    '''

    filename = config_param.get('Parameters', 'OutputName')

    if not filename:
        filename = 'id_' + os.path.basename(args.identification)

    elif not os.path.splitext(filename)[1]:
        filename += '.xls'
    
    return filename


def writeDataFrame(df, path):
    '''
    Description: Function used to export pandas dataframe with the tags added

    Input:
        - df: Pandas dataframe that will be exported
        - path: String containing the path to infile. The new file will be saved in the
        same folder.
    '''

    # Build output file path
    output_path = os.path.join(os.path.dirname(path))
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    # Get output file name
    filename = getOutputFilename()

    output_file = os.path.join(output_path, filename)

    # Get output columns
    if config_param.get('Parameters', 'OutputColumns'):
        output_columns = [name.strip() for name in config_param.get('Parameters', 'OutputColumns').split(',') if name.strip() in df.columns]
    
    else:
        output_columns = df.columns

    # Handle errors in exception case
    try:
        df.to_excel(output_file, index=False, columns=output_columns)
    
    except:
        log_str = f'Error when writing {str(Path(output_file))}'
        logging.info(log_str)
        
        # Log error class and message
        exctype, value = sys.exc_info()[:2]
        log_str = f'{exctype}: {value}'
        logging.info(log_str)

        sys.exit()
    
    log_str = f'{str(Path(output_file))} was written'
    logging.info(log_str)



##################
# Main functions #
##################

def main(args):
    '''
    main function
    '''

    # Get table with feature information
    feature_table = readInfile(args.feature, 0)

    # Get table with identifications
    identification_table = readInfile(args.identification, 0)

    # Merge both tables based on experimental mass
    merged_table = mergeTable(feature_table, identification_table)

    # Export merged table
    writeDataFrame(merged_table, args.identification)



if __name__=="__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(
        description='FeatureInfo',
        epilog='''
        Example:
            python FeatureInfo.py
        
        '''
    )

    # Set default values
    default_config_path = os.path.join(os.path.dirname(__file__), "config/configFeatureInfo/configFeatureInfo.ini")

    # Parse arguments corresponding to path files
    parser.add_argument('-if', '--feature', help="Path to input file with feature information", type=str, required=True)
    parser.add_argument('-id', '--identification', help="Path to input file with identifications", type=str, required=True)
    parser.add_argument('-c', '--config', help="Path to configTable.ini file", type=str, default=default_config_path)
 
    # Parse arguments corresponding to parameters
    parser.add_argument('-o', '--output', help="Name of output table", type=str)
    parser.add_argument('-oc', '--outCol', help='Name/Index of columns present in output table. By default, all columns will be displayed', type=str)
    parser.add_argument('-n', '--digits', help="Number of digits", type=str)

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()


    # parse config with user selection
    config_param = configparser.ConfigParser(inline_comment_prefixes='#')
    config_param.read(Path(args.config))

    # Parameters introduced in the execution replace those in the .ini file

    if args.digits:
        config_param.set('Parameters', 'N_Digits', args.digits)

    if args.output:
        config_param.set('Parameters', 'OutputName', args.output)
    
    if args.outCol:
        config_param.set('Parameters', 'OutputColumns', args.outCol)



    # logging debug level. By default, info level
    if args.identification:
        log_file = outfile = os.path.splitext(args.identification)[0] + '_log.txt'
        log_file_debug = outfile = os.path.splitext(args.identification)[0] + '_log_debug.txt'
    
    else:
        log_file = outfile = 'log.txt'
        log_file_debug = outfile = 'log_debug.txt'


    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')