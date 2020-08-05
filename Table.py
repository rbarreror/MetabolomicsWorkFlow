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
import re

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


def extractColumnNames(user_string, all_columns):
    '''
    Input:
        - user_string: String with column names/index (base 1) separated by ','
        - all_columns: Pandas index object, containing the name of the columns

    Output:
        - String Numpy array containing the name of the columns selected by the user
    '''

    user_list = [element.strip() for element in user_string.split(',')]

    # If there is a column with a non-numeric character...
    if any([True if re.search('\D', element) else False for element in user_list]):
        
        columns_out = [element for element in user_list if element in all_columns]
    
    # If all columns are numbers...
    else:

        columns_out = [all_columns[int(element)-1] for element in user_list]
    
    return np.array(columns_out, dtype=str)


def sortIndexes(serie_values):
    '''
    Input:
        - serie_values: Pandas Serie object containing values of the first comparing column and its associated index
    
    Output:
        - index_sorted: List of integers corresponding to the indexes of the values in serie_values, sorted according
        to those values
    '''

    # Get list of tuples: (index, value)
    index_value = [element for element in zip(serie_values.index, serie_values)]
    
    # Sort list of tuples by value
    index_value_sorted = sorted(index_value, key=lambda x: x[1])
    
    # Get indexes
    index_sorted = [element[0] for element in index_value_sorted]

    return index_sorted


def getConserveValues(conserve_values_i, conserve_values_j):
    '''
    Input:
        - conserve_values_i: String Numpy array with the values to be conserved in the upper row
        - conserve_values_j: String Numpy array with the values to be conserved in the downer row
    
    Output:
        - final_conserve: List of strings with the updated upper row values in conserved columns
    '''

    final_conserve = [element[0] if element[1] in element[0].split(' // ') else ' // '.join(element) \
        for element in zip(conserve_values_i, conserve_values_j)]

    return final_conserve


def fusionFunction(df, compare_column, conserve_column, sorted_indexes):
    '''
    Input:
        - df: Pandas dataframe corrresponding to the input table
        - compare_column: String Numpy array with the name of the columns to be compared
        - conserve_column: String Numpy array with the name of the columns to be conserved
        - sorted_indexes: List of integers corresponding to the sorted indexes of the dataframe. The dataframe will
        be iterated following these indexes.
    
    Output:
        - df: Pandas dataframe processed.
    '''

    logging.info("Start collapsing table")

    # List to store dropped indexes
    dropped_indexes = []
    
    # Lambda function used to get values from df given the index and the column names
    getValues = lambda index, columns: np.array([df.at[index, column] if not pd.isna(df.at[index, column]) else "" for column in columns])

    for i, index_i in enumerate(sorted_indexes):

        if index_i in dropped_indexes:
            continue
        
        compare_values_i = getValues(index_i, compare_column)

        for j, index_j in enumerate(sorted_indexes):

            if (j <= i) or index_j in dropped_indexes:
                continue

            compare_values_j = getValues(index_j, compare_column)

            # If values used to sort indexes are different, leave for loop
            if compare_values_i[0] != compare_values_j[0]:
                break

            # Else if some other values are different, move to next iteration
            elif any(compare_values_i != compare_values_j):
                continue

            # Else, proceed with the fusion
            else:
                conserve_values_i = getValues(index_i, conserve_column)
                conserve_values_j = getValues(index_j, conserve_column)

                final_conserve_values = getConserveValues(conserve_values_i, conserve_values_j)
                
                # Add to conserve_values_i to update for next iterations inside the second loop
                conserve_values_i = final_conserve_values

                # Update the upper row
                df.loc[index_i, conserve_column] = final_conserve_values

                # Drop j row
                df.drop(labels=index_j, axis=0, inplace=True)

                # Add dropped index to array
                dropped_indexes.append(index_j)
    
    logging.info("End collapsing table")
    
    return df


def collapsingTable(df):
    '''
    Input:
        - df: Pandas dataframe corresponding to input table
    
    Output:
        - df: Pandas dataframe with collapsed rows following user parameters
    '''

    # Get numpy array with columns used to compare
    compare_column = extractColumnNames(config_param.get('Parameters', 'ComparedColumns'), df.columns)

    # Get numpy array with "conserved" columns
    conserve_column = extractColumnNames(config_param.get('Parameters', 'ConservedColumns'), df.columns)

    # Get indexes sorted by the first comparing column values
    sorted_indexes = sortIndexes(df.loc[:, compare_column[0]])

    # Make the fusion
    df = fusionFunction(df, compare_column, conserve_column, sorted_indexes)

    return df


def getOutputFilename():
    '''
    Output:
        - filename: String containing the name of the output file
    '''

    filename = config_param.get('Parameters', 'OutputName')

    if not filename:
        filename = 'table_' + os.path.basename(args.infile)

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
    # output_path = os.path.join(os.path.dirname(path), "results")
    output_path = os.path.join(os.path.dirname(path))
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    # Get output file name
    filename = getOutputFilename()

    output_file = os.path.join(output_path, filename)

    # Get output columns
    if config_param.get('Parameters', 'OutputColumns'):
        output_columns = extractColumnNames(config_param.get('Parameters', 'OutputColumns'), df.columns)
    
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

    # Read input file
    df = readInfile(args.infile, 0)

    # Proceed with row fusion
    df = collapsingTable(df)

    # Write output dataframe
    writeDataFrame(df, args.infile)


if __name__=="__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Table',
        epilog='''
        Example:
            python Table.py
        
        '''
    )

    # Set default values
    default_config_path = os.path.join(os.path.dirname(__file__), "config/configTable/configTable.ini")
    default_compare_columns = "1, 3, 4, 5, 7, 8, 9"
    default_conserve_columns = "2, 6"
 
    # Parse arguments corresponding to path files
    parser.add_argument('-i', '--infile', help="Path to input file", type=str, required=True)
    parser.add_argument('-c', '--config', help="Path to configTable.ini file", type=str, default=default_config_path)
 
    # Parse arguments corresponding to parameters
    parser.add_argument('-o', '--output', help="Name of output table", type=str)
    parser.add_argument('-oc', '--outCol', help='Name/Index of columns present in output table. By default, all columns will be displayed', type=str)

    parser.add_argument('-e' '--compare', help='Name/Index of columns compared during row fusion' ,\
        type=str, default=default_compare_columns)
    parser.add_argument('-r' '--conserve', help='Name/Index of columns whose value will be conserved during row fusion' ,\
        type=str, default=default_conserve_columns)

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()


    # parse config with user selection
    config_param = configparser.ConfigParser(inline_comment_prefixes='#')
    config_param.read(Path(args.config))

    # Parameters introduced in the execution replace those in the .ini file
    if args.e__compare:
        config_param.set('Parameters', 'ComparedColumns', args.e__compare)

    if args.r__conserve:
        config_param.set('Parameters', 'ConservedColumns', args.r__conserve)

    if args.output:
        config_param.set('Parameters', 'OutputName', args.output)
    
    if args.outCol:
        config_param.set('Parameters', 'OutputColumns', args.outCol)


    # logging debug level. By default, info level
    if args.infile:
        log_file = outfile = args.infile[:-4] + '_log.txt'
        log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
    
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