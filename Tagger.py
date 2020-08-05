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
from multiprocessing import Pool, cpu_count

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


def readFoodTable(path):
    '''
    Input:
        - path: String containing the path with the food containing table
    
    Output:
        - food_list: Numpy array of strings containing all food compounds in the table
    '''

    logging.info(f"Reading Food Table: {path}")

    df = pd.read_csv(path, header=0, sep="\t", dtype=str)
    food_list = np.array(df.iloc[:, 0].drop_duplicates(keep='first', inplace=False))

    return food_list


def getNameColumnIndex(column_names):
    '''
    Input:
        - column_names: Pandas series containing the names of the columns in the infile table
    
    Output:
        - An integer indicating the position of the Name column
    '''

    return int(np.where(column_names == "Name")[0][0])


def foodTaggerBatch(df, food_list):
    '''
    Input:
        - df: Pandas dataframe containing batch of the total content
        - food_list: String Numpy Array with all food compounds extracted from the database
    
    Output:
        - df: Pandas dataframe with the "Food" tag added in a new column
    '''

    # Get numpy array with compound names in the dataframe
    compound_names = np.array(df.loc[:, 'Name']) 

    # Tag corresponding compounds using food list
    food_tag_from_db = ["Food" if compound in food_list else "" for compound in compound_names]

    # Tag compounds that fits regular expression
    food_tag_from_regex = ["Food" if re.search('^[Ee]nt-', compound) else "" for compound in compound_names]

    # Combine Food tags
    food_tag = ["Food" if "Food" in tag else "" for tag in zip(food_tag_from_db, food_tag_from_regex)]

    # Add Food tag column to the dataframe
    name_column_index = getNameColumnIndex(df.columns)
    df.insert(name_column_index+1, "Food", food_tag, True)
    
    return df


def foodTagger(df, n_cores):
    '''
    Input:
        - df: Pandas dataframe containing the whole content
        - n_cores: Integer indicating the number of cores used in the multiprocessing
    
    Output:
        - df_out: Pandas dataframe with the whole content and the added tag
    '''
    
    logging.info("Start food tagging")

    # Split dataframe so that each batch is processed by one core
    df_split = np.array_split(df, n_cores)
    
    # Get numpy array with food compounds in database
    food_list = readFoodTable(args.foodList)
        
    # Create list of tuples. Each tuple contains arguments received by foodTaggerBatch in each subprocess
    subprocess_args = [(df_i, food_list) for df_i in df_split]
        
    with Pool(n_cores) as p:

        logging.info("Tagging food compounds")
        result = p.starmap(foodTaggerBatch, subprocess_args)
        df_out = pd.concat(result)
    
    logging.info("Finished food tagging")

    return df_out


def readDrugTable(path):
    '''
    Input:
        - path: String containing the path to the Drug database
    
    Output:
        - drug_list: String Numpy Array containing drugs name extracted from the database
    '''

    logging.info(f"Reading Drug Table: {path}")

    # Import drug table as a pandas dataframe
    df = pd.read_csv(path, header=0, sep="\t", dtype=str)

    # Extract drug list name from df as a numpy array
    drug_list = np.array(df.iloc[:, 0])

    return drug_list


def drugTaggerBatch(df, drug_list):
    '''
    Input:
        - df: Pandas dataframe containing a batch of the whole infile dataframe
        - drug_list: String Numpy Array containing all drug compounds in the database
    
    Output:
        - df: Pandas dataframe with the drug tag added in a new column
    '''

    # Get numpy array with compound in input table
    compound_names = np.array(df.loc[:, 'Name']) 

    # Tag corresponding compounds
    drug_tag = ["Drug" if compound in drug_list else "" for compound in compound_names]

    # Add Drug tag column to the dataframe
    name_column_index = getNameColumnIndex(df.columns)
    df.insert(name_column_index+1, "Drug", drug_tag, True)
    
    return df


def drugTagger(df, n_cores):
    '''
    Input:
        - df: Pandas dataframe containing the whole infile content
        - n_cores: Integer indicating the number of cores used in the multiprocessing

    Output: df_out: Pandas dataframe containing the whole infile content with the Drug Tag 
    added in a new column
    '''

    logging.info("Start drug tagging")

    # Split dataframe so that each batch is processed by one core
    df_split = np.array_split(df, n_cores)

    # Get numpy array with drug list
    drug_list = readDrugTable(args.drugList)

    # Create list with parameters received by each drugTaggerBatch function in each subprocess
    subprocess_args = [(df_i, drug_list) for df_i in df_split]

    with Pool(n_cores) as p:

        logging.info("Tagging drug compounds")
        result = p.starmap(drugTaggerBatch, subprocess_args)
        df_out = pd.concat(result)
    
    logging.info("Finished drug tagging")

    return df_out


def halogenatedTaggerBatch(df, halogen_regex):
    '''
    Input:
        - df: Pandas dataframe corresponding to a batch of the infile table
        - halogen_regex: String corresponding to the regular expression used to identify halogenated compounds
    
    Output:
        - df: Pandas dataframe with Halogenated tag added in a new column
    '''

    # Get numpy array with compound in input table
    compound_names = np.array(df.loc[:, 'Name']) 

    # Tag corresponding compounds
    halogenated_tag = ["x" if re.search(halogen_regex, compound) else "" for compound in compound_names]

    # Add Drug tag column to the dataframe
    name_column_index = getNameColumnIndex(df.columns)
    df.insert(name_column_index+1, "Halogenated", halogenated_tag, True)
    
    return df


def halogenatedTagger(df, n_cores):
    '''
    Input:
        - df: Pandas dataframe containing the whole content present in infile table
        - n_cores: Integer indicating the number of cores used in the multiprocessing
    
    Output:
        - df_out: Pandas dataframe with halogenated tag added in a new column
    '''

    logging.info("Start halogenated compounds tagging")

    # Split dataframe so that each batch is processed by one core
    df_split = np.array_split(df, n_cores)

    # Get string with the regular expression used to identify halogenated compounds
    halogen_regex = config_param.get('Parameters', 'HalogenatedRegex')

    # Create list with parameters received by halogenatedTaggerBatch in each subprocess
    subprocess_args = [(df_i, halogen_regex) for df_i in df_split]

    with Pool(n_cores) as p:

        logging.info("Tagging halogenated compounds")
        result = p.starmap(halogenatedTaggerBatch, subprocess_args)
        df_out = pd.concat(result)
    
    logging.info("Finished halogenated compounds tagging")

    return df_out    


def getOutputFilename():
    '''
    Output:
        - filename: String containing the name of the output file
    '''

    filename = config_param.get('Parameters', 'OutputName')

    if not filename:
        filename = 'tagged_' + os.path.basename(args.infile)

    elif not os.path.splitext(filename)[1]:
        filename += '.xls'
    
    return filename


def getOutputColumns(df_columns):
    '''
    Input:
        - df_columns: Pandas series containing the name of the columns in the output table
    
    Output:
        - selected_columns: List of strings with the name of the columns selected by the user
    '''

    selected_columns = config_param.get('Parameters', 'OutputColumns')

    if selected_columns:
        selected_columns = [column.strip() for column in selected_columns.split(',') if column.strip() in df_columns]
    
    else:
        selected_columns = list(df_columns)
    
    return selected_columns


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
    output_columns = getOutputColumns(df.columns)

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



#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''

    # Number of cores used
    n_cores = cpu_count() - 1
    logging.info(f"Using {n_cores} cores")

    # Read infile
    df = readInfile(args.infile, 0)

    # Check user selection
    if re.search('(?i)true', config_param['TagSelection']['Food']):
        df = foodTagger(df, n_cores)
    
    if re.search('(?i)true', config_param['TagSelection']['Drug']):
        df = drugTagger(df, n_cores)
    
    if re.search('(?i)true', config_param['TagSelection']['Halogenated']):
        df = halogenatedTagger(df, n_cores)
    
    # Export dataframe as excel file
    writeDataFrame(df, args.infile)
    


if __name__=="__main__":

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Tagger',
        epilog='''
        Example:
            python Tagger.py
        
        '''
    )

    # Set default values
    default_config_path = os.path.join(os.path.dirname(__file__), "config/configTagger/configTagger.ini")
    default_food_list_path = os.path.join(os.path.dirname(__file__), "Data/food_database.tsv")
    default_drug_list_path = os.path.join(os.path.dirname(__file__), "Data/drug_database.tsv")
    # default_microbial_compound_list_path = os.path.join(os.path.dirname(__file__), "Data/TaggerData/microbial_compound_database.tsv")

    # Parse arguments corresponding to path files
    parser.add_argument('-i', '--infile', help="Path to input file", type=str, required=True)
    parser.add_argument('-c', '--config', help="Path to configTagger.ini file", type=str, default=default_config_path)
    parser.add_argument('-fL', '--foodList', help="Path to food compounds list", type=str, default=default_food_list_path)
    parser.add_argument('-dL', '--drugList', help="Path to drug compounds list", type=str, default=default_drug_list_path)
    # parser.add_argument('-mL', '--microbialList', help="Path to microbial compounds list", type=str, default=default_microbial_compound_list_path)

    parser.add_argument('-o', '--output', help="Name of output table", type=str)
    parser.add_argument('-oc', '--outCol', help='Name of columns present in output table. By default, all columns will be displayed', type=str)

    # Parser arguments indicating which tags are going to be added
    parser.add_argument('-f', '--food', help="Add food tag to compounds", action='store_true', default=False)
    parser.add_argument('-d', '--drug', help="Add drug tag to compounds", action='store_true', default=False)
    parser.add_argument('-m', '--microbial', help="Add 'microbial compound' tag to compounds", action='store_true', default=False)
    parser.add_argument('-ha', '--halogenated', help="Add 'halogenated compound' tag to compounds", action='store_true', default=False)

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    # parse config with user selection
    config_param = configparser.ConfigParser(inline_comment_prefixes='#')
    config_param.read(Path(args.config))

    # Parameters introduced in the execution replace those in the .ini file
    if args.food:
        config_param.set('TagSelection', 'Food', str(args.food))

    if args.drug:
        config_param.set('TagSelection', 'Drug', str(args.drug))

    if args.microbial:
        config_param.set('TagSelection', 'MicrobialCompound', str(args.microbial))

    if args.halogenated:
        config_param.set('TagSelection', 'Halogenated', str(args.halogenated))
    
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