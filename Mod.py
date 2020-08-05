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
import pandas as pd
import numpy as np
import re
import csv
from pygoslin.parser.Parser import LipidParser
from multiprocessing import Pool, cpu_count

import pdb

###################
# Local functions #
###################

def readLipidList(infile_path):
    '''
    Description: readLipidList reads goslin lipid list in csv format. The fucntion
    will create a list with lipid names in first columns, and all its synonyms,
    iterating over all rows.

    Input:
        - infile_path: Path where csv lipid list can be found
    Output:
        - lipid_list: List of strings with all lipid names
    '''

    with open(infile_path, 'r') as infile:

        line = infile.readline() # skip title row
        lipid_reader = csv.reader(infile, delimiter=',', quotechar='"')

        lipid_list = []
        for row in lipid_reader:
            if len(row) > 0 and len(row[0]) > 0:
                lipid_list.append(row[0])
                lipid_list.extend([r for r in row[6:] if len(r) > 0])
    
    return lipid_list


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


def removeRow(df_row, name_col_index, regex):
    '''
    Input:
        - df_row: Row of pandas dataframe received in the apply
        - name_col_index: Index of compound name
        - regex: Regular expression applied to compound name string
    Output:
        - True if regex is found and False otherwise.
    '''
    
    compound_name = df_row.iat[name_col_index]

    return bool(re.search(regex, compound_name))


def splitCompoundField(compound_name, regex_sep):
    '''
    Description: parseCompoundList receives the string present in the field
    corresponding to the compound name. It splits the string into the different
    compounds names using separator regex given by the user.

    Input:
        - compound_name: String with compound names that may be separated
        - regex_sep: String with regular expression with possible separators
    
    Output:
        - compound_list: List of different compounds
        - separator: String with the separator found for the compounds
    '''

    match_sep = re.search(regex_sep, compound_name)

    # If there is match, make split. Else, return a list with the single compound
    if match_sep is not None:
        separator = match_sep.group()
        compound_list = re.split(separator, compound_name)
    
    else:
        separator = ''
        compound_list = [compound_name]
    
    return compound_list, separator



def parseCompoundList(compound_list, regex, replace):
    '''
    Description: parseCompoundList applies the transformation to each compound in the list

    Input:
        - compound_list: List of strings, with the name of the compounds
        - regex: String with the regular expression applied
        - replace: String with the value that will replace the recognized pattern
    Output:
        - List of strings with the parsed compounds
    '''

    return [re.sub(regex, replace, compound.strip()) for compound in compound_list]


def parserCompound(df_row, name_col_index, regex, replace, regex_sep):
    '''
    Description: parserCompound apply a regular expression to the row received, replacing the
    recognized patter by the value given in 'replace'
    
    Input:
        - df_row: Series (row) from the pandas dataframe
        - name_col_index: Integer with the name column index (0-based)
        - regex: String with the regular expression
        - replace: String with the value that replaces the pattern
        - regex_sep: Separator between compounds within a field
    Output:
        - df_row_out: Series (row) transformed
    '''

    # Get compound name
    compound_name = df_row.iat[name_col_index]

    # Obtain a list with different compounds names within the field (one or more)
    compound_list, separator = splitCompoundField(compound_name, regex_sep)

    # Apply the transformation to each compound of the field (they are in the list)
    parsed_compound_list = parseCompoundList(compound_list, regex, replace)

    parsed_compound_name = separator.join(parsed_compound_list)
    
    df_row_out = df_row.copy()
    df_row_out.iat[name_col_index] = parsed_compound_name

    return df_row_out


def parserTable(df_row, name_col_index, regex_sep):
    '''
    Description: The function applies each regular expression to the row received (df_row). These
    regular expressions are stored in config_regex. With a loop we iterate over each regex and its
    replacement. Then, we call parserCompound function to apply  each transformation.

    Input:
        - df_row: Pandas series with the row being processed
        - name_col_index: Integer corresponding to the index of the column name
        - regex_sep: String corresponding to the compound separator in name field
    Output:
        - df_row_out: Pandas series corresponding to the modified row
    '''

    df_row_out = df_row.copy()

    # Iterate over each regular expresion
    for regex_section in config_regex.sections():

        # It gives a list with the two values of the section
        regex, replace = [config_regex[regex_section][option] for option in config_regex[regex_section]._options()]

        # Parse the row using parserCompound
        df_row_out = parserCompound(df_row_out, name_col_index, regex, replace, regex_sep) 
    
    return df_row_out


def isPeptide(aa_list):
    '''
    Description: isPeptide receives a list of strings obtained from splitting a compound name by the aminoacid
    separator given by the user. If all strings are aminoacids it is returned True. Otherwise, it returns False.
    '''

    aminoacids = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", \
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]

    return np.all([aa in aminoacids for aa in aa_list])


def sortPeptides(df_row, name_col_index, aa_sep):
    '''
    Description: sortPeptides split compound name by aminoacid separator and asserts that it is
    a peptide. If so, peptide aminoacids are sorted and joined by the separator. It is returned
    the processed row, with the sequence of both, the sorted peptide, and the original one, separated
    by the label '#####'
    '''

    df_row_out = df_row.copy()

    compound_name = df_row.iat[name_col_index] # Extract compound name

    # Split peptide in aminoacids
    aa_list, separator = splitCompoundField(compound_name, aa_sep)

    if isPeptide(aa_list):
        aa_sorted_list = sorted(aa_list)
        compound_name_out = separator.join(aa_sorted_list) + '#####' + separator.join(aa_list)
        df_row_out.iat[name_col_index] = compound_name_out

    return df_row_out


def lipidPreProcess(compound):
    '''
    Description: Compound names is parsed removing the information that cannot be processed by Goslin.
    
    Input:
        - compound_out: String with the name of the compounds
    Output:
        - compound_out: String with the parsed name of the compounds
    '''

    # Remove information between [ ]
    compound_out = re.sub(r'\[\w+\]', '', compound)

    # Remove 'i-' 'a-' information. There must be a '(' or '/' before it, and a number after it: TG(i-13:0/i-14:0/8:0)
    compound_out = re.sub(r'(?<=[\(/])[ia]-(?=\d)', '', compound_out)

    return compound_out


def getHeaderGroup(lipid):
    '''
    Description: Get head group from lipid, and in case it begins with Lyso, replace
    by L
    
    Input:
        - lipid: Lipid object from Goslin library
    Output:
        - head_group: String with head group
    '''

    head_group = lipid.lipid.head_group
    head_group = re.sub(r'(?i)^Lyso', 'L', head_group)

    return head_group


def getCarbonAtoms(lipid):
    '''
    Input:
        - lipid: Lipid object from Goslin library
    Output:
        - Integer with total number of carbon atoms in fatty acids present in lipids
    '''

    fa_num_carbon_list = np.array([fa.num_carbon for fa in lipid.lipid.fa_list])
    return np.sum(fa_num_carbon_list)


def getDoubleBonds(lipid):
    '''
    Input:
        - lipid: Lipid object from Goslin library
    Output:
        - Integer with total number of double bonds in fatty acids present in lipids
    '''

    fa_num_double_bonds_list = np.array([fa.num_double_bonds for fa in lipid.lipid.fa_list])
    return np.sum(fa_num_double_bonds_list)


def getFaBondType(lipid):
    '''
    Description: getFaBondTyppe recognizes if any of the fatty acids in lipid has a plasmanyl (O-)
    or plasmenyl (P-) bond type. O- is associated to number 2 and P- to number 3 in Goslin. 1 would
    be ESTHER bond and 0 is undefined.

    Input:
        - lipid: Lipid object from Goslin library
    Output:
        - 'O-' if any 2
          'P-' if any 3
          '' else    
    '''

    fa_bond_type_list = np.array([fa.lipid_FA_bond_type.value for fa in lipid.lipid.fa_list])

    if 2 in fa_bond_type_list:
        return 'O-'
    
    elif 3 in fa_bond_type_list:
        return  'P-'
    
    else:
        return ''


def lipidCandidate(compound_name, lipid):
    '''
    Description: lipidCandidate returns True if lipid string is at the beginning
    of the compound name string.
    '''

    if re.match('^'+lipid, compound_name):
        return True
    
    else:
        return False



def isGoslinLipid(compound_name, lipid_list):
    '''
    Description: isGoslinLipid returns True is compound_name is in goslin lipid_list. Otherwise
    it will return False
    '''

    return np.any([lipidCandidate(compound_name, lipid) for lipid in lipid_list])


def parserLipidCompound(compound, lipid_list):
    '''
    Description: parserLipidCompound parses the lipid name using Goslin library. First,
    it applies a filter using a regular expression, to avoid creating a goslin lipid object with all
    compounds, which will raise an error in most of the cases, making the code too slow. If the filter
    is passed, it creates the goslin lipid object, and then extract the required information from it.

    Input:
        - compound: String the compound name
    Output:
        - compound_out: String with the parsed (if possible) compound name
    '''
    
    # Apply a filter using LipidRegex in parameters.ini to make code faster
    if not isGoslinLipid(compound, lipid_list):
        return compound

    # Pre-process compound name so that it can be recognized by goslin
    pre_proc_compound = lipidPreProcess(compound)

    try:
        # Create (if possible) Goslin lipid object from compound name
        lipid_parser = LipidParser()
        lipid = lipid_parser.parse(pre_proc_compound)

        # If lipid has no fatty  acid, return compound
        if lipid.lipid.fa_list == []:
            return compound

        # Get head group
        head_group = getHeaderGroup(lipid)

        # Get total number of carbon atoms
        n_carbon_atoms = getCarbonAtoms(lipid)

        # Get total number of double bonds
        n_double_bonds = getDoubleBonds(lipid)

        # Get FA bond type (plasmanyl/plasmenyl)
        fa_bond_type = getFaBondType(lipid)

        # Build lipid name using extracted information
        compound_out = head_group + '(' + fa_bond_type + str(n_carbon_atoms) + ':' + str(n_double_bonds) + ')'

        return compound_out

    except:
        # If gosling cannot parse the compound, it is returned without any change
        return compound


def parserLipidTable(df_row, name_col_index, regex_sep, lipid_list):
    '''
    Description: parserLipidTable is the function applied over the pandas dataframe.
    It receives each row, and process the compound name of lipids.

    Input:
        - df_row: Pandas series corresponding to a row of the dataframe
        - name_col_index: Integer corresponding to index of column name
        - regex_sep: String with compound separator in name field
        - lipid_regex: String with the regular expression used to identify lipids to be processed
    Output:
        - df_row_out: Output pandas series with the parsed name
    '''

    df_row_out = df_row.copy()

    # Get compound name
    compound_name = df_row.iat[name_col_index]

    # Obtain a list with different compounds names within the field (one or more)
    compound_list, separator = splitCompoundField(compound_name, regex_sep)

    # Parse each compound of the list, and join by the separator
    parsed_compound_list = [parserLipidCompound(compound, lipid_list) for compound in compound_list]
    parsed_compound_name = separator.join(parsed_compound_list)

    df_row_out.iat[name_col_index] = parsed_compound_name

    return df_row_out


def subProcessFuncRegex(df_i, name_col_index, regex_sep, aa_sep):
    '''
    Description: Function executed using starmap() method from Pool. It receives chunks
    of dataframe, each of which are processed using apply. Dataframes are processed using
    regular expressions in the case of parserTable, and peptide aminoacids are sorted
    alphabetically in the case of sortPeptide.
    '''
    
    df_i_out = df_i.copy()

    # Parse dataframe using parserTable, which will iterates over the regular expressions
    df_i_out = df_i_out.apply(func=parserTable, axis=1, args=(name_col_index, regex_sep))

    # Process peptides names, so that peptides with equal compoisition can be fused during the fusion
    df_i_out = df_i_out.apply(func=sortPeptides, axis = 1, args=(name_col_index, aa_sep))

    return df_i_out


def subProcessFuncLipid(df_i, name_col_index, regex_sep, lipid_list):
    '''
    Description: Function executed using starmap() method from Pool. It receives chunks of dataframe,
    processed using apply with parserLipidTable, which is going to process lipids with Goslin.
    '''
    
    df_i_out = df_i.copy()

    # Parse name of compound lipids in the dataframe using goslin (parserLipidTable)
    df_i_out = df_i_out.apply(func=parserLipidTable, axis=1, args=(name_col_index, regex_sep, lipid_list))

    return df_i_out


def sortIndexByName(df, name_col_index):
    '''
    Description: sortIndexByName sorts row indexes of pandas dataframe (df) using as reference
    the alphanumeric order of compound names. In this sense, indexes of rows with equal compound
    name will be together. Thanks to this, the comparison among rows during the fusion can be
    made faster, as we can compare one row with the following ones until the compound names 
    are different.

    Input:
        - df: Pandas dataframe with all data
        - name_col_index: integer with the position of the column containing compound names
    Output:
        - index_out: List of integers with the ordered indexes
    '''

    # Extract compound name column and row indexes
    compound_name = df.iloc[:, name_col_index]
    row_index = df.index.values

    # Create a list of tuples, where each tuple contains the row index and the compound name (lower)
    index_name_tuple_list = [(index, name.lower().replace('-', '')) for index, name in zip(row_index, compound_name)]

    # Sort tuples by second element, which is the compound name
    sorted_index_name_tuple_list = sorted(index_name_tuple_list, key=lambda x: x[1])

    # Extract indexes sorted by compound name
    index_out  = [index for index, name in sorted_index_name_tuple_list]

    return index_out


def getIndex(element, column_names):
    '''
    Description: getIndex receives an element (0-based index or name) associated to a column table. The function
    will return the 0-based index of that column.

    Input:
        - element: string with the 0-based index or name of the column
        - column_names: strings numpy array with the names of the columns
    Output:
        - out_element: Integer with the 0-based index of the column
    '''

    if re.match(r'^\d+$', element):
        out_element = int(element)
    
    else:
        out_element = int(np.where(column_names == element)[0][0])
    
    return out_element


def getColumnList(string_indexes, name_col_index, column_names):
    '''
    Description: getColumnList returns an array with the column indexes that are indicated
    by the user, without repeating the one with compound names

    Input:
        - string_indexes: String with comma-separated numbers, corresponding to the 0-based
        column index
        - name_col_index: Integer corresponding to the 0-based index of compound name column
        - column_names: strings numpy array with the names of the columns
    Output:
        - index_array: Numpy array with integers corresponding to 0-based integers
    '''

    # Split comma-separated elements given by the user
    index_array = np.array([getIndex(element.strip(), column_names) for element in string_indexes.split(',')])

    # Remove index corresponding to column name
    index_array = index_array[index_array != name_col_index]

    return index_array

def getTagColumnNames(tag_str, column_names):
    '''
    Input:
        - tag_str: String containing tag containing column names separated by comma (',')
        - column_names: Pandas series containing all column names in infile
    Output:
        - tag_list: Numpy array containing tag containing column names
    '''

    tag_list = np.array([tag.strip() for tag in tag_str.split(',') if tag.strip() in column_names])

    return tag_list


def combineConservedValues(conserved_values_i, conserved_values_j):
    '''
    Description: combineConservedValues receives two arrays of equal length. Each element 
    corresponds to a value of the row in one conserved field. The function will combine
    values of each field among both rows, i and j. In other words, value of j is added to i 
    (separated by ' // '), unless value of j is present in i.

    Input:
        - conserved_values_i: Array with values of upper row
        - conserved_values_j: Array with values of downer row
    Output:
        - conserved_out: List of strings with combined values for each field
    '''

    # Convert both arrays in str numpy  arrays
    conserved_values_i_str = np.array(conserved_values_i, dtype=str)
    conserved_values_j_str = np.array(conserved_values_j, dtype=str)

    # Zip function creates a list of tuple, where each tuple has i and j value for a given field.
    # It is iterated over different fields
    conserved_out = [field_i if field_j == '' else field_j if field_i == '' \
        else field_i + ' // ' + field_j if field_j not in field_i.split(' // ') else field_i \
        for field_i, field_j in zip(conserved_values_i_str, conserved_values_j_str)]
    
    return conserved_out


def fuseTable(df, name_col_index):
    '''
    Description: fuseTable compare each row of the pandas dataframe with downer rows, following
    the row order given by 'sorted_index'. If both rows have the same values for the comparing
    columns (given by the user), the downer row is removed/dropped. Values of the downer row present
    in conserved columns (also given by the user) are added to the values of the upper row
    (combineConservedValues). Using 'sorted_index' we iterate the rows following the alphanumeric 
    order of compound names, so we avoid comparing each row with all the rest. When upper row has
    finds a downer row with different compound name, we jump to next upper row.
    '''

    # Store column names in a numpy array
    column_names = np.array(df.columns, dtype=str)

    # Get information of the column with compound names (Index and column name)
    name_col_name = column_names[name_col_index]

    # Get lists with indexes and names corresponding to columns that will be compared
    compared_columns = getColumnList(args.compareCol, name_col_index, column_names)
    compared_columns_names = column_names[compared_columns]

    # Get lists with indexes and names corresponding to columns whose values will be conserved in the fusion
    conserved_columns = getColumnList(args.conserveCol, name_col_index, column_names)
    conserved_columns_names = column_names[conserved_columns]

    # Get list with the name of the columns containing tags of the compound
    tag_columns_names = getTagColumnNames(args.tagCol, df.columns)

    # Add columns containing tags to the set of conserving columns
    conserved_columns_names = np.concatenate((conserved_columns_names, tag_columns_names))

    # Get list with index of rows in the order over which they are going to be iterated
    sorted_index = sortIndexByName(df, name_col_index)

    # List with dropped-row indexes, to avoid their iteration
    removed_index = []

    # Loop upper row index
    for i, index_i in enumerate(sorted_index):

        # If upper row index was dropped, move to the next
        if index_i in removed_index:
            continue

        # Get string with compound name and array string with compared values, corresponding to upper row
        compound_name_i_original = str(df.at[index_i, name_col_name])
        compared_values_i = np.array([df.at[index_i, col] for col in compared_columns_names], dtype=str)

        # Variable to control when to move upper row to the next index (when i and j compound names are different)
        # continue_i = False

        # Loop over down row index
        for j, index_j in enumerate(sorted_index):
            
            # If down row index is below upper row index, or if it was dropped, move to the next
            if (j <= i) or (index_j in removed_index):
                continue

            # Get string with compound name and array string with compared values, corresponding to down row
            compound_name_i = compound_name_i_original  # In case compound_name_i was modified in previous iterations (peptides), we use the original
            compound_name_j = str(df.at[index_j, name_col_name])
            compared_values_j = np.array([df.at[index_j, col] for col in compared_columns_names], dtype=str)

            # If both are peptides, take only sorted name
            if ('#####' in compound_name_i) and ('#####' in compound_name_j):
                compound_name_i = compound_name_i.split('#####')[0]
                compound_name_j = compound_name_j.split('#####')[0]

            # If compound names are different, break downer loop, and move upper index to the next
            if compound_name_i.lower().replace('-', '') != compound_name_j.lower().replace('-', ''):
                # continue_i = True
                break

            # If any value in comparing field is different, move downer index to next
            elif np.any(compared_values_i != compared_values_j):
                continue

            # If all values in comparing field are the same, make the row fusion
            else:
                # Get conserved values
                conserved_values_i = ['' if pd.isna(df.at[index_i, col]) else df.at[index_i, col] for col in conserved_columns_names]
                conserved_values_j = ['' if pd.isna(df.at[index_j, col]) else df.at[index_j, col] for col in conserved_columns_names]
                
                # Combine conserved values of i and j, and store them in upper row
                df.loc[index_i, conserved_columns_names] = combineConservedValues(conserved_values_i, conserved_values_j)

                # Drop downer row and append its index to the removed_index list
                df.drop(axis=0, index=index_j, inplace=True)
                removed_index.append(index_j)

        # If i and j compound names were different, continue_i is True, so upper index move to the next
        # if continue_i:
        #    continue
    
    return df


def originalPeptides(df, name_col_index):
    '''
    Description: originalPeptides takes column name from the dataframe, and iterated over all of them. If label
    '#####' is in the compound name, it is taken the part of the name after it. By this way, it is taken the
    original peptide name, removing the sorted one.
    '''

    # Take compound name column as a list
    all_names = list(df.iloc[:, name_col_index])

    # Iterate over all names. If ##### is present, split and take the second part
    compound_name_out = [compound_name if '#####' not in compound_name else compound_name.split('#####')[1] for compound_name in all_names]

    # Restore compound name column
    df.iloc[:, name_col_index] = compound_name_out

    return df


def getOutFileName(infile):
    '''
    Description: getOutFileName generate a string with the output filename. If user did not specified the output
    filename in OutputName of parameters.ini, the name will be 'parsed_'+infile. Besides, the function tests if
    the output file name given by the user has extension. If not, .xls will be used.
    '''

    outfile_name, outfile_ext = os.path.splitext(config_param['Parameters']['OutputName'])

    if not outfile_name:
        outfile = 'mod_' + infile
    
    elif not outfile_ext:
        outfile = outfile_name + '.xls'
    
    else:
        outfile = outfile_name + outfile_ext
    
    return outfile


def getOutColumns(column_names):
    '''
    Description: getOutColumns receives a numpy array with the names of the columns. It returns
    the name of those columns selected by the user.
    '''

    out_columns = config_param['Parameters']['OutputColumns']

    if out_columns: 
        out_columns_index = getColumnList(string_indexes=out_columns, name_col_index=None, column_names=column_names)
        out_columns_name = column_names[out_columns_index]
    
    else:
        out_columns_name = column_names
    
    return out_columns_name


def writeDataFrame(df, path):
    '''
    Description: The function will write the padas dataframe in a 
    result folder using pandas.to_excel method.

    Input:
        - df: Pandas dataframe that is going to be written
        - path: Infile path given by the user
    '''

    # Build path of output file, including dir and name
    # output_path = os.path.join(os.path.dirname(path), 'results')
    output_path = os.path.join(os.path.dirname(path))
    output_filename = getOutFileName(os.path.basename(path))
    output_path_filename = os.path.join(output_path, output_filename)

    # Get output columns
    out_columns_name = getOutColumns(np.array(df.columns))


    logging.info(f'Writing output table in {output_path_filename}')

    # If result folder does not exist, create it
    if not os.path.isdir(output_path):
        logging.info('Creating "results" folder')
        os.mkdir(output_path)
    
    # Handle errors in exception case
    try:
        df.to_excel(output_path_filename, index=False, columns=out_columns_name)
    
    except:
        log_str = f'Error when writing {str(Path(output_path_filename))}'
        logging.info(log_str)
        
        # Log error class and message
        exctype, value = sys.exc_info()[:2]
        log_str = f'{exctype}: {value}'
        logging.info(log_str)

        sys.exit()
    
    log_str = f'{str(Path(output_path_filename))} was written'
    logging.info(log_str)

    return True


##################
# Main functions #
##################

def main(args):
    '''
    Main function
    '''

    # Number of cores. This should be a user variable
    n_cores = cpu_count() - 1 
    logging.info(f"Using {n_cores} cores")

    # Store compound name column index
    name_col_index = args.column      # Column containing compound names (0-based)
    header_index = args.row           # Row at which table begins (0-based)

    regex_remove = config_param['Parameters']['RemoveRow']
    regex_sep = config_param['Parameters']['Separator']
    aa_sep = config_param['Parameters']['AminoAcidSeparator']

    # Read goslin lipid list in csv format
    lipid_list = readLipidList(args.liplist)

    # Read table as pandas data frame
    df = readInfile(args.infile, header_index)

    # Remove rows given by RemoveRow regular expression
    logging.info("Removing rows identified by RemoveRow parameter")
    remove_row_bool = df.apply(func=removeRow, axis=1, args=(name_col_index, regex_remove))
    df.drop(axis=0, index=np.where(remove_row_bool)[0], inplace=True)


    # Split dataframe so that each one is processed by one core
    df_split = np.array_split(df, n_cores)

    # Create list of tuples. Each tuple contains arguments received by subProcessFunction
    subprocess_args = [(df_i, name_col_index, regex_sep, aa_sep) for df_i in df_split]

    with Pool(n_cores) as p: 
        logging.info(f'Applying regular expression from {os.path.basename(args.regex)} and sorting peptide aminoacids alphabetically')
        result = p.starmap(subProcessFuncRegex, subprocess_args)
        df_processed = pd.concat(result)

    # Fuse rows with the same value for the selected columns. Make a fusion before goslin lipid processing to make the code faster
    logging.info(f'Collapsing rows after metabolite name parsing')
    fused_df1 = fuseTable(df_processed, name_col_index)

    # For each peptide, take only the original part
    logging.info(f'Peptides post-processing (replace alphabetically sorted name by one of the original names)')
    fused_df1 = originalPeptides(fused_df1, name_col_index)


    # Split dataframe so that each is processed by one core
    df_split = np.array_split(fused_df1, n_cores)

    # Create list of tuples. Each tuple contains arguments received by subProcessFunctionLipid
    subprocess_args = [(df_i, name_col_index, regex_sep, lipid_list) for df_i in df_split]

    with Pool(n_cores) as p:
        logging.info(f'Parsing lipid names using Goslin')
        result = p.starmap(subProcessFuncLipid, subprocess_args)
        df_processed2 = pd.concat(result)

    # Fuse rows with the same value for the selected columns
    logging.info(f'Collapsing rows after lipid processing')
    fused_df2 = fuseTable(df_processed2, name_col_index)

    # Write output dataframe
    writeDataFrame(fused_df2, args.infile)

    

if __name__ == '__main__':
    
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Mod',
        epilog='''
        Example:
            python Mod.py
        
        '''
    )

    # Set default values
    default_config_regex = os.path.join(os.path.dirname(__file__), "config/configMod/regex.ini")
    default_config_parameters = os.path.join(os.path.dirname(__file__), "config/configMod/parameters.ini")
    default_lipid_list = os.path.join(os.path.dirname(__file__), "Data/goslinLipidList.csv")

    default_column_index = 5    # Column containing compound names (0-based)
    default_header_index = 0    # Row at which table begins (0-based)
    default_compare_column = "0, 5" # Columns used to compare rows during fusion
    default_conserve_column = "1"   # Columns whose value is conserved during fusion
    default_tag_column = "Food, Drug, Halogenated, Microbial"


    # Parse arguments corresponding to input, .ini and lipid list paths
    parser.add_argument('-i', '--infile', help='Path to input file', required=True, type=str)
    parser.add_argument('-re', '--regex', help='Path to custom regex.ini file', default=default_config_regex, type=str)
    parser.add_argument('-pr', '--param', help='Path to custom parameters.ini file', default=default_config_parameters, type=str)
    parser.add_argument('-ll', '--liplist', help='Path to goslin lipid list csv file', default=default_lipid_list, type=str)


    # Parameters corresponding to parameters.ini  (not all of them are in parameters.ini file)
    parser.add_argument('-n', '--name', help='Name of output table', type=str)
    parser.add_argument('-p', '--column', help='Column index of compound names (0-based)', default=default_column_index, type=int)
    parser.add_argument('-r', '--row', help='Row of column headers, at which the table starts (0-based)', default=default_header_index, type=int)
    parser.add_argument('-s', '--separator', help='Characters used to separate compound within a field (accept regex)', type=str)
    parser.add_argument('-aas', '--aa_separator', help='Characters used to separate aminoacids in peptides', type=str)
    parser.add_argument('-rm', '--rmRow', help='Regular expression used in Name field to identify rows that will be dropped', type=str)
    parser.add_argument('-cmp', '--compareCol', help='Index/Name of columns (0-based) that will be compared to make the row fusion (e.g. 0,5)',\
        default=default_compare_column, type=str)
    parser.add_argument('-cns', '--conserveCol', help='Index/Name of columns (0-based) whose values will be conserved during the row fusion (e.g. 1)',\
        default=default_conserve_column, type=str)
    parser.add_argument('-tag', '--tagCol', help='Name of columns containing tags of the compounds (e.g. FoodTag, DrugTag). Their values will be conserved',\
        default=default_tag_column, type=str)
    parser.add_argument('-oc', '--outCol', help='Index/Name of columns present in output table. By default, all columns will be displayed (e.g. 0,2,5)', type=str)

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()


    # parse config with regular expressions
    config_regex = configparser.ConfigParser(inline_comment_prefixes='#')
    config_regex.read(Path(args.regex))

    # parse config with parameters
    config_param = configparser.ConfigParser(inline_comment_prefixes='#')
    config_param.read(Path(args.param))

    # Parameters introduced in the execution replace those in the .ini file
    if args.name is not None:
        config_param.set('Parameters', 'OutputName', str(args.name))

    if args.separator is not None:
        config_param.set('Parameters', 'Separator', str(args.separator))

    if args.aa_separator is not None:
        config_param.set('Parameters', 'AminoAcidSeparator', str(args.aa_separator))

    if args.rmRow is not None:
        config_param.set('Parameters', 'RemoveRow', str(args.rmRow))

    if args.outCol is not None:
        config_param.set('Parameters', 'OutputColumns', str(args.outCol))


    # logging debug level. By default, info level
    if args.infile:
        log_file = outfile = os.path.splitext(args.infile)[0] + '_log.txt'
        log_file_debug = outfile = os.path.splitext(args.infile)[0] + '_log_debug.txt'
    
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