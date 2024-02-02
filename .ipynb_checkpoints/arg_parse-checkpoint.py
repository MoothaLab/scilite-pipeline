import argparse
import csv
import os
import pandas as pd

import utils

def map_seq_to_edit_type(seq_to_type_fp):
    ''' Create a dictionary that maps target read sequences to edit types
    
    Args:
        - seq_to_type_fp: str, the path to the csv file containing the 
                               target read sequence-to-edit type mapping
    
    Returns:
        - edit_type_map: dict {str -> str}, maps a target sequence to an edit type
    '''

    seq_to_type_df = pd.read_csv(seq_to_type_fp, 
                                 skiprows=1, 
                                 names=['edit_type', 'r1_sequence'],
                                 header=None)
    seq_to_type_df['r1_sequence'] = [x.strip().upper() for x in seq_to_type_df['r1_sequence'].str.upper()]
    seq_to_type_df['edit_type'] = [x.strip() for x in seq_to_type_df['edit_type']]
    edit_type_map = dict(seq_to_type_df.to_numpy())

    return edit_type_map

def map_bc1_to_cond(bc1_to_cond_fp):
    ''' Create a dictionary that maps barcode 1 to sample conditions
    
    Args:
        - bc1_to_cond_fp: str, the path to the csv file containing the
                               barcode 1-to-condition mapping
    
    Returns:
        - condition_map: dict {str -> str}, maps a barcode 1 to a condition
    '''

    bc1_to_cond_df = pd.read_csv(bc1_to_cond_fp, 
                                 skiprows=1, 
                                 names=['barcode_1', 'condition'],
                                 header=None)
    condition_map = dict(bc1_to_cond_df[['barcode_1', 'condition']].to_numpy())

    return condition_map

def get_bc2(bc2_list_fp):
    ''' Create a list that contains all unique barcode 2 
    
    Args:
        - bc2_list_fp: str, the path to the csv file containing the 
                       unique barcode 2
                       (first row contains column names)
    
    Returns: 
        - bc2_list: list of strs, a list containing all unique barcode 2
    '''

    bc2_df = pd.read_csv(bc2_list_fp, 
                         skiprows=1, 
                         names=['barcode_2'],
                         header=None)
    bc2_list = list(bc2_df['barcode_2'].to_numpy())

    return bc2_list

def get_bc_tracer(bc_tracer_list_fp):
    ''' Create a list that contains all unique tracer barcodes 
    
    Args:
        - bc_tracer_list_fp: str, the path to the csv file containing the 
                          unique tracer barcodes
                          (first row contains column names)
    
    Returns: 
        - bc_tracer_list: list of strs, a list containing all unique tracer barcodes
    '''

    bc_tracer_df = pd.read_csv(bc_tracer_list_fp, 
                            skiprows=1, 
                            names=['barcode_tracer'],
                            header=None)
    bc_tracer_list = list(bc_tracer_df['barcode_tracer'].to_numpy())

    return bc_tracer_list

def process_barcode_and_seq_mappings(fq_process_kwargs, exp_kwargs):
    ''' Read user-specified barcode 1-condition mapping and target sequence-edit type mappings
        as well as generate substitution-tolerated mapping for barcodes and target sequences
    
    Args:
        - fq_process_kwargs: dict, keyword arguments specific for processing fastq reads
                                   and need to be appended with more arguments
        - exp_kwargs: dict, keyword arguments specific for the experiment setup
                            that have already been parsed
     
     Returns:
        - fq_process_kwargs: dict, keyword arguments specific for the experiment setup
    '''

    # read and parse the barcode 1-to-condition mapping and the barcode 2 list
    condition_map = map_bc1_to_cond(exp_kwargs['bc1_to_cond_map_fpath'])
    bc2_list = get_bc2(exp_kwargs['bc2_list_fpath'])

    # generate mappings between substituted barcodes and the original barcode
    # the substituted barcodes are within a predefined edit distance of the 
    # original barcode (default edit distance limit = 1 base pair)
    bc1_ed_map = utils.gen_sub_map(list(condition_map.keys()), n_max_subs=1)
    bc2_ed_map = utils.gen_sub_map(bc2_list, n_max_subs=1)

    # append those mappings to fq_process_kwargs since they are required
    # in the fastq read processing
    fq_process_kwargs['condition_map'] = condition_map
    fq_process_kwargs['bc1_ed_map'] = bc1_ed_map
    fq_process_kwargs['bc2_ed_map'] = bc2_ed_map
    fq_process_kwargs['edit_type_map'] = exp_kwargs['edit_type_map']

    return fq_process_kwargs

def parse_one_tracer_construct(line, edit_to_seq_map, in_path):
    ''' Parse one line of experimental set-up regarding an tracer construct
    
    Args:
        - line: list of strs, the line in the parsed .CSV file regarding the tracer construct
        - edit_to_seq_map: dict {str -> str}, maps an edit type to its corresponding 
                                              expected target sequence
        - in_path: str, the path to the directory where the experimental setting .CSV file
                        as well as all other input files should be located
     
     Returns:
        - construct_info: dict, keyword arguments specific for the tracer construct being parsed
    '''
    construct_info = {}
    construct_name = line[0].strip()
    tracer_seq = edit_to_seq_map[construct_name]

    construct_info['construct_name'] = construct_name
    tracer_bc_fns = [ fn.strip() for fn in line[1].strip().split(',') ]
    tracer_bc_pos = [ int(pos.strip()) for pos in line[2].strip().split(',') ]
    n_tracer_bc = len(tracer_bc_fns)
    for bc_idx in range(n_tracer_bc):
        construct_info['bc_tracer_'+str(bc_idx)+'_pos'] = [tracer_bc_pos[bc_idx*2], tracer_bc_pos[bc_idx*2+1]]
        tracer_bc_fn =  tracer_bc_fns[bc_idx]
        if os.path.isabs(tracer_bc_fn) == False:
            tracer_bc_fn = os.path.abspath(os.path.join(in_path, tracer_bc_fn))
        construct_info['bc_tracer_'+str(bc_idx)+'_list'] = utils.gen_sub_map(get_bc_tracer(tracer_bc_fn), n_max_subs=1, immutable_pos=[])
    construct_info['n_tracer_bc'] = n_tracer_bc

    construct_info['bc1_pos'] = [ int(pos.strip()) for pos in line[3].strip().split(',') ]
    construct_info['bc2_pos'] = [ int(pos.strip()) for pos in line[4].strip().split(',') ]
    construct_info['umi_pos'] = [ int(pos.strip()) for pos in line[5].strip().split(',') ]

    if_used_in_filtering = line[6].strip().upper()
    construct_info['umi_for_filtering'] = 'Y' in if_used_in_filtering

    construct_info['seq_len_set'] = [len(tracer_seq)]
    construct_info['seq_ed_map'] = utils.gen_sub_map([tracer_seq], n_max_subs=1, immutable_pos=[])

    return construct_info
    
def parse_one_edit_construct(line, edit_to_seq_map):
    ''' Parse one line of experimental set-up regarding an target-edit construct
    
    Args:
        - line: list of strs, the line in the parsed .CSV file regarding the tracer construct
        - edit_to_seq_map: dict {str -> str}, maps an edit type to its corresponding 
                                              expected target sequence
                             
     Returns:
        - construct_info: dict, keyword arguments specific for the target-edit construct being parsed
    '''
    construct_info = {}
    construct_name = line[0].strip()

    construct_info['construct_name'] = construct_name
    construct_info['bc1_pos'] = [ int(pos.strip()) for pos in line[3].strip().split(',') ]
    construct_info['bc2_pos'] = [ int(pos.strip()) for pos in line[4].strip().split(',') ]
    construct_info['umi_pos'] = [ int(pos.strip()) for pos in line[5].strip().split(',') ]

    if_used_in_filtering = line[6].strip().upper()
    construct_info['umi_for_filtering'] = 'Y' in if_used_in_filtering

    immutable_pos = [int(pos.strip()) for pos in line[1].strip().split(',')] if line[1].strip() else []
    edit_types = [edit.strip() for edit in line[2].split(',') if edit != '']
    edit_seqs = []
    seq_lens = []
    for edit_type in edit_types:
        edit_seq = edit_to_seq_map[edit_type]
        edit_seqs.append(edit_seq)
        seq_lens.append(len(edit_to_seq_map[edit_type]))

    construct_info['immutable_pos'] = immutable_pos
    construct_info['edit_types'] = edit_types
    construct_info['seq_len_set'] = list(set(seq_lens))
    construct_info['seq_ed_map'] = utils.gen_sub_map(edit_seqs, n_max_subs=1, immutable_pos=immutable_pos)
    return construct_info

def parse_constructs(construct_arg_lines, edit_to_seq_map, in_path):
    ''' Parse lines of experimental set-up information regarding an target-edit or tracer construct
    
    Args:
        - construct_arg_lines: list of strs, the lines in the parsed .CSV file regarding all constructs
        - edit_to_seq_map: dict {str -> str}, maps an edit type to its corresponding 
                                              expected target sequence
        - in_path: str, the path to the directory where the experimental setting .CSV file
                        as well as all other input files should be located
     
     Returns:
        - all_construct_info: dict, keyword arguments specific for all constructs in this experiment set-up
        - tracer_edit_types: list of strs, contains all edit types that are actually tracer construct names
    '''
    all_construct_info = {}
    tracer_edit_types = []
    # set up boolean flags to indicate whether current line being parsed
    # should belong to a tracer construct or a target-edit sequence construct
    # since they have different fields
    reading_tracer_construct = False
    reading_edit_construct = False
    for line in construct_arg_lines:
        if line[0] == 'Lineage Tracing Sequence Name':
            reading_tracer_construct = True
            reading_edit_construct = False
            continue
        elif line[0] == 'Sequence Name':
            reading_edit_construct = True
            reading_tracer_construct = False
            continue
        elif line[0] == '':
            reading_tracer_construct = False
            reading_edit_construct = False
            continue
        elif line[0] == 'Edit Type to Mutant Type Map':
            return all_construct_info, tracer_edit_types

        # parse the current line as an tracer or target-edit sequence construct
        if reading_tracer_construct == True:
            construct_info = parse_one_tracer_construct(line, edit_to_seq_map, in_path)
            construct_name = construct_info['construct_name']
            all_construct_info[construct_name] = construct_info
            tracer_edit_types.append(construct_name)
        elif reading_edit_construct == True:
            construct_info = parse_one_edit_construct(line, edit_to_seq_map)
            construct_name = construct_info['construct_name']
            all_construct_info[construct_name] = construct_info

    return all_construct_info, tracer_edit_types

def parse_mappings(map_arg_lines):
    ''' Parse the edit to mutant type map and the mutant to heteroplasmy map
    
    Args:
        - map_arg_lines: list of strs, each str contains one arg_name - arg_val pair
    
    Returns:
        - mut_type_map: dict {str -> str}, contains the parsed edit_type -> mutant_type pairs
        - het_map: dict {str -> list of strs}, contains the parsed heteroplasmy_type -> mutant_type pairs
    '''

    mut_type_map = {}
    het_map = {}
    reading_mut_type_map = False
    reading_het_map = False
    for line in map_arg_lines:

        # decide whether the current line being parse contains 
        # edit-to-mut or mut-to-het information
        if line[0] == 'Mutant Type':
            reading_mut_type_map = True
            continue
        elif line[0] == 'Heteroplasmy Type':
            reading_het_map = True
            continue
        elif line[0] == '':
            reading_mut_type_map = False
            reading_het_map = False
            continue
        elif line[0] == 'Edit Type to Mutant Type Map' or line[0] == 'Mutant Type to Heteroplasmy Map':
            continue

        # append one entry to the map if applicable
        if reading_mut_type_map == True:
            edit_types = [edit.strip() for edit in line[1].split(',') if edit != '']
            mut_type = line[0].strip()
            for edit_type in edit_types:
                mut_type_map[edit_type] = mut_type
        elif reading_het_map == True:
            mut_types = [mut.strip() for mut in line[1].split(',')  if mut != '']
            het_type = line[0].strip()
            het_map[het_type] = mut_types

    return mut_type_map, het_map

def parse_exp_args(exp_arg_names, exp_arg_lines, in_path):
    ''' Parse the experiment-specific arguments
    
    Args:
        - exp_arg_names: list of strs, each str contains one arg_name
        - exp_arg_lines: list of strs, each str contains one arg_name - arg_val pair
        - in_path: str, the path to the folder where all input files are located
    
    Returns:
        - args: dict {str -> *}, contains the parsed arg_name - arg_val pairs
    '''

    arg_vals = [line[1] for line in exp_arg_lines] # first line of the CSV file contains column names
    for i in range(4):
        fp = arg_vals[i]
        if os.path.isabs(fp) == False:
            fp = os.path.abspath(os.path.join(in_path, fp))
        arg_vals[i] = fp
    exp_args = dict(zip(exp_arg_names, arg_vals))

    exp_args['required_fields'] = [field.strip() for field in exp_args['required_fields'].split(',')]
    return exp_args

def parse_input_arg_file(arg_csv_fp):
    ''' Parse the CSV file containing required input arguments
    
    Args:
        - arg_csv_fp: str, the path to the CSV file containing required arguments
    
    Returns:
        - exp_kwargs: dict {str -> *}, containing the parsed arguments
    '''

    # read the CSV file into memory
    lines = []
    with open(arg_csv_fp, newline='') as in_file:
        argreader = csv.reader(in_file, delimiter=',')
        for row in argreader:
            lines.append(row)

    arg_csv_fp = os.path.abspath(arg_csv_fp)
    in_path = '/'.join(arg_csv_fp.split('/')[:-1])

    exp_arg_names = ['seq_to_edit_map_fpath', 
                     'bc1_to_cond_map_fpath', 'bc2_list_fpath',
                     'fastq_directory_path', 'fastq_suffix',
                     'required_fields']
    exp_arg_lines = lines[1 : len(exp_arg_names)+1]
    other_lines = lines[len(exp_arg_names)+1:]
    exp_args = parse_exp_args(exp_arg_names, exp_arg_lines, in_path)
    exp_kwargs = exp_args

    edit_to_seq_map = map_seq_to_edit_type(exp_kwargs['seq_to_edit_map_fpath'])
    edit_type_map = {v: k for (k, v) in edit_to_seq_map.items() }
    all_construct_info, tracer_edit_types = parse_constructs(other_lines, edit_to_seq_map, in_path)
    exp_kwargs['edit_type_map'] = edit_type_map
    exp_kwargs['construct_info'] = all_construct_info
    exp_kwargs['tracer_edit_types'] = tracer_edit_types

    # parse the following lines containing the edit_type -> mut_type, 
    # and the mutant_type -> heteroplasmy_type mappings
    mut_type_map, het_map = parse_mappings(other_lines)
    exp_kwargs['het_map'] = het_map
    exp_kwargs['mut_type_map'] = mut_type_map

    mut_types_for_filtering = []
    for (cid, cinfo) in all_construct_info.items():
        if cinfo['umi_for_filtering'] == True:
            etypes = cinfo['edit_types']
            for etype in etypes:
                mut_types_for_filtering.append(mut_type_map[etype])
    exp_kwargs['mut_types_for_filtering'] = list(set(mut_types_for_filtering))

    return exp_kwargs

def create_parser():
    ''' Create a parser for processing command-line arguments
    
    Returns:
        - parser: ArgumentParser obj, the parser for commandline arguments
    '''

    parser = argparse.ArgumentParser(description='Process reads and calculate heterplasmy for a single SCI-LITE experiment')

    parser.add_argument('-e', '--experimentSettingFile', 
                        help='File path to the experiment argument spreadsheet (saved as .csv file)',
                        required=True)
    parser.add_argument('-o', '--outputPath', 
                        help='Path to where the output files where be located (defaults to the directory containing the experiment argument spreadsheet)',
                        default=None)
    parser.add_argument('-m', '--maxCells', 
                        help='Specifies the maximum number of cells allowed to be preserved in the filtering-by-knee-plot process. The final resulting heteroplasmy spreadsheet may have fewer cells than specified, but never more (default=no limit)',
                        default=None)

    parser.add_argument('-p', '--isPairedEnd', 
                        help='Specifies that the input FASTQs contain paired-end reads (default=False)',
                        action='store_true')
    parser.add_argument('-n', '--isWellIDFromRead', 
                        help='Specifies that the well ID should be extracted from FASTQ read entries (the last space-delimited element from the first line of each FASTQ record), instead of the FASTQ file name (default=False)',
                        action='store_true')
    parser.add_argument('-r', '--reverseComplementIDs', 
                        help='Specifies the ID sequences extracted from the FASTQ reads (barcode 1s, barcode 2s, etc.) should be reverse-complemented to match their expected sequences defined in bc1_map.csv, bc2_list.csv, etc. (default=False)',
                        action='store_true')
    parser.add_argument('-s', '--saveIntermediates',
                        help='Specifies that intermediate read-level and UMI-level spreadsheets will also be saved as output files (default=False)',
                        action='store_true')
    parser.add_argument('-k', '--noPlotKnee', 
                        help='Specifies that the knee plot used for filtering out cells with too few UMIs will NOT be visualized and saved as an output file – “knee_plot.pdf” (default=False)',
                        action='store_true')
    parser.add_argument('-c', '--compressOutputFiles', 
                        help='Specifies that output .csv files will be compressed into .gz files (default=False)',
                        action='store_true')
    parser.add_argument('--phase_range', type=int, default=1,
                        help='Specifies the maximum length of heterogeneity spacers used in the library construction (default is 1 for no heterogeneity spacers). For example, the experiments in the SCI-LITE paper that used heterogeneity spacers required --phase_range=8, while those that didn\'t used the default --phase_range=1.')
    
    parser.add_argument('--read_dataframe', type=str, default=None,
                    help='Specifies the path to an intermediate read-level dataframe to skip the FASTQ read parsing step. Useful for tuning the cell filtering by knee plot (default=None).')

    return parser
