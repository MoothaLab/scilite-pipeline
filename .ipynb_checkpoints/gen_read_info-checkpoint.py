import numpy as np
import pandas as pd
import pysam
import re, glob, os, itertools

import utils

def find_fastq_paths(directory_path, fastq_suffix):
    ''' Find the path to all fastq files with the same suffix under a specified directory
    
    Args:
        - directory_path: str, the path of the directory containing all desired fastq files
        - fastq_suffix: str, the shared suffix of desired fastq files
        
    Returns:
        - fastq_fps: list of strs, the list containing paths to all desired fastq files
    '''
    
    # append a forward slash to the end of the directory path if it is not there yet
    if directory_path[-1] != '/' and fastq_suffix[0] != '/':
        directory_path += '/'
    
    # find the filenames of all desired fastq files under the directory
    fastq_fns = [fn for fn in glob.glob(directory_path + '*' + fastq_suffix) if 'Undetermined' not in fn]

    test_if_well_id_in_fn = ( len(fastq_fns[0].split('_')) > 2 )
    
    # extract the alphabetical and the numerical parts of the well id indpendently
    if test_if_well_id_in_fn:
        well_ids = [os.path.basename(fn).split('_')[0] for fn in fastq_fns]
        well_id_alphs = [''.join([ele for ele in well_id if not ele.isdigit()]) for well_id in well_ids]
        well_id_nums = [int(''.join([ele for ele in well_id if ele.isdigit()])) for well_id in well_ids]
        del well_ids

        # sort the fastq filenames 
        # by the alphabetical part of the well id first
        # then the numerical part of the well id
        fastq_fps = [fn for _, _, fn in sorted(zip(well_id_alphs, well_id_nums, fastq_fns),
                                               key=lambda x: (x[0], x[1]))]
    else:
        fastq_fps = fastq_fns

    return fastq_fps

def get_pe_filepath(r1_fpath, r1_dif_fname, r2_dif_fname):
    ''' Get filenames for a pair of (paired-end) fastq files
    
    Args: 
        - r1_fpath: str, the file path for the fastq file containing read1
        - r1_dif_fname: str, the substring of r1_fpath that is different from its corresponding read2 file path
        - r2_dif_fname: str, the substring of r2_fpath that is different from its corresponding read1 file path
        
    Returns:
        - r1_fpath: str, the file path for the fastq file containing read1
        - r2_fpath: str, the file path for the fastq file containing read2
    '''

    r2_fpath = r1_fpath.replace(r1_dif_fname, r2_dif_fname)

    return r1_fpath, r2_fpath

def extract_id_info(read, read_stat, 
                    if_paired_end_read, if_reverse_complement_id_seq, offset, 
                    bc1_pos, bc2_pos, umi_pos):
    ''' Extract the read-specific information from a read
    
    Args:
        - read: PySAM FastqProxy object, a single entry in the fastq file that contains the read-specific info
                                         assumed to be read 1 in single-read sequencing and read 2 in paired-end
        - read_stat: dict {str -> str}, contains information regarding a single read
        - if_paired_end_read: bool, indicating if the experiment produces paired-ends
                                    if set to TRUE, the id sequences (barcode1/2/UMI) should be located
                                    on read2 and thus have their locations unaffected by heterogeneity
                                    spacers being spiked in
        - if_reverse_complement_id_seq: bool, reverse-complement the barcodes if set to TRUE
        - offset: int, the length of heterogeneity spacer getting spiked into the current read
                       as the results, the positions of id sequences need to be shifted by this length
        - bc1_pos: list of ints, the 0-index positions of barcode 1 on the read
        - bc2_pos: list of ints, the 0-index positions of barcode 2 on the read
        - umi_pos: list of ints, the 0-index positions of umi sequence on the read
    
    Returns:
        - read_stat: dict {str -> str}, now contains the barcodes, umi and read sequence of a single read
    '''

    # in paired-read experiments, the id sequences (barcode1/2/UMI) should be located on read 2 and thus
    # have their locations unaffected by heterogeneity spacers spiked into the beginning of read 1
    if if_paired_end_read == True:
        offset = 0

    # locate the id sequences on the read
    read_stat['bc1'] = read.sequence[offset+bc1_pos[0] : offset+bc1_pos[1]]
    read_stat['bc2'] = read.sequence[offset+bc2_pos[0] : offset+bc2_pos[1]]
    read_stat['umi'] = read.sequence[offset+umi_pos[0] : offset+umi_pos[1]]

    # reverse-complement the id sequences 
    # if the provided reference barcodes are expected to be their reverse-complements
    if if_reverse_complement_id_seq == True:
        read_stat['bc1'] = utils.gen_reverse_complement(read_stat['bc1'])
        read_stat['bc2'] = utils.gen_reverse_complement(read_stat['bc2'])
        read_stat['umi'] = utils.gen_reverse_complement(read_stat['umi'])

    return read_stat

def find_construct(read_stat, target_seq, info_seq, construct_info, tracer_edit_types,
                   edit_type_map, bc1_ed_map, bc2_ed_map, 
                   if_paired_end, if_reverse_complement_id_seq,
                   phase_range=1):
    ''' Find out which construct does the read belong to
    
    Args:
        - read_stat: dict {str -> str}, contains information regarding a single read
        - target_seq: PySAM FastqProxy object, a single entry in the fastq file that contains the targeted sequence
                                               assumed to be read 1 in both single-read and paired-read experiments
        - info_seq: PySAM FastqProxy object, a single entry in the fastq file that contains the read-specific info
                                             assumed to be read 1 in single-read sequencing and read 2 in paired-end
        - construct_info: dict, contains construct-specific information such as where the id sequences will be located
                          see below for its structure:
                                {'seq_ed_map': seq_ed_map (dict {str -> str}, maps the substituted targeted sequence (key) 
                                                                                to the original (value)),
         Stefano Monti                        'bc1_pos': bc1_pos (list of ints, the 0-index position of barcode 1 on the read),
                                 'bc2_pos': bc2_pos (list of ints, the 0-index position of barcode 2 on the read),
                                 'umi_pos':  umi_pos (list of ints, the 0-index position of umi sequence on the read),
                                 'bc_tracer_pos': bc_tracer_pos (list of ints, the 0-index position of the tracer barcode on the read 
                                                           if there are tracer constructs present in the experiment),
                                 'bc_tracer_list': bc_tracer_list (dict {str -> str}, contains all expected tracer barcodes),             
                                 'immutable_pos': immutable_pos (list of ints, the 0-index position of targeted edit on the read sequence)
                                }
        - tracer_edit_types: list of strs, contains all edit types belonging to any tracer construct
        - edit_type_map: dict {str -> str}, maps target sequence to edit types
        - bc1_ed_map: dict {str -> str}, maps the substituted barcode 1 (key) to the original (value)
        - bc2_ed_map: dict {str -> str}, maps the substituted barcode 2 (key) to the original (value)
        - if_paired_end: bool, indicates if the experiment is paired-ended
        - if_reverse_complement_id_seq: bool, indicates if the id sequences like barcode 1/2/UMI should be reverse-complemented 
                                                to match their corresponding references
        - phase_range: int, indicates the number of insertions at the start of the read to tolerate when trying to match the edit type. 
                            Generally only greater than 1 to accommodate libraries prepared with heterogeneity spacers. 
                            (default: 1)
        
    Returns:
        - read_stat: dict {str -> str}, now contains the barcodes, umi and read sequence of a single read
    '''

    found_stat = {'bc1': 'UNKNOWN',
                  'bc2': 'UNKNOWN'}

    # iterate through each construct involved in the experiment
    for construct_name in construct_info.keys():

        # get the corresponding dictionary containing construct-specific information
        cinfo = construct_info[construct_name]

        # try different sequence lengths & heterogeneity spacer lengths
        for seq_len, offset in itertools.product(cinfo['seq_len_set'], range(phase_range)):

            # locate the target sequence assuming it belongs to the current construct
            cur_seq = target_seq.sequence[offset:offset+seq_len]
            # check if the located target sequence is within a reasonable edit distance
            # of the pre-defined expected target sequences
            cur_seq_corr = cinfo['seq_ed_map'].get(cur_seq, 'UNKNOWN')

            # if the located target sequence is deemed as expected, gather other read-level
            # information from the current read according to the assumptions of the current construct
            if cur_seq_corr != 'UNKNOWN':

                read_stat['seq'] = cur_seq
                read_stat['seq_corr'] = cur_seq_corr
                read_stat['het_spacer_offset'] = offset
                read_stat['edit_type'] = edit_type_map.get(cur_seq_corr, 'UNKNOWN')

                # extract bc1/2/umi from info_seq
                read_stat = extract_id_info(info_seq, read_stat, 
                                            if_paired_end, if_reverse_complement_id_seq,
                                            offset,
                                            cinfo['bc1_pos'], cinfo['bc2_pos'], cinfo['umi_pos'])

                # extract the tracer barcode if the current read belongs to an tracer construct
                if read_stat['edit_type'] in tracer_edit_types:
                    all_bc_tracers = []
                    for bc_idx in range(cinfo['n_tracer_bc']):
                        bc_tracer_name = 'bc_tracer_'+str(bc_idx)
                        tracer = target_seq.sequence[offset+cinfo[bc_tracer_name+'_pos'][0] : offset+cinfo[bc_tracer_name+'_pos'][1]]
                        tracer = cinfo[bc_tracer_name+'_list'].get(tracer, 'UNKNOWN')
                        all_bc_tracers.append(tracer)
                        if if_reverse_complement_id_seq == True:
                            tracer = utils.gen_reverse_complement(tracer)
                    read_stat['bc_tracer'] = '+'.join(all_bc_tracers)
                else:
                    read_stat['bc_tracer'] = 'UNKNOWN'

                # stop the search here and return extracted read-level information
                return read_stat

            else:
                test_stat = {}
                test_stat = extract_id_info(info_seq, read_stat, 
                                            if_paired_end, if_reverse_complement_id_seq,
                                            offset,
                                            cinfo['bc1_pos'], cinfo['bc2_pos'], cinfo['umi_pos'])
                bc1_corr =  bc1_ed_map.get(test_stat['bc1'], 'UNKNOWN')
                if bc1_corr != 'UNKOWN':
                    found_stat['bc1'] = test_stat['bc1']

                bc2_corr =  bc2_ed_map.get(test_stat['bc2'], 'UNKNOWN')
                if bc2_corr != 'UNKOWN':
                    found_stat['bc2'] = test_stat['bc2']

    read_stat['bc1'] = found_stat['bc1']
    read_stat['bc2'] = found_stat['bc2']

    return read_stat

def parse_one_read_pe(well_id, r1, r2, **kwargs):
    ''' Parse one single read produced in a paired-end sequencing experiment
        and extract relevant information from it
    
    Args:
        - well_id: str,
        - r1: PySAM FastqProxy object, a single entry in the fastq file that contains read 1
        - r2: PySAM FastqProxy object, a single entry in the fastq file that contains read 2
        - **kwargs: dict, other keyword arguments specific for the experiment setup
                    see parse_all_fastq() for specific information
                        
    Returns:
        - read_stat: dict {str -> str}, contains information regarding a single read that is being parsed
        - read_id: str, contains the identifier for the parsed read
    '''

    # initialize the dictionary that will contain extracted information for the current read
    read_stat = {}
    if kwargs['if_well_id_from_read'] == True:
        well_id = r1.comment.split(' ')[-1]
    read_stat['well_id'] = well_id
    # print(read_stat['well_id'])
    read_stat['r1'] = r1.sequence
    read_stat['r2'] = r2.sequence

    # construct the read-level identfier using the well id, read 1 and read 2 sequences
    read_id = well_id + r1.sequence + r2.sequence

    # initialize construct-specific fields to UNKNOWN
    for stat_name in ['seq', 'seq_corr', 'edit_type',
                      'bc1', 'bc2', 'umi', 
                      'bc_tracer', 'het_spacer_offset']:
        read_stat[stat_name] = 'UNKNOWN'

    # try to extract sequence edit types & id sequences according to different constructs
    read_stat = find_construct(read_stat, r1, r2, 
                               kwargs['construct_info'], kwargs['tracer_edit_types'],
                               kwargs['edit_type_map'], kwargs['bc1_ed_map'], kwargs['bc2_ed_map'],
                               kwargs['if_paired_end'], kwargs['if_reverse_complement_id_seq'],
                               kwargs['phase_range'])

    # map the potentially substituted sequences to their corrected ones
    # set to UNKNOWN if the corresponding corrected sequences cannot be identified
    read_stat['bc1_corr'] = kwargs['bc1_ed_map'].get(read_stat['bc1'], 'UNKNOWN')
    read_stat['bc2_corr'] = kwargs['bc2_ed_map'].get(read_stat['bc2'], 'UNKNOWN')

    # map the corrected barcode 1 to the experimental condition
    # set to UNKNOWN if the edit type / condition cannot be identified
    read_stat['condition'] = kwargs['condition_map'].get(read_stat['bc1_corr'], 'UNKNOWN')

    for key in ['seq', 'bc1', 'bc2']:
        read_stat.pop(key, None)

    return read_stat, read_id

def parse_one_read_se(well_id, r1, **kwargs):
    ''' Parse one single read produced in a single-end sequencing experiment
        and extract relevant information from it
    
    Args:
        - well_id: str,
        - r1: PySAM FastqProxy object, a single entry in the fastq file that contains read 1
        - **kwargs: dict, other keyword arguments specific for the experiment setup
                    see parse_all_fastq() for specific information
    
    Returns:
        - read_stat: dict {str -> str}, contains information regarding a single read that is being parsed
        - read_id: str, contains the identifier for the parsed read
    '''

    # initialize the dictionary that will contain extracted information for the current read
    read_stat = {}
    if kwargs['if_well_id_from_read'] == True:
        well_id = r1.comment.split(' ')[-1]
    # print(read_stat['well_id'])
    read_stat['well_id'] = well_id
    read_stat['r1'] = r1.sequence

    # construct the read-level identfier using the well id and read sequence
    read_id = well_id +  r1.sequence

    # initialize construct-specific fields to UNKNOWN
    for stat_name in ['seq', 'seq_corr', 'edit_type',
                      'bc1', 'bc2', 'umi', 
                      'bc_tracer', 'bc_tracer', 'het_spacer_offset']:
        read_stat[stat_name] = 'UNKNOWN'

    # try to extract sequence edit types & id sequences according to different constructs
    read_stat = find_construct(read_stat, r1, r1, 
                               kwargs['construct_info'], kwargs['tracer_edit_types'],
                               kwargs['edit_type_map'], kwargs['bc1_ed_map'], kwargs['bc2_ed_map'],
                               kwargs['if_paired_end'], kwargs['if_reverse_complement_id_seq'],
                               kwargs['phase_range'])

    # map the potentially substituted sequences to their corrected ones
    # set to UNKNOWN if the corresponding corrected sequences cannot be identified
    read_stat['bc1_corr'] = kwargs['bc1_ed_map'].get(read_stat['bc1'], 'UNKNOWN')
    read_stat['bc2_corr'] = kwargs['bc2_ed_map'].get(read_stat['bc2'], 'UNKNOWN')

    # map the corrected barcode 1 to the experimental condition
    # set to UNKNOWN if the edit type / condition cannot be identified
    read_stat['condition'] = kwargs['condition_map'].get(read_stat['bc1_corr'], 'UNKNOWN')

    for key in ['seq', 'bc1', 'bc2']:
        read_stat.pop(key, None)

    return read_stat, read_id

def parse_one_fastq_pe(r1_fpath, read_dict={}, record_count=0, 
                       r1_dif_fname='_R1_', r2_dif_fname='_R2_', 
                       **kwargs):

    ''' Parse one pair of fastq files produced in a paired-end sequencing experiment
        and extract relevant information from it
        
    Args:
        - r1_fpath: str, the file path for the fastq file containing read 1
        - read_dict: dict {str -> dict}, contains all parsed reads so far
        - record_count: int, the number of parsed reads so far
        - r1_dif_fname: str, the substring of r1_fpath that is different from its corresponding read2 file path
        - r2_dif_fname: str, the substring of r2_fpath that is different from its corresponding read1 file path
        - **kwargs: dict, other keyword arguments specific for the experiment setup
                    see parse_all_fastq() for specific information 
    
    Returns:
        - read_dict: dict {str -> dict}, now also contains reads from the current fastq file pair
        - record_count: int, the number of parsed reads after the current fastq file pair was processed
    '''

    # extract the well id from the file name
    if re.match(r'[A-H]{1,2}\d{1,2}_', os.path.basename(r1_fpath)):
        well_id = os.path.basename(r1_fpath).split('_')[0]
    else:
        well_id = os.path.basename(r1_fpath).split('.')[1]

    # get the paths for both fastq files in the current pair
    r1_fpath, r2_fpath = get_pe_filepath(r1_fpath, r1_dif_fname, r2_dif_fname)

    # iterate through each read
    with pysam.FastqFile(r1_fpath) as fq1_in, pysam.FastqFile(r2_fpath) as fq2_in:

        for r1, r2 in zip(fq1_in, fq2_in):

            # extract information from each read
            read_stat, read_id = parse_one_read_pe(well_id, r1, r2, **kwargs)

            # increment the read count for the current read
            # if it is a duplicate for some read that has already been parsed
            try:
                read_dict[read_id]['read_count'] += 1
            except KeyError:
                read_stat['read_count'] = 1
                read_dict[read_id] = read_stat

            # progress bar
            if record_count and record_count%1e6 == 0:
                print('Processed {!s} reads.'.format(record_count))
            record_count += 1

    return read_dict, record_count

def parse_one_fastq_se(r1_fpath, read_dict={}, record_count=0, **kwargs):
    
    ''' Parse one pair of fastq files produced in a single-end sequencing experiment
        and extract relevant information from it
    
    Args:
        - r1_fpath: str, the file path for the fastq file containing read 1
        - read_dict: dict {str -> dict}, contains all parsed reads so far
        - record_count: int, the number of parsed reads so far
        - **kwargs: dict, other keyword arguments specific for the experiment setup
                    see parse_all_fastq() for specific information 
    
    Returns:
        - read_dict: dict {str -> dict}, now also contains reads from the current fastq file pair
        - record_count: int, the number of parsed reads after the current fastq file pair was processed
    '''

    # extract the well id from the file name
    if re.match(r'[A-H]{1,2}\d{1,2}_', os.path.basename(r1_fpath)):
        well_id = os.path.basename(r1_fpath).split('_')[0]
    else:
        well_id = os.path.basename(r1_fpath).split('.')[1]

    # iterate through each read
    with pysam.FastqFile(r1_fpath) as fq1_in:

        for r1 in fq1_in:

            # extract information from each read
            read_stat, read_id = parse_one_read_se(well_id, r1, **kwargs)

            # increment the read count for the current read
            # if it is a duplicate for some read that has already been parsed
            try:
                read_dict[read_id]['read_count'] += 1
            except KeyError:
                read_stat['read_count'] = 1
                read_dict[read_id] = read_stat

            # progress bar
            if record_count and record_count%5e6 == 0:
                print('Processed {!s} reads.'.format(record_count))
            record_count += 1

    return read_dict, record_count

def parse_all_fastq(fastq_fps, **kwargs):
    ''' Parse a list of of fastq files and extract relevant information from it
        
    Args:
        - fastq_fps: list of strs, the list containing paths to all desired fastq files
        - **kwargs: dict, other keyword arguments specific for the experiment setup
                    see below for its structure:
                        {'edit_type_map': edit_type_map (dict {str -> str}, maps target sequence to edit types),
                         'condition_map': condition_map (dict {str -> str}, maps barcode 1s to conditions),
                         'construct_info': construct_info (dict {str -> dict}, contains construct-specific info,
                         'tracer_edit_types': tracer_edit_types (list of strs, contains all edit types belonging to any tracer construct),
                         'bc1_ed_map': bc1_ed_map (dict {str -> str}, maps the substituted barcode 1 (key) to the original (value)), 
                         'bc2_ed_map': bc2_ed_map (dict {str -> str}, maps the substituted barcode 2 (key) to the original (value)),
                         'if_paired_end': if_paired_end (bool, indicates if the experiment is paired-ended),
                         'if_reverse_complement_id_seq': if_reverse_complement_id_seq (bool, indicates if the id sequences like 
                                                                                       barcode 1/2/UMI should be reverse-complemented
                                                                                       to match their corresponding references),
                         'fastq_directory_path': fastq_directory_path (str, the path to the directory that contains all fastq files),
                         'fastq_suffix': fastq_suffix (str, the suffix of file names shared by all fastq files involved in the experiment)
                        }
    
    Returns:
        - read_df: dict {str -> dict}, contains reads from all parsed fastq files
    '''

    read_dict = {}
    record_count = 0

    # parse each fastq file in the list
    for fastq_fp in fastq_fps:
        if kwargs['if_paired_end'] == True:
            read_dict, record_count =  parse_one_fastq_pe(fastq_fp, read_dict, record_count, **kwargs)
        else:
            read_dict, record_count =  parse_one_fastq_se(fastq_fp, read_dict, record_count, **kwargs)

    # strip the intermediate read-level identfiers 
    # and convert the dictionary to list to decrease memory usage
    read_dict = list(read_dict.values())
    read_dict = {col_name: [read[col_name] for read in read_dict] for col_name in read_dict[0]}

    # convert the list of unique reads into a read-level dataframe
    read_df = pd.DataFrame(read_dict).astype({'well_id': 'category',
                                              'bc1_corr': 'category',
                                              'bc2_corr': 'category',
                                              'seq_corr': 'category',
                                              'edit_type': 'category',
                                              'condition': 'category'})

    return read_df
