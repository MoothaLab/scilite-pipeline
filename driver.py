import gc
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# SCI-LITE pipeline module imports
import arg_parse
import calc_heteroplasmy
import gen_read_info


def calc_heteroplasmy_from_fastq(fq_process_kwargs, het_calc_kwargs,
                                 mut_type_map, het_map):
    ''' Main driver function for running the pipeline
    
    Args:
        - fq_process_kwargs: dict, keyword arguments specific for the experiment setup
                                see below for its structure:
                                    {'edit_type_map': edit_type_map (dict {str -> str}, maps target sequence to edit types),
                                     'condition_map': condition_map (dict {str -> str}, maps barcode 1s to conditions),
                                     'construct_info': construct_info (dict {str -> dict}, contains construct-specific info -- see extract_id_info()
                                                                       for detail information),
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
        - het_calc_kwargs: dict, keyword arguments specific for the heterplasmy calculation
                            see below for its structure:
                                { 'required_fields': required_fields (list of strs, the information needed for a read to be considered identifiable),
                                  'max_cells': max_cells (int, the maximum number of cells to be included during knee plot filtering),
                                  'no_plot_knee': no_plot_knee (bool, whether the knee plot should be visualized)
                                }
        - mut_type_map: dict {str -> list of strs}, contains the parsed mutant_type -> edit_type pairs
        - het_map: dict {str -> list of strs}, contains the parsed heteroplasmy_type -> mutant_type pairs

    Returns:
        - results: list of pandas DataFrame(s), all resulted DataFrame(s) from the current pipeline run
    '''
    if fq_process_kwargs['read_dataframe_path'] != None:
        print('reading processed read-level information from '+fq_process_kwargs['read_dataframe_path']+'.......\n')
        read_df = pd.read_csv(fq_process_kwargs['read_dataframe_path'])
    else:
    # process each read and construct an read-level dataframe
        fq_fps = gen_read_info.find_fastq_paths(fq_process_kwargs['fastq_directory_path'], 
                                                fq_process_kwargs['fastq_suffix'])
        read_df = gen_read_info.parse_all_fastq(fq_fps, **fq_process_kwargs)
        print('\nAll FASTQ files have been processed, resulting in {:d} total reads.'.format(read_df.read_count.sum()))

    read_type_stat = read_df.groupby('edit_type').read_count.sum().reset_index()
    for _, row in read_type_stat.iterrows():
        print('{:d} reads were identified as edit type {}'.format(row['read_count'], row['edit_type']))
    print()

    if het_calc_kwargs['if_compress_out'] == True:
        compression = {'method': 'gzip'}
        out_suffix = '.gz'
    else:
        compression = 'infer'
        out_suffix = ''

    read_type_stat.to_csv(het_calc_kwargs['out_fpath'] + '/read_type_stats.csv' + out_suffix, 
                          index=False, compression=compression)

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    axes.bar(read_type_stat['edit_type'], np.log2(read_type_stat['read_count']))
    axes.set_xticklabels(axes.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    axes.set_xlabel('')
    axes.set_ylabel('log2(read count)')
    fig.tight_layout()
    plt.savefig(het_calc_kwargs['out_fpath'] + '/read_type_stats.png', dpi=300)
    plt.close(fig)

    read_cond_stat = read_df.groupby('condition').read_count.sum().reset_index()
    for _, row in read_cond_stat.iterrows():
        print('{:d} reads were identified as condition {}'.format(row['read_count'], row['condition']))
    print()
    read_cond_stat.to_csv(het_calc_kwargs['out_fpath'] + '/read_condition_stats.csv' + out_suffix, 
                          index=False, compression=compression)

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    axes.bar(read_cond_stat['condition'], np.log2(read_cond_stat['read_count']))
    axes.set_xticklabels(axes.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    axes.set_xlabel('')
    axes.set_ylabel('log2(read count)')
    fig.tight_layout()
    plt.savefig(het_calc_kwargs['out_fpath'] + '/read_condition_stats.png', dpi=300)
    plt.close(fig)

    # filter out reads with any unidentifiable necessary field
    filtered_df = calc_heteroplasmy.filter_unknown_reads(read_df.copy(), 
                                                         het_calc_kwargs['required_fields'], 
                                                         fq_process_kwargs['tracer_edit_types'])

    if fq_process_kwargs['if_save_intermediates'] == True:
        if len(fq_process_kwargs['tracer_edit_types']) > 0 and 'bc_tracer' not in read_df.columns:
            read_df.to_csv(het_calc_kwargs['out_fpath'] + '/read_dataframe.csv', 
                           index=False, compression=compression)
        else:
            if 'bc_tracer' in read_df.columns:
                # read_df.drop(columns='bc_tracer').to_csv(het_calc_kwargs['out_fpath'] + '/read_dataframe.csv' + out_suffix,
                #                                          index=False, compression=compression)
                read_df.to_csv(het_calc_kwargs['out_fpath'] + '/read_dataframe.csv' + out_suffix,
                               index=False, compression=compression)

    del read_df
    gc.collect()
    print('{:d} identifiable reads passed filtering.'.format(filtered_df.read_count.sum()))
    print()


    # collapse the read-level dataframe into one that contains UMI-level information
    # including deciding each UMI's mutation type by majority voting
    umi_df = calc_heteroplasmy.collapse_by_umi(filtered_df.copy(), 
                                               fq_process_kwargs['tracer_edit_types'])
    del filtered_df
    gc.collect()
    print('{:d} UMIs with at least 3 reads mapped to them passed filtering.'.format(umi_df.shape[0]))
    print()

    # map each read's edit type to its corresponding mutation
    # since each mutation can take the form of multiple edits
    umi_df = calc_heteroplasmy.map_edit_to_mut_type(umi_df, mut_type_map)

    # collapse the UMI-level dataframe into one that contains cell-level information
    cell_df = calc_heteroplasmy.collapse_by_cell(umi_df, fq_process_kwargs['tracer_edit_types'])
    if fq_process_kwargs['if_save_intermediates'] == True:
        if len(fq_process_kwargs['tracer_edit_types']) > 0 and 'bc_tracer' not in umi_df.columns:
            umi_df.to_csv(het_calc_kwargs['out_fpath'] + '/umi_dataframe.csv', 
                          index=False, compression=compression)
        else:
            if 'bc_tracer' in umi_df.columns:
                umi_df.drop(columns='bc_tracer').to_csv(het_calc_kwargs['out_fpath'] + '/umi_dataframe.csv' + out_suffix, index=False)
        cell_df.to_csv(het_calc_kwargs['out_fpath'] + '/cell_dataframe.csv', 
              index=False, compression=compression)
    del umi_df
    gc.collect()
    print('{:d} unique cells have been identified from UMI collapsing.'.format(cell_df.shape[0]))
    print()

    # filter out cells with too few UMIs mapped to them using a knee plot
    filtered_cell_df = calc_heteroplasmy.filter_cells_by_knee_plot(cell_df, 
                                                                   het_calc_kwargs['max_cells'], 
                                                                   het_calc_kwargs['mut_types_for_filtering'],
                                                                   no_plot_knee=het_calc_kwargs['no_plot_knee'],
                                                                   out_fpath=het_calc_kwargs['out_fpath'])
    del cell_df
    gc.collect()

    
    het_df = calc_heteroplasmy.calc_heteroplasmy(filtered_cell_df, het_map, fq_process_kwargs['tracer_edit_types'])
    del filtered_cell_df
    het_df.to_csv(het_calc_kwargs['out_fpath'] + '/heteroplasmy_dataframe.csv' + out_suffix, 
                  index=False, compression=compression)
    
    return het_df

def process_args(cl_args, exp_args):
    ''' Process all arguments and assign them into different key word argument structures
        to be used in downstream data analysis
        
    Args:
        - cl_args: dict {str -> *}, contains arguments filled by command line inputs
        - exp_args: dict {str -> *}, contains experiment set-up related arguments filled
                                     by file inputs
                                     
    Returns:
        - fq_process_kwargs: dict, keyword arguments specific for the experiment setup
        - het_calc_kwargs: dict, keyword arguments specific for the heterplasmy calculation
        - mut_type_map: dict {str -> list of strs}, contains the parsed mutant_type -> edit_type pairs
        - het_map: dict {str -> list of strs}, contains the parsed heteroplasmy_type -> mutant_type pairs
    '''

    fq_process_kwargs = arg_parse.process_barcode_and_seq_mappings({}, exp_args)
    fq_process_kwargs['construct_info'] = exp_args['construct_info']
    fq_process_kwargs['tracer_edit_types'] = exp_args['tracer_edit_types']

    fq_process_kwargs['if_paired_end'] = cl_args.isPairedEnd
    fq_process_kwargs['if_well_id_from_read'] = cl_args.isWellIDFromRead
    fq_process_kwargs['if_reverse_complement_id_seq'] = cl_args.reverseComplementIDs
    fq_process_kwargs['if_save_intermediates'] = cl_args.saveIntermediates
    fq_process_kwargs['fastq_directory_path'] = exp_args['fastq_directory_path']
    fq_process_kwargs['fastq_suffix'] = exp_args['fastq_suffix']
    fq_process_kwargs['phase_range'] = cl_args.phase_range
    fq_process_kwargs['read_dataframe_path'] = cl_args.read_dataframe

    het_calc_kwargs = {}
    if cl_args.maxCells:
        het_calc_kwargs['max_cells'] = int(cl_args.maxCells)
    else:
        het_calc_kwargs['max_cells'] = None
    het_calc_kwargs['no_plot_knee'] = cl_args.noPlotKnee
    het_calc_kwargs['if_compress_out'] = cl_args.compressOutputFiles
    het_calc_kwargs['required_fields'] = exp_args['required_fields']
    het_calc_kwargs['mut_types_for_filtering'] = exp_args['mut_types_for_filtering']
    if cl_args.outputPath is None:
        cl_args.outputPath = '/'.join(cl_args.experimentSettingFile.split('/')[:-1])
    elif os.path.isdir(cl_args.outputPath) == False:
        os.mkdir(cl_args.outputPath)
    het_calc_kwargs['out_fpath'] = cl_args.outputPath

    mut_type_map = exp_args['mut_type_map']
    het_map = exp_args['het_map']

    return fq_process_kwargs, het_calc_kwargs, mut_type_map, het_map

if __name__ == "__main__":
    ''' Main driver function for the pipeline '''
    parser = arg_parse.create_parser()
    cl_args = parser.parse_args()

    exp_args = arg_parse.parse_input_arg_file(cl_args.experimentSettingFile)

    print('Parsing input files......')

    fq_process_kwargs, het_calc_kwargs, mut_type_map, het_map = process_args(cl_args, exp_args)

    if het_calc_kwargs['if_compress_out'] == True:
        print(" -c option is on, output files will be compressed with gzip...this may increase execution time.")
        print()

    results = calc_heteroplasmy_from_fastq(fq_process_kwargs, het_calc_kwargs,
                                           mut_type_map, het_map)
