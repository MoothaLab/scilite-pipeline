import gc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings(action='ignore')

def filter_tracer_fail_reads(unfiltered_df, tracer_edit_types):
    ''' Filter out read-level dataframe entries that belong to any tracer construct
        but does not have an identifiable tracer barcode
    
    Args:
        - unfiltered_df:  pandas DataFrame, the read-level dataframe to be filtered (note the filtering occurs in-place)
        - tracer_edit_types: list of strs, contains all edit types belonging to any tracer construct
    
    Returns:
        - filtered_df: pandas DataFrame, the filtered read-level dataframe
    '''

    # perform the filtering for each tracer construct, and do it in place on unfiltered_df
    for tracer_type in tracer_edit_types:
        unfiltered_df['tracer_fail'] = False
        tracer_reads = unfiltered_df[unfiltered_df['edit_type']==tracer_type].copy()
        unfiltered_df.loc[tracer_reads[tracer_reads.bc_tracer.str.contains('UNKNOWN')].index, 'tracer_fail'] = True
        unfiltered_df = unfiltered_df[~unfiltered_df.tracer_fail]
    filtered_df = unfiltered_df

    return filtered_df

def filter_unknown_reads(unfiltered_df, required_fields, tracer_edit_types):
    ''' Filter out read-level dataframe entries that contains UNKNOWN values in any of the required fields
    
    Args:
        - unfiltered_df:  pandas DataFrame, the read-level dataframe to be filtered (note the filtering occurs in-place)
        - required_fields: list of strs, contains the column names for the required fields
    
    Returns:
        - filtered_df: pandas DataFrame, the filtered DataFrame
    '''

    if len(tracer_edit_types) != 0:
        print('filtering out tracer fails......')
        unfiltered_df = filter_tracer_fail_reads(unfiltered_df, tracer_edit_types)
    print('filtering out unidentifiable reads......')
    criteria = (unfiltered_df.loc[:, required_fields]!='UNKNOWN').all(axis='columns')
    filtered_df = unfiltered_df[criteria]

    return filtered_df

def map_edit_to_mut_type(umi_df, mut_type_map):
    ''' Maps the edit type to its specified corresponding mutation type
    
    Args:
        - umi_df: pandas DataFrame, contains each unique UMI molecule as an entry
        - mut_type_map: dict {str -> list of strs}, maps the mutattion type to all edit types 
                                                    that are considered to be belonging to 
                                                    the same mutation type
    
    Returns:
        - umi_df: pandas DataFrame, now appended with a new field "mut_type"
    '''

    umi_df['mut_type'] = umi_df['edit_type'].map(mut_type_map).astype('category')

    return umi_df

def collapse_by_umi(read_df, tracer_edit_types):
    ''' Collapse reads that belong to the same cell and the same molecule and generate UMI-level Dataframe
    
    Args:
        - read_df: pandas DataFrame, contains read-level information including "mut_type"
        - tracer_edit_types: list of strs, contains all edit types belonging to any tracer construct
    
    Returns:
        - umi_df: pandas DataFrame, contains UMI-level information
    '''
    print('filtering out UMIs with fewer than 3 reads......')

    # construct cell-level identifer using well_id and corrected barcodes
    read_df['cell_id'] = (read_df['well_id'].astype('str') + '-' + 
                          read_df['bc1_corr'].astype('str') + '-' + 
                          read_df['bc2_corr'].astype('str')).astype('category')

    # filter out any reads belonging to UMIs with fewer than 3 reads mapped to it 
    read_df.loc[:, 'cell_umi'] = (read_df['cell_id'].astype('str') + '_' + 
                                  read_df['umi'].astype('str')).astype('category')
    well_sup_umis_rc = read_df.groupby('cell_umi')['read_count'].sum().reset_index()
    well_sup_umis = well_sup_umis_rc.query('read_count > 2')['cell_umi']
    read_df = read_df[read_df['cell_umi'].isin(well_sup_umis)]
    if len(tracer_edit_types) != 0:
        read_df['edit_tracer_type'] = (read_df['edit_type'].astype('str') + '_' + 
                                           read_df['bc_tracer'].astype('str')).astype('category').copy()

    # majority voting of allele type for each UMI
    print('Determining UMI types by majority voting......')
    # fields to be included in the UMI-level dataframe
    umi_df_preserved_cols = ['cell_id', 'umi', 
                             'edit_type', 
                             'bc1_corr', 'bc2_corr', 
                             'condition']
    # if tracer constructs are involved in the experiment, each unique tracer barcode
    # should be considered as an independent edit type that requires voting
    if len(tracer_edit_types) != 0:
        umi_df_preserved_cols.extend(['edit_tracer_type', 'bc_tracer'])

    umi_df = (read_df.groupby(umi_df_preserved_cols, observed=True)['read_count']
              .sum()
              .reset_index())
    #count total # reads per UMI
    umi_count = (umi_df.groupby(['cell_id', 'umi'], observed=True)['read_count']
                 .sum()
                 .rename('umi_read_count'))
    umi_df = umi_df.join(umi_count, on=['cell_id', 'umi'])

    #filter UMIs based on the alleles with a large fraction of the reads for that UMI
    umi_df['frac_umi_reads'] = umi_df['read_count']/umi_df['umi_read_count']
    umi_df = umi_df[umi_df['frac_umi_reads'] > 0.66].copy()

    # include the total number of reads mapped to each UMI as well
    umi_df = umi_df.rename(columns={'read_count': 'umi_majority_allele_read_count', 
                                    'umi_read_count': 'umi_total_read_count'})

    return umi_df

def collapse_by_cell(umi_df, tracer_edit_types):
    ''' Collapse UMIs that belong to the same cell and generate cell-level Dataframe
    
    Args:
        - umi_df: pandas DataFrame, contains UMI-level information
        - tracer_edit_types: list of strs, contains all edit types belonging to any tracer construct
    
    Returns:
        - cell_df: pandas DataFrame, contains cell-level information
    '''
    print('collapsing UMIs into unique cells......')

    # group UMIs by their cell id and condition
    grouped = umi_df.groupby(['cell_id', 'condition'], observed=True)

    # generate cell-level total UMI counts
    umi_count = grouped.size().rename('umi_count')

    # generate cell-level UMI counts for each mutation type
    mut_count = grouped[['mut_type']].value_counts().rename('mut_count')
    mut_count_cols = mut_count.reset_index().pivot(index=['cell_id', 'condition'],
                                                   columns='mut_type',
                                                   values='mut_count').fillna(0)
    cell_df = pd.concat([umi_count, mut_count_cols], axis=1, join='inner')

    # if any tracer construct is present in the experiment, we will also need the number of unique tracer barcodes 
    # as well as that of all tracer UMIs present in each cell
    if len(tracer_edit_types) != 0:

        # only take the tracer UMIs for now, and group them by their cell ids and conditions
        grouped = umi_df[umi_df.edit_type.isin(tracer_edit_types)].groupby(['cell_id', 'condition'], observed=True)

        # count the number of all tracer UMIs present in each cell and append to the cell-level dataframe
        tracer_umi_count = grouped.size().rename('tracer_umi_count').to_frame()
        cell_df = cell_df.join(tracer_umi_count['tracer_umi_count'],
                                    on=['cell_id', 'condition'])

        # iterate over each unique tracer construct
        for tracer_type in tracer_edit_types:
            # count number of unique tracer barcodes belonging to the construct being considered
            tracer_umi = umi_df[umi_df.edit_type==tracer_type].copy()
            top_tracers = tracer_umi.groupby(['cell_id', 'condition'], observed=True)['bc_tracer'].value_counts()
            tracer_type_counts = top_tracers.groupby(level=[0,1]).size().rename('tracer_type_count')
            top_tracers = top_tracers.groupby(level=[0, 1], group_keys=False).head(2).rename('top2_tracer_count').reset_index()
            grouped = top_tracers.groupby(['cell_id', 'condition'], observed=True)
            top_tracer_seqs = grouped['bc_tracer'].apply(lambda x: pd.Series(x.values)).unstack()
            top_tracer_seqs = top_tracer_seqs.rename(columns={i: 'top{}_{}_sequence'.format(i+1, tracer_type) for i in range(top_tracer_seqs.shape[1])})
            top_tracer_counts = grouped['top2_tracer_count'].apply(lambda x: pd.Series(x.values)).unstack()
            cols = ["top{}_{}_count".format(i+1, tracer_type) for i in range(top_tracer_counts.shape[1])]
            top_tracer_counts = top_tracer_counts.rename(columns={i: 'top{}_{}_count'.format(i+1, tracer_type) for i in range(top_tracer_counts.shape[1])})

            cell_df = cell_df.join(tracer_type_counts, on=['cell_id', 'condition'])
            cell_df = cell_df.join(top_tracer_seqs[['top1_'+tracer_type+'_sequence', 'top2_'+tracer_type+'_sequence']],
                                    on=['cell_id', 'condition'])
            cell_df = cell_df.join(top_tracer_counts[['top1_'+tracer_type+'_count', 'top2_'+tracer_type+'_count']],
                                    on=['cell_id', 'condition'])

    return cell_df.reset_index()

def plot_knee(cell_df, knee_thresh, out_fpath, reference_cell_num=None):
    ''' Generate the knee plot to visualize the process of filtering cells 
        by their UMI counts
        
    Args:
        - cell_df: pandas DataFrame, contains cell-level information
        - knee_thresh: int, the minimum number of UMIs required for a cell to be considered valid
        - out_fpath: str, where the generated knee plot should be exported
        - reference_cell_num: int, the rank of cell-level UMI count that the knee_thresh corresponds to
    '''

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    axes.scatter(np.log10(cell_df['umi_count_rank']), np.log10(cell_df['umi_count_for_filtering']), s=5)
    if reference_cell_num is not None:
        axes.axvline(np.log10(reference_cell_num), color='darkgrey', linestyle='dotted')
    axes.axvline(np.log10(np.sum(cell_df['umi_count_for_filtering'] > knee_thresh)), color='k', linestyle='dotted')
    axes.axhline(np.log10(knee_thresh), color='red', linestyle='dotted')
    axes.set_xlabel('log10(Cell rank)')
    axes.set_ylabel('log10(UMI count)')
    fig.tight_layout()
    plt.show()
    plt.savefig(out_fpath+'/knee_plot.pdf', dpi=300)

def filter_cells_by_knee_plot(cell_df, max_cells, mut_types_for_filtering,
                              reference_cell_num=None, no_plot_knee=False, out_fpath=None):
    ''' Filter out cells that fall below the "knee" in a knee plot of cell-level UMI counts
    
    Args:
        - cell_df: pandas DataFrame, contains cell-level information
        - max_cells: int, the maximum rank of cells that can be considered as being the valid "knee" 
        - mut_types_for_filtering: list of strs, the mutation types whose UMI counts are to be considered during the filtering
        - reference_cell_num: int, the rank of cell-level UMI count that the knee_thresh corresponds to
        - no_plot_knee: bool, do NOT visualize and export the resulting knee plot if set to TRUE
        - out_fpath: str, where the generated knee plot should be exported
    
    Returns:
        - filtered_cell_df: pandas DataFrame, contains cell-level information belonging to cells that
                            meet the minimum UMI count criteria
    '''
    print('filtering out cells with too few UMIs using a knee plot......')

    # if the user did not provide a customized max_cell value
    # default shall be the size of unfiltered cell-level dataframe
    # i.e. no upper limit
    if max_cells is None:
        max_cells = cell_df.shape[0]

    cell_df['umi_count_for_filtering'] = cell_df[mut_types_for_filtering].sum(axis=1)

    # rank the cells by the count of UMIs mapped to them, in a descending order
    # i.e. cells with the rank of 1 will have the highest UMI count
    cell_df['umi_count_rank'] = cell_df.shape[0] - np.argsort(np.argsort(cell_df['umi_count_for_filtering'].to_numpy()))

    # find the "knee" by representing cells in a 2-dimensional space
    # with the x coordinate being the cell's rank in total UMI counts
    # and the y coordinate being the cell's total UMI counts
    # the "knee" is therefore where the cell that has the highest Euclidean distance from the origin is
    # the upper limit of the x coordinates that will be considered in this "knee"-finding process is 
    # defined as the value of max_cells
    rank_cutoff_mask = (cell_df['umi_count_rank'] <= max_cells)
    y_coords = np.log10(cell_df.loc[rank_cutoff_mask, 'umi_count_for_filtering'].to_numpy())
    x_coords = np.log10(cell_df.loc[rank_cutoff_mask, 'umi_count_rank'].to_numpy())
    origin_dist = np.sqrt(x_coords**2 + y_coords**2)
    knee_thresh = cell_df.loc[rank_cutoff_mask, 'umi_count_for_filtering'].to_numpy() [np.argmax(origin_dist)]

    # visualize and export the representation of cells in the 2-d space described above
    # if the user desires
    if no_plot_knee == False:
        plot_knee(cell_df, knee_thresh, out_fpath)

    # filter out cells that do not meet the minimum UMI count threshold
    # i.e. ranking below the "knee" in the knee plot
    filtered_cell_df = cell_df.loc[cell_df['umi_count_for_filtering'] >= knee_thresh].copy()
    print('{:d} cells with at least {:d} UMIs passed filtering.'.format(filtered_cell_df.shape[0], knee_thresh))

    return filtered_cell_df

def calc_heteroplasmy(cell_df, het_map, tracer_edit_types):
    ''' Calculate the heteroplasmy level of each cell
    
    Args:
        - cell_df: pandas DataFrame, contains cell-level information
        - het_map: dict {str -> list of strs}, contains the list of mutation types being 
                                               considered in each different heteroplasmy type
        - tracer_edit_types: list of strs, contains all edit types belonging to any tracer construct
    
    Returns:
        - cell_df: pandas DataFrame, contains cell-level information, and now appended with additional
                                     columns containing the desired heteroplasmy information
    '''

    # iterate over each heteroplasmy type that needs to be calculated
    for het in het_map.keys():
        mut_types = het_map[het]

        # the numerator includes all mutation types considered to be
        # contributing to the current heteroplasmy type
        numerator = cell_df.loc[:, mut_types]

        # sometimes there is only one mutation type considered to be
        # contributing to the current heteroplasmy type
        # and in this case we will need to skip the summing-over-columns
        # operation to avoid Pandas error
        if len(numerator.shape) == 1:
            numerator = numerator
        else:
            numerator = numerator.sum(axis=1)

        # if the heteroplasmy type is just for an tracer construct
        # then the denominator should just be the total UMI count
        # of the cell
        if het in tracer_edit_types:
            denominator = cell_df['umi_count'].copy()
        # otherwise, the denominator should not really consider UMI counts 
        # contributed by the tracer constructs 
        else:
            denominator = cell_df['umi_count'].copy()
            for tracer_edit_type in tracer_edit_types:
                denominator -= cell_df[tracer_edit_type]

        cell_df[het+'_het'] = numerator / denominator
        cell_df[het+'_het'] = cell_df[het+'_het'].fillna(0).replace([np.inf, -np.inf], 0)

    return cell_df
