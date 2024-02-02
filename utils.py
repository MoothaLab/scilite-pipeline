import glob
import itertools
import numpy as np
import os

def gen_reverse_complement(dna_seq):
    ''' Generate the reverse complement of a DNA sequence
    
    Args:
        dna_seq: str, the DNA sequence to be reversely complemented with
                 and contains only ATCG bases
    Returns:
        rc_seq: str, the reverse complement of the input DNA sequence
    '''
    bp_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    rc_seq = ''.join([bp_map[elt] for elt in list(dna_seq[::-1])])

    return rc_seq

def gen_sub_map(seq_list, n_max_subs=2, immutable_pos=[]):
    ''' Generate a substitution-tolerated DNA sequence dictionary that 
        maps the substituted DNA sequence to the original one
    
    Args:
        - seq_list: list of strs, the list containing all original DNA sequences (A, T, C, G, N)
                    all sequences should have the same length
        - n_max_subs: int, the maximum number of substitution tolerated
        - immutable_pos: list of ints, the list containing all immutable 
                         sequence base position indices (0-indexed)
    Returns:
        - seq_map: dict {str -> str}, maps the substituted DNA sequence (key) 
                  to the original one (value)
    '''

    immutable_pos = set(immutable_pos)
    seq_map = {}
    seq_len = len(seq_list[0])

    # input sanity check
    if seq_len < n_max_subs:
        raise Exception('Number of substitutions allowed ({!s}) must be less '
                        'than or equal to the sequence length ({!s})'.format(n_max_subs, seq_len))

    # generate all possible index combinations for n-base substitutions
    edit_poss = list(itertools.combinations(np.arange(seq_len), # 0, 1, 2, ..., n_max_subs
                                            n_max_subs))
    # generate all possible A, T, C, G, N base permutations for n-base substitutions
    edit_bases = list(set(itertools.permutations(['A', 'C', 'G', 'T', 'N'], n_max_subs)))

    # iterate through all sequences
    for seq in seq_list:
        seq_base_list = list(seq)

        # iterate through all possible (i-th index, j-th base) pairs generated
        for edit_pos, edit_base in itertools.product(edit_poss, edit_bases):

            # create a copy of the original sequence
            seq_tmp = list(seq_base_list)

            # mutates the sequence
            # and maps the mutated sequence to the original sequence
            # omit the base mutation if it occurred at any immutable base
            for idx, base in zip(edit_pos, edit_base):
                if idx in immutable_pos:
                    continue
                try:
                    seq_tmp[idx] = base
                except:
                    raise Exception(seq_tmp, edit_pos, edit_base, idx, base)
            seq_map[''.join(seq_tmp)] = seq

    return seq_map
