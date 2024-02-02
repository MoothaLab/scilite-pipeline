# SCI-lite pipeline

## Introduction
This repository contais source code for the pipeline in the following paper: Anna Kotrys , Timothy Durham , Xiaoyan Guo , Venkata Vantaku , Sareh Parangi, and Vamsi Mootha. Large-scale, single cell analysis of mtDNA heteroplasmy dynamics reveals key role of environment-dependent selection.

## Prerequisites
The code runs on Python 3.9.15
- pysam=0.20.0 
- numpy 
- pandas=1.5.2 
- matplotlib

## Quick Start
Download the repository:

`git clone https://github.com/MoothaLab/scilite-pipeline.git scilite_pipeline`

**(optional)** create a new folder for the SCI-lite experiment (eg. `20230811_lineage_tracing`) that you are trying to analyze, and place the argument spreadsheet as well as (encouraged) other required files (such as FASTQ files, barcode mapping files, etc.) in there.
Please see this [tutorial](https://docs.google.com/presentation/d/1nqRbsgIkjIsaFwHQFOxVYje8Xda_JCNxyz9kAJ9B6Ow/edit?usp=sharing) for specifications about the required inputs.

In the following quick start guide, we assume that your folders were structured like this:
```
current working directory
├── scilite_pipeline
│   ├── arg_parse.py
│   ├── calc_heteroplasmy.py
│   ├── driver.py
│   ├── gen_read_info.py
│   ├── utils.py
│   ├── environment.yml
├── 20230811_lineage_tracing
│   ├── exp_args_20230811.csv
│   ├── fastq
│   │   ├── A1_S1_L001_R1_001.fastq.gz
│   │   ├── A1_S1_L001_R2_001.fastq.gz
│   │   ├── ......
│   ├── read1_map.csv
│   ├── barcode1_map.csv
│   ├── barcode2_list.csv
│   ├── Cellecta-SEQ-CloneTracker-XP-Barcode-Libraries-100_x_BC14_list.csv
└── ├── Cellecta-SEQ-CloneTracker-XP-Barcode-Libraries-100K_x_BC30_list.csv
```

The `environment.yml` file provides a quick way of creating a new anaconda environment and installing all dependecies required:

```
conda env create --name scilite --file=scilite_pipeline/environment.yml
conda activate scilite
```

Please remember to save the argument spreadsheet as a .csv file, instead of an .xlsx file!

The pipeline can then be executed by the following command:

```
python ./scilite_pipeline/driver.py -e ./20230811_lineage_tracing/test_args_20230811.csv -p -m 4000 --phase_range 8
```

By default, the pipeline will output a spreadsheet containing single-cell heteroplasmy information `heteroplasmy_dataframe.csv` and a visualization of cell-level UMI counts vs cell rank by UMI counts `knee_plot.pdf` in your experiment folder `20230811_lineage_tracing`.

You can also specify the files to be output to a different location by the `-o [output_location]` option.

The output files can be compressed into `.gz` files if the `-c` option was specified.

Intermediate information aggregated at the read-level (`read_dataframe.csv`) and the UMI-level (`umi_dataframe.csv`) can be additionally included in the output by specifying the `-s` option.

Similarly, the `-k` option can be specify to prevent the knee plot visualization from being included in the output.

The `-m [maximum_number_of_valid_cell]` option can be calibrated to set a hard threshold for finding the appropriate threshold (the "knee") in the cell filtering process using a knee plot, regardless whether `-k` was specified or not. The default setting for this argument is to not impose any hard threshold.

The `--phase_range [maximum_length_of_heterogeneity_spacer]` can be changed to indicate the expected length of heterogeneity spacer if it was used in the experiment. All SCI-lite experiments in the 2023 paper should have this argument set to `--phase_range 8` if heterogeneity spacer were used.

The `-p` option was used to specify that the experiment was a paired-read one.

Please see the helpstrings of the pipeline and the [tutorial](https://docs.google.com/presentation/d/1nqRbsgIkjIsaFwHQFOxVYje8Xda_JCNxyz9kAJ9B6Ow/edit?usp=sharing) for more detailed descriptions of all the command arguments and options.

## Default Edit -> Mutation -> Heteroplasmy Type Mapping of LHON constructs in Kotrys et al 2023 paper
![alt text](https://github.com/MoothaLab/scilite-pipeline-dev//blob/main/fig_edit_mut_het_mapping.png?raw=true)

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/). You are free to share or adapt the material for non-commercial purposes.

If you find this work helpful, please cite:
```
@article{kotrys2023nature,
  title = {Large-scale, single cell analysis of mtDNA heteroplasmy dynamics reveals key role of environment-dependent selection},
  author = {Kotrys, Anna V. and
            Durham, Timothy J. and
            Guo, Xiaoyan A. and
            Vantaku, Venkata R. and
            Parangi, Sareh  and
            Vamsi, Mootha K.},
  journal = {Nature},
  elocation-id = {...},
  doi = {...},
  publisher = {Springer Nature},
  URL = {...},
  eprint = {...},
  year = {2023},
}
```
