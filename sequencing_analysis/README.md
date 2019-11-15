# Usage
The file `fastq_utils.py` contains functions for analyzing raw FASTQ files from Illumina MiSeq or MiniSeq, given a known set of "template" sequences -- known file sequences that are expected to be observed. The function `analyze_fastq()` in `fastq_utils.py` performs clustering of reads and attempts to assign a cluster to each expected file. A complete explanation of the clustering and assignment algorithm is given in the SI of the accompanying manuscript.

Example FASTQ and template sequence files are included. For example, the following command may be run within an interactive Python session or from a Python script:

`
import fastq_utils

fastq_utils.analyze_fastq(
    path_forward = 'reads_forward.fastq',
    path_reverse = 'reads_reverse.fastq',
    path_templates = 'templates.txt',
    output_prefix = 'TEST',
    sample_id = 'TEST'
)
`

