# Automated Bioformatics Pipeline for Near-Real-Time Public-Health Responses to New HIV Diagnoses

This software is described in further detail in the manuscript:

> M. Howison, et al. 2022. Informing near-real-time public-health responses
> to new HIV diagnoses in a statewide HIV epidemic

Through a unique academic-public-health partnership, we developed an automated,
open-source pipeline and applied it to prospective, routine analysis of statewide
molecular HIV data in near-real-time in the State of Rhode Island. This collaboration
informed public-health actions to optimize disruption of HIV transmission. This
repository provides the source code for the automated bioinformatics pipeline,
for other public health teams to use and learn from.

## Installation

This pipeline is designed to run in a Linux-based compute cluster environment.
It has been tested on Red Hat Enterprise Linux Server release 7.9 and with the
[SLURM](https://www.schedmd.com/) batch system. (It may be possible to adapt the
pipeline to run with a different batch system by modifying the call to the
SLURM `srun` command in `template-qc/scons/setup` and `template-phylo/scons/setup`
but we have not tested any other batch systems at this time.)

The pipline uses the [scons](https://scons.org/) build system to orchestrate
analyses and provide checkpointing and restarting capabilities.

Dependencies are managed with [Anaconda Python](https://www.anaconda.com/), which
you must install first (you can use the minimal `miniconda` installer if you prefer).
After installing conda, create a new environment called `rtp` from the included
package list with:

    conda create --name rtp --channel kantorlab --file conda.yaml

Once this conda environment is ready, additional R package dependencies must be
installed with:

    conda activate rtp
    Rscript -e 'install.packages(read.table("R-packages.txt")$V1, repos = c("https://cran.r-project.org"))'

Before you run the pipeline, you must first define the following environment variables
for directories where the analysis results (`RTP_ANALYSES`), cached/intermediate results
(`RTP_CACHE`), and input datasets (`RTP_DATASETS`) reside. For example, you could
use subdirectories in your current directory with:

    export RTP_ANALYSES=$PWD/analyses
    export RTP_CACHE=$PWD/cache
    export RTP_DATASETS=$PWD/datasets

## Input data

The pipeline is designed to operate on a **dataset** that you place in the `RTP_DATASETS`
directory. Typically, the dataset would be updated or generated on a near-real-time schedule,
such as every month or week. The dataset must be named using the schema `D{n}_{YYYYMMDD}_V{m}`
where `n` is an accession number starting at 1, `YYYYMMDD` is the most recent sequencing date
for sequences in the dataset, and `m` is a versioning number starting at 1. For a typical
run of the pipeline, the `V1` version will contain **all available sequences for all indivdiuals**
and is used for quality control. The `V2` will contain **one sequence (the oldest available sequence)
per individual** and is used for phylogenetic analysis. For example, say this is your first run
of the pipeline, the latest sequencing date of your data is 2022-10-11, and you plan to run the
quality control and phylogenetic modules. You would create the following directories:

    mkdir $RTP_DATASETS/D1_20221011_V1
    mkdir $RTP_DATASETS/D1_20221011_V2

Then you would place the following required CSV files in those directories. Each file is either
keyed on `StudyID` (the identifier for a distinct individual) or `SequenceID` (the identifier
for a distinct sequence). The `SequenceID` **must have the format** `{StudyID}_{n}` where
`StudyID` is the individual to whom the sequence belongs and `n` is a sequence acession number
starting at 1.

**patients.csv** - a table with one row per StudyID and fields:
    * `StudyID`
    * `Gender` (coded with 1=Male, 2=Female, 3=Transgender)
    * `Ethnicity` (coded with 1=Hispanic, 2=Non-Hispanic)
    * `Race` (coded with 1=White, 2=Black, 3=Asian, 4=Other)
    * `HIVDxDate` (date of HIV diagnosis; *MM/DD/YYYY* format)
    * `HIVDxYear` (year of HIV diagnosis; *YYYY* format)
    * `PrimaryRiskFactor` (a comma-separated list of risk factor labels for the individual)
    * `MSM` (0/1 indicator if MSM is a risk factor for the individual)
    * `EverAtACI` (0/1 indicator if the individual has been previously incarcerated)
    * `EverSubstanceUse` (0/1 indicator if the individual has a known history of controlled substance use)
    * `EverPsychotic` (0/1 indicator if the individual has a previous diagnosis of a mental health condition)
    * `EverIDU` (0/1 indicator if the individual has a known history of intraveneous drug use)
    * `YearOfLastNegativeTest` (year of last known HIV negative test result; *YYYY* format)
    * `CountryOfBirth` (country name or abbreviation)
    * `AgeAtDx` (age at the date of HIV diagnosis)
    * `HIVDx6mo` (0/1 indicator if the HIV diagnosis date was in the previous 6 months)
    * `HIVDx12mo` (0/1 indicator if the HIV diagnosis date was in the previous 12 months)
    * `HIVDx18mo` (0/1 indicator if the HIV diagnosis date was in the previous 18 months)

**sequences.csv** - a table with one row per SequenceID (including **all sequences for quality
control analyses** but **only the earliest available sequence per StudyID for phylogenetic
analyses**) and fields:
    * `SequenceID`
    * `StudyID`
    * `AgeAtSeq`
    * `Year` (year of sequencing; *YYYY* format)
    * `Date` (date of sequencing; *YYYY-MM-DD* format)
    * `TreatmentStatus`
    * `Length` (nucleotide length of the sequence)

**new_seq_ids.csv** - the SequenceID values from the **sequences.csv** file that are newly acquired
since the previuos analysis, with fields:
    * `StudyID`
    * `SequenceID`
    * `Date` (date of sequencing; *YYYY-MM-DD* format)

**sequences.fa** - a FASTA file containing the nucleotide sequences listed in **sequences.csv**.
The FASTA header text for each sequence **must be equal to the** `SequenceID`.

**sierra.json** - the result of processing **sequences.fa** with the
(sierrapy)[https://pypi.org/project/sierrapy/] client for the Stanford HIVdb Sierra algorithm.

**cd4.csv** - a long table (multiple observations per StudyID) of clinical CD4 results with fields:
    * `StudyID`
    * `Date` (date of clinical result; *MM/DD/YYYY* format)
    * `CD4` (clinical value)

**pvl.csv** - a long table (multiple observations per StudyID) of clinical viral load results with fields:
    * `StudyID`
    * `Date` (date of clinical result; *MM/DD/YYYY* format)
    * `PVL` (clinical value)

## Running analyses

To run an analysis, you make a copy of the appropriate template directory for the
module you would like to run (quality control or phylogenetic analysis). An analysis
directory **must be named with the format** `analysis-{n}` where `n` is an analysis
accession number starting at 1. For example, if your first analysis is to run quality
control on the V1 dataset described above, you would copy the `template-qc` like this:

    cp -R template-qc analysis-1

After running the quality control analysis, you could setup a phylogenetic analysis
on the V2 dataset described above like this:

    cp -R template-phylo analysis-2

### Running the quality control (QC) module

The quality control module analyzes the output of the Sierra algorithm and performs
a pairwise distance analysis of sequences.

After making a copy of the QC template above, change to that new analysis directory
with:

    cd analysis-1

The `SConstruct` file in this directory describes the datasets and components used
in the analysis. Edit the `SConstruct` file and fill in the name of your dataset,
`D1_20221011_V1`, on line 7 where the template has `INSERT_DATASET_NAME_HERE`. Fill
in the integer value for the starting `SequenceID` (explained above in the *Input data*
section) where the template has `INSERT_START_SEQUENCE_ID_HERE`. Save and close the
`SConstruct` file.

The entire analysis can be run by executing the `./build` wrapper from this directory.
Each command in the analysis will be automatically run as a SLURM batch job by the
scons build system. You can test the commands that would be run by scons without
actually running them using:

    ./build -n

You can parallelize the build to use concurrent SLURM batch jobs with the `-j`
argument. For example, to use up to 4 concurrent batch jobs, run the analysis
with:

    ./build -j4

When scons finishes running all of the tasks in the `SConstruct` file, it will print
the message:

    scons: done building targets.

You can now find the quality control results in the directory `$RTP_ANALYSES/1/reports/qc`.
There are four reports you can review:

**D1_20221110_V1.Apobec.html** - the number of Apobec mutations identified by the Sierra
algorithm in each sequence.

**D1_20221110_V1.StopCodon.html** - the number of stop codons identified by the Sierra
algorithm in each sequence.

**D1_20221110_V1.Unusual.html** - the number of unusual mutations identified by the Sierra
algorithm in each sequence. These are defined by the Stanford HIVdb as "amino acids with an
overall group M prevalence <0.01% that do not have a mutation penalty score in the HIVdb
genotypic resistance interpretation program and that are not signature APOBEC mutations"
(https://hivdb.stanford.edu/page/release-notes/). 

**mafft.D1_20221110_V1.SameSubtypeLow.html** - sequence pairs that have the same subtype
ordered by the smallest pairwise genetic distance. Sequence pairs with genetic distance
less than 0.5% are very similar and may require further evaluation to confirm that they are
from distinct individuals (e.g. that a StudyID has not been misassigned to one of the
sequences).

### Running the clustering/phylogenetic (phylo) module

The clustering/phylogenetic module performs phylogenetic analysis using the following
five commonly-used phylogenetic methods: RAxML, IQ-TREE, FastTree, FastTree (with
the alternative likelihood ratio test), and MEGA. Phylogenetic results are clustered
using ClusterPicker. The pipeline also performs distance-only cluster analysis using
HIV-TRACE. 

After making a copy of the phylo template above, change to that new analysis directory
with:

    cd analysis-2

The `SConstruct` file in this directory describes the datasets and components used
in the analysis. Edit the `SConstruct` file and fill in the name of your dataset,
`D1_20221011_V2`, on line 7 where the template has `INSERT_DATASET_NAME_HERE`. Fill
in the integer value for the starting `SequenceID` (explained above in the *Input data*
section) where the template has `INSERT_START_SEQUENCE_ID_HERE`. Save and close the
`SConstruct` file.

The entire analysis can be run by executing the `./build` wrapper from this directory.
Each command in the analysis will be automatically run as a SLURM batch job by the
scons build system. You can test the commands that would be run by scons without
actually running them using:

    ./build -n

You can parallelize the build to use concurrent SLURM batch jobs with the `-j`
argument. For example, to use up to 4 concurrent batch jobs, run the analysis
with:

    ./build -j4

When scons finishes running all of the tasks in the `SConstruct` file, it will print
the message:

    scons: done building targets.

You can now find the clustering/phylogenetic results in the directory
`$RTP_ANALYSES/2/reports/summary`, which will contain a PDF report
`D1_20221011_V2_Report.pdf`.

## License

Copyright 2018, Brown University, Providence, RI. All Rights Reserved.

Copyright 2019-2020, Innovative Policy Lab d.b.a. Research Improving People's Lives ("RIPL"), Providence, RI. All Rights Reserved.

See [LICENSE.txt](https://github.com/kantorlab/hiv-real-time-phylogeny/blob/main/LICENSE.txt) for full terms of use.
