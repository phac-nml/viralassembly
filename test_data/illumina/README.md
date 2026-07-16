# SARS-CoV-2 Illumina Amplicon Test Dataset

Reads from: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR39441883&display=metadata

Remarks:

The sample was downsampled to achieve the 50 thousand reads and 25 reads variation inluded in the tests. Downsampling was completed using:

```
seqtk sample -s100 <name> <number of reads> > <output_name>
```

Reference is:

- MN908947.3

# SARS-CoV-2 Illumina WGS Test Dataset

Reads from: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR23802888&display=metadata

Abstract:

> Rapid identification of the rise and spread of severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) variants of concern currently remains critical for monitoring of the efficacy of diagnostics, therapeutics, vaccines, and control strategies. A wide range of SARS-CoV-2 next-generation sequencing (NGS) methods have been developed over the last years, but cross-sequence technology benchmarking studies are scarce. In the current study, 26 clinical samples were sequenced using five protocols: AmpliSeq SARS-CoV-2 (Illumina), EasySeq RC-PCR SARS-CoV-2 (Illumina/NimaGen), Ion AmpliSeq SARS-CoV-2 (Thermo Fisher), custom primer sets (Oxford Nanopore), and capture probe-based viral metagenomics (Roche/Illumina). Studied parameters included genome coverage, depth of coverage, amplicon distribution, and variant calling.The median SARS-CoV-2 genome coverage of samples with cycle threshold (Ct) values of 30 and lower ranged from 81.6 to 99.8 for, respectively, the Oxford Nanopore protocol and Illumina Ampliseq protocol. Correlation of coverage with PCR Ct-values varied and was dependent on the protocol. Amplicon distribution signatures differed across the methods, with peak differences of up to 4 log10 at disbalanced positions in samples with high viral loads (Ct-values <= 23). Phylogenetic analyses of consensus sequences to some extent showed clustering dependent on the workflow, illustrating the limitations of cluster detection when combining results of different platform technologies. The proportion of SARS-CoV-2 reads in relation to background sequences, as a (cost-)efficiency metric, was highest for the EasySeq protocol. The hands-on time was lowest when using EasySeq and ONT protocols, with the latter additionally having the shortest sequence runtime.In conclusion, the studied protocols differed on a variety of the studied metrics. This study provides data that can assist laboratories when selecting protocols for their specific setting.

Remarks:

The sample was downsampled to achieve the 50 thousand reads and 25 reads variation inluded in the tests. Downsampling was completed using:

```
seqtk sample -s100 <name> <number of reads> > <output_name>
```

Reference is:

- MN908947.3
