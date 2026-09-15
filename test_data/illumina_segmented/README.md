# Rift Valley Fever Segmented Illumina Test Dataset

Reads from: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR25168471&display=metadata

Abstract:

> To manage zoonotic emerging infectious diseases (EIDs), we need to understandtheir origins, their diversity, importance in different communities, and their drivers. Currentgenomic technologies demonstrated great potential and suitability for clinical diagnostics or thereal-time detection and surveillance of epidemics. Here we present a simple, low-cost viralmetatranscriptomic protocol adapted to the smallest Illumina sequencer (iSeq100) capable todetect any DNA and RNA virus on very low concentration using low-throughput sequencing.Despite the comparably much lower throughput of the iSeq100 plastfom relative to the usualIllumina counterparts, the mNGS protocol showed high sensitivity and specificity when appliedto both diagnostic and pathogen discovery in different sample types. This configuration may beimplemented in low-settings or low budget laboratories. The metatranscriptomic workflow isalso suitable as portable (on-site) diagnostic sequencing platform capable to deliver criticalinformation in a clinically relevant turnaround time in case of Disease X and valuable insight intoepidemic transmission and pathogen evolution.

Remarks:

The sample was also downsampled to achieve the 25 reads variation inluded in the tests. Downsampling was completed using:

```
seqtk sample -s100 <name> 25 > <output_name>
```

Reference is from:

- NC_014395.1
- NC_014396.1
- NC_014397.1
