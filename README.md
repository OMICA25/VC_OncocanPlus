##VC_OncocanPlus
VC_OncoCanPlus is a Nextflow-based pipeline designed to detect germline and somatic variants from dog (CanFam3.1) targeted sequencing data.
The workflow performs:

Read QC (FASTQC)
BWA-MEM alignment
Read Group assignment (Picard)
On‑target filtering
Coverage computation
Variant calling using Mutect2, VarScan, and FreeBayes
Consensus variant extraction (bcftools isec)
ANNOVAR annotation
Gene symbol enrichment via BioMart
MultiQC summary report

This pipeline is compatible with both local EC2 runs and Omica Platform execution.
