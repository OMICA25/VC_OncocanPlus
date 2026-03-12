#!/bin/bash -ue
set -euo pipefail

# Use dedicated temp directory
export TMPDIR=${TMPDIR:-$PWD/tmp}
mkdir -p "$TMPDIR"

# Align, convert, sort
bwa mem s3://omica-pipeline-data/examples/CanFam3.1ref.fa Run4.IonCode_1517_BL27.fastq     | samtools view -b -     | samtools sort          -@ 3          -m 1G          -T "$TMPDIR/Run4.IonCode_1517_BL27.sort"          -o Run4.IonCode_1517_BL27.sorted.bam

# Index
samtools index -@ 3 Run4.IonCode_1517_BL27.sorted.bam
