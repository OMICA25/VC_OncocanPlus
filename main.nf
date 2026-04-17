/*
VC_OncoCanPlus v1.0 — main.nf
Author: Antonia Noce

Target:
Targeted plasma cfDNA sequencing panel for dog (CanFam3.1),
Ion Torrent single-end reads.

Description:
Nextflow pipeline for somatic variant calling (SNVs and small indels)
from plasma cfDNA targeted sequencing data. The workflow includes
quality control, alignment, on-target filtering, variant calling
(VarScan and Mutect2), variant normalization, annotation with ANNOVAR,
and gene symbol enrichment.

Execution:
Containerized execution using Docker through Nextflow.
 */

nextflow.enable.dsl=2

// ------------------------------
// PARAMETERS
// ------------------------------
params.reads           = params.reads           ?: null
params.fasta           = params.fasta           ?: null
params.bed             = params.bed             ?: null
params.outdir          = params.outdir          ?: "./results"
params.annovar_home    = params.annovar_home    ?: null
params.annovar_db      = params.annovar_db      ?: null
params.isec_n          = params.isec_n          ?: "+2"
params.biomart_host    = params.biomart_host    ?: "https://may2021.archive.ensembl.org"
params.biomart_dataset = params.biomart_dataset ?: "clfamiliaris_gene_ensembl"
params.keep_intermediates = params.keep_intermediates ?: false

// Fail fast on missing critical parameters
assert params.reads, "ERROR: params.reads not set"
assert params.fasta, "ERROR: params.fasta not set"
assert params.bed,   "ERROR: params.bed not set"
assert params.annovar_home, "ERROR: params.annovar_home not set"
assert params.annovar_db,   "ERROR: params.annovar_db not set"


// Reference files (staged as Nextflow file objects)
fasta_file  = file(params.fasta)
fasta_fai   = file(params.fasta + ".fai")
fasta_dict  = file(params.fasta.replaceFirst(/\.[^.]+$/, "") + ".dict")
bed_file    = file(params.bed)
bwa_bwt = file(params.fasta + ".bwt")
bwa_ann = file(params.fasta + ".ann")
bwa_amb = file(params.fasta + ".amb")
bwa_pac = file(params.fasta + ".pac")
bwa_sa  = file(params.fasta + ".sa")

//Annovar files
annovar_bin_dir = file(params.annovar_home)
annovar_db_dir  = file(params.annovar_db)
//file for adding gene_symbol
add_gene_symbols_py = file("s3://omica-pipeline-data/examples/add_gene_symbols.py")
panel_csv = file("s3://omica-pipeline-data/examples/Target_panel_gene_symbol.csv")
// ------------------------------
// CHANNELS (Single-End FASTQ)
// ------------------------------
Channel
    .fromPath(params.reads, checkIfExists:true)
    .map { f ->
        def base = f.getName().replaceFirst(/\.fastq(\.gz)?$/, "")
        tuple([id: base], f)
    }
    .set { READS }

// Reference genome + BWA index (staged by Nextflow)
Channel
  .fromPath(
    "${params.fasta}{,.bwt,.ann,.amb,.pac,.sa}",
    checkIfExists: true
  )
  .collect()
  .set { REF }




// ------------------------------
// FASTQC
// ------------------------------
process FASTQC {
  tag { meta.id }
  input:
    tuple val(meta), path(read)
  output:
    tuple val(meta),
          path("${meta.id}_fastqc.html"),
          path("${meta.id}_fastqc.zip")
  script:
  """
  set -euo pipefail
  fastqc ${read} --outdir .
  """
}

// ------------------------------
// MULTIQC (only FASTQC + flagstat + stats)
// ------------------------------
process MULTIQC {
  tag "multiqc"
  publishDir "${params.outdir}/qc", mode: 'copy'
  input:
    path qc_files
  output:
    path "multiqc_report.html"
  script:
  """
  set -euo pipefail
  multiqc ${qc_files.join(' ')} -o .
  """
}

// ------------------------------
// ALIGNMENT (classic bwa mem)
// -----------------------------
process ALIGNMENT {
  label 'heavy_sort'
  tag { meta.id }
  publishDir "${params.outdir}/bam", mode: 'copy'

  cpus 3
  memory '12 GB'

  input:
    tuple val(meta),
          path(read),
          path(fasta_file),
          path(bwa_bwt),
          path(bwa_ann),
          path(bwa_amb),
          path(bwa_pac),
          path(bwa_sa)


  output:
    tuple val(meta),
          path("${meta.id}.sorted.bam"),
          path("${meta.id}.sorted.bam.bai")

  script:
  """
  set -euo pipefail


  bwa mem ${fasta_file} ${read} \\
     | samtools sort \\
        -@ ${task.cpus} \\
        -m 1G \\
        -o ${meta.id}.sorted.bam

  samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
  """
}
// ------------------------------
// ADD READ GROUPS
// ------------------------------
process ADD_READ_GROUP {
  label 'heavy_sort'
  tag { meta.id }
  publishDir "${params.outdir}/bam/rg_sorted", mode:'copy', overwrite:true

  cpus 3
  memory '8 GB'

  input:
    tuple val(meta), path(bam), path(bai)

  output:
    tuple val(meta),
          path("${meta.id}.RG.sorted.bam"),
          path("${meta.id}.RG.sorted.bam.bai")

  script:
  """
  set -euo pipefail

  # --- Use dedicated large TMPDIR to avoid /tmp exhaustion ---
  export TMPDIR=\${TMPDIR:-\$PWD/tmp}
  mkdir -p "\$TMPDIR"

  # --- Add read groups ---
  picard AddOrReplaceReadGroups \
      I=${bam} \
      O=${meta.id}.RG.bam \
      RGID=${meta.id} \
      RGLB=lib1 \
      RGPL=IONTORRENT \
      RGPU=unit1 \
      RGSM=${meta.id} \
      VALIDATION_STRINGENCY=LENIENT \
      CREATE_INDEX=false

  # --- Memory‑safe sort ---
  samtools sort \
      -@ ${task.cpus} \
      -m 1G \
      -T "\$TMPDIR/${meta.id}.rg" \
      -o ${meta.id}.RG.sorted.bam \
      ${meta.id}.RG.bam

  # --- Index ---
  samtools index -@ ${task.cpus} ${meta.id}.RG.sorted.bam

  # Debug info
  samtools view -H ${meta.id}.RG.sorted.bam | grep '^@RG' || true
  """
}

// ------------------------------
// ON-TARGET BED FILTERING
// ------------------------------
process BAM_FILTER {
  tag { meta.id }
  publishDir "${params.outdir}/bam/filtered", mode:'copy'
  input:
    tuple val(meta), path(bam), path(bai), path(bed_file)
  output:
    tuple val(meta),
          path("${meta.id}.on_target.bam"),
          path("${meta.id}.on_target.bam.bai")
  script:
  """
  set -euo pipefail
  samtools view -b -L ${bed_file} ${bam} > ${meta.id}.on_target.bam
  samtools index ${meta.id}.on_target.bam
  """
}

// ------------------------------
// SAMTOOLS FLAGSTAT
// ------------------------------
process SAMTOOLS_FLAGSTAT {
  tag { meta.id }
  publishDir "${params.outdir}/bam", mode:'copy'
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    tuple val(meta), path("${meta.id}.flagstat.txt")
  script:
  """
  set -euo pipefail
  samtools flagstat ${bam} > ${meta.id}.flagstat.txt
  """
}

// ------------------------------
// SAMTOOLS STATS
// ------------------------------
process SAMTOOLS_STATS {
  tag { meta.id }
  publishDir "${params.outdir}/bam", mode:'copy'
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    tuple val(meta), path("${meta.id}.samtools.stats")
  script:
  """
  set -euo pipefail
  samtools stats ${bam} > ${meta.id}.samtools.stats
  """
}

// ------------------------------
// VARSCAN CALL (reheader + normalize + PASS)
// ------------------------------
process VARSCAN_CALL {

  tag { meta.id }
  publishDir "${params.outdir}/vcf/varscan", mode: 'copy'

  input:
    tuple val(meta),
          path(bam),
          path(bai),
          path(bed),
          path(fasta)

  output:
    tuple val(meta),
          path("${meta.id}.varscan.norm.pass.vcf.gz"),
          path("${meta.id}.varscan.norm.pass.vcf.gz.tbi")

  script:
"""
set -euo pipefail

# 1) Raw calling (BED restricted, VarScan-compatible mpileup)
samtools mpileup \
  -f ${fasta} \
  -l ${bed} \
  -B \
  -q 1 \
  -Q 1 \
  ${bam} \
| varscan mpileup2cns \
    --min-var-freq 0.01 \
    --p-value 0.05 \
    --strand-filter 1 \
    --output-vcf 1 \
> ${meta.id}.varscan.vcf

# 2) Add ##contig via reheader using FASTA .fai
bcftools reheader -f ${fasta}.fai \
    ${meta.id}.varscan.vcf \
    -o ${meta.id}.varscan.reheader.vcf

# 3) Sort, bgzip, and index
bcftools sort -Oz \
    -o ${meta.id}.varscan.sorted.vcf.gz \
    ${meta.id}.varscan.reheader.vcf
tabix -f -p vcf ${meta.id}.varscan.sorted.vcf.gz

# 4) Normalize
bcftools norm -f ${fasta} --check-ref w -m -both \
    -Oz -o ${meta.id}.varscan.norm.vcf.gz \
    ${meta.id}.varscan.sorted.vcf.gz
tabix -f -p vcf ${meta.id}.varscan.norm.vcf.gz

# 5) Keep only PASS variants
bcftools view -f .,PASS -Oz \
    -o ${meta.id}.varscan.norm.pass.vcf.gz \
    ${meta.id}.varscan.norm.vcf.gz
tabix -f -p vcf ${meta.id}.varscan.norm.pass.vcf.gz

# 6) Cleanup intermediates
rm -f \
  ${meta.id}.varscan.vcf \
  ${meta.id}.varscan.reheader.vcf \
  ${meta.id}.varscan.sorted.vcf.gz* \
  ${meta.id}.varscan.norm.vcf.gz* || true
"""
}
// ------------------------------
// MUTECT2 CALL (reheader + normalize + PASS)
// ------------------------------
process MUTECT2_CALL {
    tag { meta.id }
    publishDir "${params.outdir}/vcf/Mutect2", mode: 'copy', overwrite: true
    cpus 4
    memory '6 GB'

    input:
        tuple val(meta), path(bam), path(bai), path(bed), path(fasta), path(fai), path(dict)

    output:
        tuple val(meta),
              path("${meta.id}.mutect2.norm.pass.vcf.gz"),
              path("${meta.id}.mutect2.norm.pass.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    # 1) Call variants (SNPs+INDELs)
    gatk Mutect2 \
        -R "${fasta}" \
        -I "${bam}" \
        -L "${bed}" \
        -tumor "${meta.id}" \
        --initial-tumor-lod 1.0 \
        --tumor-lod-to-emit 3.0 \
        --pcr-indel-model CONSERVATIVE \
        -O "${meta.id}.mutect2.vcf.gz"

    # 2) Apply Mutect2 filters
    gatk FilterMutectCalls \
        -R "${fasta}" \
        -V "${meta.id}.mutect2.vcf.gz" \
        -O "${meta.id}.mutect2.filtered.vcf.gz"

    # 3) Sort + reheader
    bcftools sort -O z -o "${meta.id}.mutect2.sorted.vcf.gz" "${meta.id}.mutect2.filtered.vcf.gz"
    tabix -f -p vcf "${meta.id}.mutect2.sorted.vcf.gz"

    bcftools reheader -f "${fai}" \
        "${meta.id}.mutect2.sorted.vcf.gz" \
        -o "${meta.id}.mutect2.reheader.vcf.gz"
    tabix -f -p vcf "${meta.id}.mutect2.reheader.vcf.gz"

    # 4) Normalize (keep SNPs+INDELs)
    bcftools norm -f "${fasta}" -m -both --check-ref w \
        -Oz -o "${meta.id}.mutect2.norm.vcf.gz" \
        "${meta.id}.mutect2.reheader.vcf.gz"
    tabix -f -p vcf "${meta.id}.mutect2.norm.vcf.gz"

    # 5) PASS only
    bcftools view -f .,PASS \
        -Oz -o "${meta.id}.mutect2.norm.pass.vcf.gz" \
        "${meta.id}.mutect2.norm.vcf.gz"
    tabix -f -p vcf "${meta.id}.mutect2.norm.pass.vcf.gz"

    """
}
// ------------------------------
// FREEBAYES CALL (reheader + normalize + PASS)
// ------------------------------
process FREEBAYES_CALL {
  tag { meta.id }
  publishDir "${params.outdir}/vcf/freebayes", mode:'copy'

  input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

  output:
    tuple val(meta),
          path("${meta.id}.freebayes.norm.pass.vcf.gz"),
          path("${meta.id}.freebayes.norm.pass.vcf.gz.tbi")

  script:
  """
  set -euo pipefail

  # 1) Raw FreeBayes calling
  freebayes -f "${fasta}" --min-alternate-fraction 0.01 "${bam}" \
      > "${meta.id}.freebayes.vcf"

  # 2) Reheader FIRST: inject ##contig entries from .fai
  bcftools reheader -f "${fai}" "${meta.id}.freebayes.vcf" \
      -o "${meta.id}.freebayes.reheader.vcf.gz"

  # 3) Sort & re-index
  bcftools sort -O z \
      -o "${meta.id}.freebayes.reheader.sorted.vcf.gz" \
      "${meta.id}.freebayes.reheader.vcf.gz"

  mv "${meta.id}.freebayes.reheader.sorted.vcf.gz" \
     "${meta.id}.freebayes.reheader.vcf.gz"

  tabix -f -p vcf "${meta.id}.freebayes.reheader.vcf.gz"

  # 4) Normalize
  bcftools norm -f "${fasta}" --check-ref w -m -both \
      -Oz -o "${meta.id}.freebayes.norm.vcf.gz" \
      "${meta.id}.freebayes.reheader.vcf.gz"

  tabix -f -p vcf "${meta.id}.freebayes.norm.vcf.gz"

  # 5) PASS-only
  bcftools view -f .,PASS -Oz \
      -o "${meta.id}.freebayes.norm.pass.vcf.gz" \
      "${meta.id}.freebayes.norm.vcf.gz"

  tabix -f -p vcf "${meta.id}.freebayes.norm.pass.vcf.gz"

  # Cleanup
  rm -f \
      "${meta.id}.freebayes.vcf" \
      "${meta.id}.freebayes.reheader.vcf.gz" \
      "${meta.id}.freebayes.reheader.vcf.gz.tbi" \
      "${meta.id}.freebayes.norm.vcf.gz" \
      "${meta.id}.freebayes.norm.vcf.gz.tbi" \
      || true
  """
}
// ------------------------------
// CONSENSUS (bcftools isec)
// ------------------------------
process ISEC_CONSENSUS {
  tag { meta.id }
  publishDir "${params.outdir}/vcf/consensus", mode:'copy', overwrite:true
  cpus 2
  memory '4 GB'

  input:
    tuple val(meta), path(vcfs)   // list: varscan_vcf, varscan_tbi, mutect2_vcf, mutect2_tbi, freebayes_vcf, freebayes_tbi

  output:
    tuple val(meta),
          path("${meta.id}.consensus.isec.vcf.gz"),
          path("${meta.id}.consensus.isec.vcf.gz.*i")
script:
  // Compute these in Groovy so we don't rely on Bash vars for interpolation
  def id          = meta.id
  def isecN       = (params.isec_n ?: '+2')
  // Join the staged files into one literal space-separated string
  def vcfList     = vcfs.collect{ it.toString() }.join(' ')
  def templateVcf = "${id}.mutect2.norm.pass.vcf.gz"

  """
  set -euo pipefail

  # 1) Ensure all files exist
  for f in ${vcfList}; do
    [ -f "\$f" ] || { echo "[ERROR] Missing file: \$f" >&2; exit 1; }
  done

  # 2) Build a list with ONLY the VCFs (*.vcf.gz), exclude index files
  vcf_list=""
  first=1
  for f in ${vcfList}; do
    case "\$f" in
      *.vcf.gz)
        if [ "\$first" -eq 1 ]; then
          vcf_list="\$f"; first=0
        else
          vcf_list="\$vcf_list \$f"
        fi
        ;;
    esac
  done

  [ -n "\$vcf_list" ] || { echo "[ERROR] No VCF files provided to ISEC_CONSENSUS" >&2; exit 1; }

  # 3) Compute consensus sites present in at least N callers
  bcftools isec -c all -n ${isecN} \$vcf_list -p "${id}.isec"

  # 4) Mutect2 template VCF must exist
  [ -f "${templateVcf}" ] || { echo "[ERROR] Missing template VCF: ${templateVcf}" >&2; exit 1; }

  # 5) Build the final consensus VCF by filtering Mutect2 with consensus sites
  if [ ! -s "${id}.isec/sites.txt" ]; then
    # No consensus sites: emit an empty, valid VCF
    echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" \
      | bgzip -c > "${id}.consensus.isec.vcf.gz"
  else
    bcftools view -R "${id}.isec/sites.txt" "${templateVcf}" \
      | bcftools sort -Oz -o "${id}.consensus.isec.vcf.gz"
  fi

  # 6) Index the consensus VCF
  tabix -f -p vcf "${id}.consensus.isec.vcf.gz" \
    || bcftools index -f "${id}.consensus.isec.vcf.gz"
  """
}
//-----------------------------
// ANNOVAR ANNOTATION (full + nonsyn-only + nonsyn VCF)
// ------------------------------
process ANNOVAR_ANNOTATE_AFDP_TLOD {
  tag { meta.id }
  publishDir "${params.outdir}/anno", mode: 'copy', overwrite: true
  cpus 3
  memory '6 GB'

  input:
    tuple val(meta),
          path(mut_vcf),        // Mutect2 VCF
          path(mut_vcf_index),  // Mutect2 index 
          path(annovar_bin_dir),
          path(annovar_db_dir)

  output:
    tuple val(meta),
          path("${meta.id}.annovar.full_multianno.txt")

  script:
  def id = meta.id

  """
  set -euo pipefail
  
# Ensure all ANNOVAR scripts are executable
 chmod +x ${annovar_bin_dir}/*.pl
  ###############################################
  # 1) Mutect2 VCF -> AVinput (used for annotation)
  ###############################################
  perl ${annovar_bin_dir}/convert2annovar.pl \
      -format vcf4 "${mut_vcf}" \
      > "${id}.avinput"

  ###############################################
  # 2) ANNOVAR annotation
  ###############################################
 perl  "${annovar_bin_dir}"/table_annovar.pl \
      "${id}.avinput" "${annovar_db_dir}" \
      -buildver CanFam3.1_ensembl \
      -out "${id}" \
      -remove \
      -protocol refGene \
      -operation g \
      -nastring . \
      -polish \
      -otherinfo

  mv "${id}.CanFam3.1_ensembl_multianno.txt" \
     "${id}.annovar.full_multianno.txt"

###############################################
# 3) Extract Mutect2 AF/DP/TLOD
###############################################
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TLOD\t[%AF]\t[%DP]\n' \
  "${mut_vcf}" > "${id}.mutect2.info.txt"

###############################################
# 4) Merge Mutect2 AF/DP/TLOD into multianno
###############################################
awk -F'\t' -v OFS='\t' '
  NR==FNR {
    key = \$1":"\$2":"\$3":"\$4
    af[key]   = \$6
    dp[key]   = \$7
    tlod[key] = \$5
    next
  }

  FNR==1 {
    print \$0, "AF", "DP", "TLOD"
    next
  }

  {
    key = \$1":"\$2":"\$4":"\$5
    print \$0, (af[key]?af[key]:"."), (dp[key]?dp[key]:"."), (tlod[key]?tlod[key]:".")
  }
' "${id}.mutect2.info.txt" "${id}.annovar.full_multianno.txt" \
  > "${id}.annovar.full_multianno.with_info.txt"

mv "${id}.annovar.full_multianno.with_info.txt" \
   "${id}.annovar.full_multianno.txt"
"""
}

//--------------------------
// ADD_GENE_SYMBOLS_FROM_CSV
// ------------------------------
process ADD_GENE_SYMBOLS_FROM_CSV {
  tag { meta.id }
  publishDir "${params.outdir}/anno", mode: 'copy', overwrite: true
  cpus 1
  memory '2 GB'

  input:
    tuple val(meta), path(multianno_full), path(add_gene_symbols_py), path(panel_csv)

  output:
    tuple val(meta),
          path("${meta.id}.annovar.full_multianno.with_symbols.txt")

  script:
  """
  set -euo pipefail

  # Run the Python script: it writes the target file and prints its path
  python3 ${add_gene_symbols_py} \
    "${multianno_full}" \
    "${panel_csv}" \
    "${meta.id}.annovar.full_multianno.with_symbols.txt" \
    > symbol_output_path.txt

  # Read back the printed path (one line)
  read out_path < symbol_output_path.txt

  # Ensure the expected file path exists (a no-op if already correct)
if [ ! -s "\${out_path}" ]; then
    echo "ERROR: Expected output not found: \${out_path}" >&2
    exit 1
fi

# Normalize the filename (again, a no-op if names already match)
if [ "\${out_path}" != "${meta.id}.annovar.full_multianno.with_symbols.txt" ]; then
    mv "\${out_path}" "${meta.id}.annovar.full_multianno.with_symbols.txt"
fi
"""
}


//--------------------------
// COLLAPSE_TRANSCRIPTS
// ------------------------------
process COLLAPSE_TRANSCRIPTS {
  tag { meta.id }
  publishDir "${params.outdir}/anno", mode: 'copy', overwrite: true
  cpus 1
  memory '2 GB'

  input:
    tuple val(meta), path(multianno)

  output:
    tuple val(meta),
          path("${meta.id}.annovar.gene_summary.tsv")

  script:
  """
  set -euo pipefail

  awk -F'\\t' '
    NR==1 {
      header=\$0
      for (i=1;i<=NF;i++) {
        if (\$i=="Gene.refGene") gcol=i
        if (\$i=="ExonicFunc.refGene") ecol=i
        if (\$i=="GeneSymbol") scol=i
      }
      next
    }

    {
      gene = (scol ? \$scol : \$gcol)
      effect = \$ecol

      # Assign severity score
      score=(effect=="stopgain"?3:(effect=="frameshift"?2:(effect=="nonsynonymous SNV"?1:0)))

      # Keep highest severity
      if (score > best[gene]) {
        best[gene]=score
        line[gene]=\$0
      }
    }

    END {
      print header
      for (g in line) print line[g]
    }
  ' "${multianno}" > "${meta.id}.annovar.gene_summary.tsv"
  """
}

//--------------------------
// WORKFLOW
// ------------------------------
workflow {

  // 0) FASTQC on raw reads
  fastqc_out = READS | FASTQC

  // Collect FASTQC artifacts for MultiQC (html + zip)
  fastqc_files = fastqc_out
                   .map { meta, html, zip -> [html, zip] }
                   .flatten()

  // 1) Alignment (classic bwa mem)
aligned = ALIGNMENT(
  READS.map { meta, read ->
    tuple(
      meta,
      read,
      fasta_file,
      bwa_bwt,
      bwa_ann,
      bwa_amb,
      bwa_pac,
      bwa_sa
    )
  }
)
  // 2) Add Read Groups
  rg_bams = aligned | ADD_READ_GROUP

  // 3) On-target BED filtering
  filtered = rg_bams
               .map { meta, rbam, rbai -> tuple(meta, rbam, rbai, bed_file) }
               | BAM_FILTER

  // 4) QC & Coverage (coverage hist NOT sent to MultiQC by choice)
  flagstat_out = filtered | SAMTOOLS_FLAGSTAT
  stats_out    = filtered | SAMTOOLS_STATS

  // Aggregate QC files for MultiQC (FASTQC + flagstat + stats)
  flagstat_files = flagstat_out.map { meta, f -> f }
  stats_files    = stats_out.map    { meta, s -> s }
  qc_files       = fastqc_files.mix(flagstat_files).mix(stats_files).collect()

  // 5) Variant callers (consume filtered RG BAMs)
  varscan_out =  filtered
    .map { meta, bam, bai ->
      tuple(meta, bam, bai, bed_file, fasta_file) }
    | VARSCAN_CALL

  mutect2_out = filtered
    .map { meta, fbam, fbai -> tuple(meta, fbam, fbai, bed_file, fasta_file, fasta_fai, fasta_dict)}
    | MUTECT2_CALL

  freebayes_out = filtered
                    .map { meta, fbam, fbai -> tuple(meta, fbam, fbai, fasta_file, fasta_fai) }
                    | FREEBAYES_CALL

  // 6) Consensus inputs: pass VCFs + their .tbi indexes to ISEC_CONSENSUS
// -------------------------------
// 7) Group PASS VCFs per sample and include both VCF + index
// -------------------------------

// Convert caller outputs to (sample_id, [vcf, tbi])
varscan_pairs = varscan_out.map { meta, vcf, idx ->
    tuple(meta.id, [vcf, idx])
}

mutect2_pairs = mutect2_out.map { meta, vcf, idx ->
    tuple(meta.id, [vcf, idx])
}

freebayes_pairs = freebayes_out.map { meta, vcf, idx ->
    tuple(meta.id, [vcf, idx])
}

 // Mix all callers → group by sample → flatten into a single VCF list
merged_vcfs_per_sample =
    varscan_pairs
        .mix(mutect2_pairs)
        .mix(freebayes_pairs)
        .groupTuple()    // produces: (id, [ [v1,idx1], [v2,idx2], [v3,idx3] ])
        .map { id, pairs ->
            def files = pairs.flatten()   // [var_vcf, var_idx, mut_vcf, mut_idx, fb_vcf, fb_idx]
            tuple([id:id], files)
        }

// Run consensus on grouped VCFs
consensus = merged_vcfs_per_sample | ISEC_CONSENSUS

// =========================
// 8) ANNOVAR + AF/DP/TLOD (Mutect2-only)
// =========================
// mutect2_out upstream emits: (meta, vcf, idx)
anno_af_out = mutect2_out
  .map { meta, vcf, idx ->
      tuple(
        meta,                // val(meta)
        vcf,                 // path(mut_vcf)
        idx,                 // path(mut_vcf_index) 
        annovar_bin_dir,
        annovar_db_dir

      )
  }
  | ANNOVAR_ANNOTATE_AFDP_TLOD

// =========================
// 9) Add GeneSymbol (CSV map)
// =========================
// ADD_GENE_SYMBOLS_FROM_CSV expects: tuple(val(meta), path(multianno_full), path(add_gene_symbols_py),path(panel_csv))
symbols_out = anno_af_out
  .map { meta, full_multianno ->
      tuple(
        meta,
        full_multianno,
        add_gene_symbols_py,
        panel_csv         // "path/Target_panel_gene_symbol.csv"
      )
  }
  | ADD_GENE_SYMBOLS_FROM_CSV

// =========================
// 10) Collapse to gene-level
// =========================
// COLLAPSE_TRANSCRIPTS expects: tuple(val(meta), path(multianno))
gene_summary = symbols_out
  .map { meta, full_with_symbols ->
      tuple(meta, full_with_symbols)
  }
  | COLLAPSE_TRANSCRIPTS


    //
    // === 7) MultiQC
    //
    MULTIQC(qc_files)
}
