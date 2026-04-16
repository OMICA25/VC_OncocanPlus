/*
 * VC_OncoCanPlus v1.0 - main.nf
 * Author: OMICA25
 * Target: CanFam3.1 panel sequencing (SE)
 * Execution: AWS EC2 (local executor) using conda env `ngs_env` and local ANNOVAR
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

// Index paths
fasta_fai  = file(params.fasta + ".fai")
fasta_dict = file(params.fasta.replaceFirst(/\\.[^.]+$/, "") + ".dict")

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
    tuple val(meta), path(read), path(ref_files)

  output:
    tuple val(meta),
          path("${meta.id}.sorted.bam"),
          path("${meta.id}.sorted.bam.bai")

  script:
  """
  set -euo pipefail

  export TMPDIR=\${TMPDIR:-\$PWD/tmp}
  mkdir -p "\$TMPDIR"

  echo "=== Reference files present ==="
  ls -lh
  echo "==============================="

  fasta=\$(ls *.fa)

  bwa mem \$fasta ${read} \\
    | samtools view -b - \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -m 1G \\
        -T "\$TMPDIR/${meta.id}.sort" \\
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
    tuple val(meta), path(bam), path(bai), val(bed)
  output:
    tuple val(meta),
          path("${meta.id}.on_target.bam"),
          path("${meta.id}.on_target.bam.bai")
  script:
  """
  set -euo pipefail
  samtools view -b -L ${bed} ${bam} > ${meta.id}.on_target.bam
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
// BEDTOOLS COVERAGE HISTOGRAM
// (NOT fed into MultiQC per your choice)
// ------------------------------
process BEDTOOLS_COVERAGE_HIST {
  tag { meta.id }
  publishDir "${params.outdir}/coverage", mode:'copy'
  input:
    tuple val(meta), path(bam), val(bed)
  output:
    tuple val(meta), path("${meta.id}.coverage.hist.tsv")
  script:
  """
  set -euo pipefail
  bedtools coverage -a ${bed} -b ${bam} -hist > ${meta.id}.coverage.hist.tsv
  """
}
// ------------------------------
// VARSCAN CALL (reheader + normalize + PASS)
// ------------------------------
process VARSCAN_CALL {
  tag { meta.id }
  publishDir "${params.outdir}/vcf/varscan", mode:'copy'

  input:
    tuple val(meta), path(bam), path(bai), val(bed), val(fasta)

  output:
    tuple val(meta),
          path("${meta.id}.varscan.norm.pass.vcf.gz"),
          path("${meta.id}.varscan.norm.pass.vcf.gz.tbi")

  script:
  """
  set -euo pipefail

  # 1) Raw calling (BED restricted)
  samtools mpileup -f ${fasta} -l ${bed} ${bam} \
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

  # 6) Cleanup
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
  publishDir "${params.outdir}/vcf/mutect2", mode:'copy'

  input:
    tuple val(meta), path(bam), path(bai), val(bed), val(fasta), path(fai), path(dict)

  output:
    tuple val(meta),
          path("${meta.id}.mutect2.norm.pass.vcf.gz"),
          path("${meta.id}.mutect2.norm.pass.vcf.gz.tbi")

  script:
  """
  set -euo pipefail

  # 1) Raw Mutect2 calling
  gatk Mutect2 \
      -R ${fasta} \
      -I ${bam} \
      -L ${bed} \
      -tumor ${meta.id} \
      -O ${meta.id}.mutect2.vcf.gz

  # 2) Mutect2 filters
  gatk FilterMutectCalls \
      -R ${fasta} \
      -V ${meta.id}.mutect2.vcf.gz \
      -O ${meta.id}.mutect2.filtered.vcf.gz

  # 3) Sort + bgzip + index
  bcftools sort -Oz \
      -o ${meta.id}.mutect2.sorted.vcf.gz \
      ${meta.id}.mutect2.filtered.vcf.gz
  tabix -f -p vcf ${meta.id}.mutect2.sorted.vcf.gz

  # 4) Reheader to add ##contig
  bcftools reheader -f ${fai} \
      ${meta.id}.mutect2.sorted.vcf.gz \
      -o ${meta.id}.mutect2.reheader.vcf.gz
  tabix -f -p vcf ${meta.id}.mutect2.reheader.vcf.gz

  # 5) Normalize
  bcftools norm -f ${fasta} --check-ref w -m -both \
      -Oz -o ${meta.id}.mutect2.norm.vcf.gz \
      ${meta.id}.mutect2.reheader.vcf.gz
  tabix -f -p vcf ${meta.id}.mutect2.norm.vcf.gz

  # 6) Keep only PASS
  bcftools view -f .,PASS -Oz \
      -o ${meta.id}.mutect2.norm.pass.vcf.gz \
      ${meta.id}.mutect2.norm.vcf.gz
  tabix -f -p vcf ${meta.id}.mutect2.norm.pass.vcf.gz

  # 7) Cleanup
  rm -f \
    ${meta.id}.mutect2.vcf.gz* \
    ${meta.id}.mutect2.filtered.vcf.gz* \
    ${meta.id}.mutect2.sorted.vcf.gz* \
    ${meta.id}.mutect2.reheader.vcf.gz* \
    ${meta.id}.mutect2.norm.vcf.gz* || true
  """
}

// ------------------------------
// FREEBAYES CALL (reheader + normalize + PASS)
// ------------------------------
process FREEBAYES_CALL {
  tag { meta.id }
  publishDir "${params.outdir}/vcf/freebayes", mode:'copy'

  input:
    tuple val(meta), path(bam), path(bai), val(fasta), path(fai)

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
process ANNOVAR_ANNOTATE {
  tag { meta.id }
  publishDir "${params.outdir}/anno",
             mode:'copy',
             overwrite:true,
             saveAs: { filename ->
               if( !params.keep_intermediates ) {
                 if ( filename.endsWith('.annovar.txt') ||
                      filename.endsWith('.annovar.nonsynonymous.tsv') ) {
                   return null
                 }
               }
               return filename
             }

  cpus 2
  memory '4 GB'

  input:
    tuple val(meta), path(vcf), path(vcf_index), val(annovar_home), val(annovar_db)

  output:
    tuple val(meta),
          path("${meta.id}.annovar.txt"),
          path("${meta.id}.annovar.nonsynonymous.tsv"),
          path("${meta.id}.consensus.isec.nonsynonymous.vcf.gz"),
          path("${meta.id}.consensus.isec.nonsynonymous.vcf.gz.tbi")

  script:
  """
  set -euo pipefail

  # 0) Keep only SNVs before annotation
  bcftools view -Ov -v snps ${vcf} > ${meta.id}.consensus.isec.snps.vcf

  # 1) Convert to ANNOVAR input
  ${annovar_home}/convert2annovar.pl -format vcf4 \
      ${meta.id}.consensus.isec.snps.vcf \
      > ${meta.id}.avinput

  # 2) Run table_annovar
  ${annovar_home}/table_annovar.pl \
      ${meta.id}.avinput ${annovar_db} \
      -buildver CanFam3.1_ensembl \
      -out ${meta.id} \
      -remove \
      -protocol refGene \
      -operation g \
      -nastring . \
      -polish \
      -otherinfo

  # 3) Rename output
  mv ${meta.id}.CanFam3.1_ensembl_multianno.txt ${meta.id}.annovar.txt

  # 4) Non-synonymous only
  awk -F'\t' 'NR==1 || /nonsynonymous SNV/ {print}' ${meta.id}.annovar.txt \
      > ${meta.id}.annovar.nonsynonymous.tsv

  # 5) BED of nonsynonymous SNVs (0-based)
  awk -F'\t' 'NR>1 && /nonsynonymous SNV/ {print \$1"\\t"(\$2-1)"\\t"\$2}' \
      ${meta.id}.annovar.txt \
      > ${meta.id}.nonsynonymous.sites.bed

  # 6) Generate nonsynonymous VCF
  bcftools view -R ${meta.id}.nonsynonymous.sites.bed ${vcf} \
      -Oz -o ${meta.id}.consensus.isec.nonsynonymous.vcf.gz
  tabix -f -p vcf ${meta.id}.consensus.isec.nonsynonymous.vcf.gz
  """
}

// ------------------------------
// COLLAPSE TRANSCRIPTS
// ------------------------------
process COLLAPSE_TRANSCRIPTS {
  tag { meta.id }
  publishDir "${params.outdir}/anno",
             mode:'copy',
             overwrite:true,
             saveAs: { filename ->
               if( !params.keep_intermediates ) {
                 if ( filename.endsWith('.annovar.gene_summary.tsv') ) {
                   return null
                 }
               }
               return filename
             }

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
    NR==1 { header=\$0; next }
    {
      gene=\$7; effect=\$9;
      score=(effect=="stopgain"?3:(effect=="frameshift"?2:(effect=="nonsynonymous SNV"?1:0)));
      if(score>best[gene]) { best[gene]=score; line[gene]=\$0 }
    }
    END {
      print header;
      for(g in line) print line[g];
    }
  ' ${multianno} > ${meta.id}.annovar.gene_summary.tsv
  """
}

// ------------------------------
// BUILD GENE MAP FROM BIOMART
// ------------------------------
process BUILD_GENE_MAP_BIOMART {
  tag { meta.id }
  publishDir "${params.outdir}/anno",
             mode:'copy',
             overwrite:true,
             saveAs: { filename ->
               if (!params.keep_intermediates && filename.endsWith('.gene_map.tsv'))
                 return null
               filename
             }

  cpus 1
  memory '1 GB'

  input:
    tuple val(meta), path(annovar_full), val(biomart_host), val(biomart_dataset)

  output:
    tuple val(meta), path("${meta.id}.gene_map.tsv")

  script:
    def id   = meta.id
    def full = annovar_full
    def host = biomart_host
    def ds   = biomart_dataset

  """
  set -euo pipefail

  # 1) Extract ENS IDs from Gene.refGene
  awk -F'\\t' '
    NR==1 {
      gcol = 0;
      for (i=1;i<=NF;i++) if (\$i=="Gene.refGene") gcol=i;
      if (!gcol) { print "ERROR: Gene.refGene column missing" > "/dev/stderr"; exit 1 }
      next;
    }
    {
      n = split(\$gcol, A, /[,; ]+/)
      for (i=1;i<=n;i++) if (A[i] ~ /^ENS/) ids[A[i]] = 1
    }
    END { for (k in ids) print k }
  ' ${full} | sort -u > ${id}.ens_ids.txt

  if [ ! -s ${id}.ens_ids.txt ]; then
    : > ${id}.gene_map.tsv
    exit 0
  fi

  # 2) Build comma-separated list
  ids_csv=\$(tr '\\n' ',' < ${id}.ens_ids.txt | sed 's/,\$//')

  # 3) Create BioMart XML query (Groovy interpolates ${ds})
  cat > biomart_query.xml <<XML
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1">
  <Dataset name="${ds}" interface="default">
    <Filter name="ensembl_gene_id" value="\${ids_csv}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>
XML

  # 4) Query BioMart
  curl -s --get --data-urlencode query@biomart_query.xml \
       ${host}/biomart/martservice \
       > ${id}.gene_map.tsv || true

  # 5) If nothing returned, map ENSID → ENSID
  if [ ! -s ${id}.gene_map.tsv ]; then
    awk '{print \$0"\\t"\$0}' ${id}.ens_ids.txt > ${id}.gene_map.tsv
  fi
  """
}

// ------------------------------
// ADD GENE SYMBOLS
// (Your original working version, Option A)
// -----------------------------
process ADD_GENE_SYMBOLS {
  tag { meta.id }
  publishDir "${params.outdir}/anno", mode: 'copy', overwrite: true
  cpus 1
  memory '1 GB'

  input:
    tuple val(meta),
          path(annovar_full),          // ${meta.id}.annovar.txt
          path(annovar_ns),            // ${meta.id}.annovar.nonsynonymous.tsv
          path(gene_map)               // ${meta.id}.gene_map.tsv

  output:
    tuple val(meta),
          path("${meta.id}.annovar.with_symbols.txt"),
          path("${meta.id}.annovar.nonsynonymous.with_symbols.tsv")

  script:
    // Capture Nextflow values as Groovy variables
    def id   = meta.id
    def full = annovar_full
    def ns   = annovar_ns
    def mapf = gene_map

    // Use triple-single quotes so Groovy doesn't parse the AWK regex.
    // Inject NF vars exactly once as shell assignments via concatenation.
    '''
set -euo pipefail

# ---------- Bring Nextflow values into shell vars (ONE TIME) ----------
id="''' + id + '''"
full="''' + full + '''"
ns="''' + ns + '''"
mapf="''' + mapf + '''"

# ---------- Normalize CRLF in mapping table (no-op on LF) ----------
perl -0777 -pe 's/\r\n/\n/g' "$mapf" > "$mapf.lf"
mv "$mapf.lf" "$mapf"

annotate_symbols () {
  in_file="$1"
  out_file="$2"
  mapfile="$3"

  # AWK stays literal; no Groovy interpolation inside a triple-single block
  awk -F'\t' -v OFS='\t' -v MAP="$mapfile" '
    BEGIN {
      mapped=0; unmapped=0;
      # Load ENS -> Symbol map
      while ((getline line < MAP) > 0) {
        n = split(line, f, /\t/)
        key = f[1]
        val = (n >= 2 ? f[2] : key)
        if (key != "") M[key] = val
      }
      close(MAP)
    }

    # Robust trim based on POSIX [[:space:]]
    function trim(s) {
      sub(/^[[:space:]]+/, "", s)
      sub(/[[:space:]]+$/, "", s)
      return s
    }

    FNR==1 {
      gcol=0
      for (i = 1; i <= NF; i++) if ($i=="Gene.refGene") gcol=i
      if (!gcol) { print "ERROR: Gene.refGene column not found" > "/dev/stderr"; exit 1 }
      print $0, "GeneSymbol"
      next
    }

    {
      genes = $gcol
      n = split(genes, a, /[;, ]+/)
      out = ""
      for (i = 1; i <= n; i++) {
        id = trim(a[i])
        sub(/:.*/, "", id)
        sym = (id in M ? M[id] : id)
        out = (out == "" ? sym : out "," sym)
        if (id in M) mapped++; else unmapped++
      }
      print $0, out
    }

    END {
      # could print stats to stderr if desired
      # print "ADD_GENE_SYMBOLS: mapped=" mapped ", unmapped=" unmapped > "/dev/stderr"
    }
  ' "$in_file" > "$out_file"
}

# Run twice: full and nonsyn tables
annotate_symbols "$full" "${id}.annovar.with_symbols.txt" "$mapf"
annotate_symbols "$ns"   "${id}.annovar.nonsynonymous.with_symbols.tsv" "$mapf"
'''
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
  READS.map { meta, read -> tuple(meta, read) }
       .combine(REF)
       .map { tuple(meta, read), ref_files -> tuple(meta, read, ref_files) }
)

  // 2) Add Read Groups
  rg_bams = aligned | ADD_READ_GROUP

  // 3) On-target BED filtering
  filtered = rg_bams
               .map { meta, rbam, rbai -> tuple(meta, rbam, rbai, params.bed) }
               | BAM_FILTER

  // 4) QC & Coverage (coverage hist NOT sent to MultiQC by choice)
  flagstat_out = filtered | SAMTOOLS_FLAGSTAT
  stats_out    = filtered | SAMTOOLS_STATS
  coverage_out = filtered
                   .map { meta, fbam, fbai -> tuple(meta, fbam, params.bed) }
                   | BEDTOOLS_COVERAGE_HIST

  // Aggregate QC files for MultiQC (FASTQC + flagstat + stats)
  flagstat_files = flagstat_out.map { meta, f -> f }
  stats_files    = stats_out.map    { meta, s -> s }
  qc_files       = fastqc_files.mix(flagstat_files).mix(stats_files).collect()

  // 5) Variant callers (consume filtered RG BAMs)
  varscan_out = filtered
                  .map { meta, fbam, fbai -> tuple(meta, fbam, fbai, params.bed, params.fasta) }
                  | VARSCAN_CALL

  mutect2_out = filtered
                  .map { meta, fbam, fbai -> tuple(meta, fbam, fbai, params.bed, params.fasta, fasta_fai, fasta_dict) }
                  | MUTECT2_CALL

  freebayes_out = filtered
                    .map { meta, fbam, fbai -> tuple(meta, fbam, fbai, params.fasta, fasta_fai) }
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

  // 7) ANNOVAR annotate (full + nonsyn VCF)
  anno_out = consensus
               .map { meta, vcf, tbi -> tuple(meta, vcf, tbi, params.annovar_home, params.annovar_db) }
               | ANNOVAR_ANNOTATE

  // 8) Collapse transcripts to gene-level summary (from full multianno)
  gene_summary = anno_out
                   .map { meta, ann_txt, ann_ns_tsv, ns_vcfgz, ns_tbi -> tuple(meta, ann_txt) }
                   | COLLAPSE_TRANSCRIPTS

  // 9) Build EnsemblID -> GeneSymbol map via BioMart (from full multianno)
  gene_map = anno_out
               .map { meta, ann_txt, ann_ns_tsv, ns_vcfgz, ns_tbi ->
                 tuple(meta, ann_txt, params.biomart_host, params.biomart_dataset)
               }
               | BUILD_GENE_MAP_BIOMART

  // 10) Add Gene Symbols to both ANNOVAR tables (join by sample id)
  with_symbols =
    gene_map
      .map { meta, mapTsv -> tuple(meta.id, mapTsv) }
      .join(
        anno_out.map { meta, annTxt, annNsTsv, nsVcf, nsTbi ->
          tuple(meta.id, [meta, annTxt, annNsTsv])
        }
      )
      .map { id, mapTsv, right ->
        def (metaA, annTxt, annNsTsv) = right
        tuple(metaA, annTxt, annNsTsv, mapTsv)
      }
      | ADD_GENE_SYMBOLS

  // 11) MultiQC (SAFE set: fastqc + flagstat + stats)
  MULTIQC(qc_files)
}
