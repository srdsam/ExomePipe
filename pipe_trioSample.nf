
/*
 * GGC Exome Analysis Pipeline
 * Run Pipeline Test:
 * nextflow run pipe.nf -profile slurm \
 * -c /efs/sam/configScripts/slurm.config \
 * --inputBed /efs/sam/Macrogen_HN00115050/SureSelect_v6.bed \
 * --inputFastqs /efs/sam/Macrogen_HN00115050/fastq \
 * --inputHpo /efs/sam/Macrogen_HN00115050/hpo/x.hpo
 * --inputPed /efs/sam/Macrogen_HN00115050/ped
 */

// Each params.x can be overwritten (e.g. --inputBed 2390DE23492.bed)

// Genotype and Phenotype data files
params.inputBed = "*.bed"
params.inputFastqs = "./fastq"
params.inputHpo = "./hpo"
params.inputPed = "*.ped"
params.targetProband = ""

// Config files for Exomiser and Slurm
params.inputSlurmConfig = "/efs/sam/configScripts/slurm.config"

// Name for publishing directories
params.outdir = "results"
params.outdirSarek = "results-sarek"
params.outdirExomiser = "results-exomiser"

// Convert to file obj
tProband = params.targetProband
inputPed_path = params.inputPed
inputBed_file = file(params.inputBed)
inputSlurm_file = file(params.inputSlurmConfig)
inputHPO_files = file(params.inputHpo)
fastqPairs_ch = Channel.fromFilePairs( "$params.inputFastqs/*_{1,2}.fastq.gz")


log.info """\
G G C   E X O M E    A N A L Y S I S    v 1.0 
==============================================
MultiQC Analysis  : $params.outdirSarek
Exomiser Analysis : $params.outdirExomiser
Variant Calling   : snpEff, HaplotypeCaller
Prioritser        : hiPhive

Runtime ~ 6h 45m

Parameters

- profile
-- inputBed
-- inputFastqs
-- inputHpo
-- inputSlurmConfig "/efs/sam/configScripts/slurm.config"
-- outdirSarek
-- outdirExomiser

"""

/* STAGE 1
 * Running nf-core/sarek
 */

// Renames FastQ files from _1 to _R1
process renameFastQs {

    input:
        tuple val(name), file(reads) from fastqPairs_ch
    output:
        tuple val(name), stdout, file("$name/${reads[0].baseName.replace('_1', '_R1_')}.gz"), file("$name/${reads[1].baseName.replace('_2', '_R2_')}.gz") into pairedFastQ_ch

    script:

        """
        mkdir $name
        mv ${reads[0]} $name/${reads[0].baseName.replace('_1', '_R1_')}.gz
        mv ${reads[1]} $name/${reads[1].baseName.replace('_2', '_R2_')}.gz
        pwd
        """
}

// Removes headers and chr tag from BED file
process cleanBed {

    input:
    file bed from inputBed_file

    output:
    file 'SureSelect_v6_nochr.bed' into bed_ch

    script:
    """
    tail -n +4 $bed > intermediary.bed
    sed 's/^chr//' intermediary.bed > SureSelect_v6_nochr.bed
    """
}

// Runs nf-core/sarek and publish results
process runPipe {
    publishDir "$params.outdir/$params.outdirSarek", mode: 'copy'
    errorStrategy 'finish'

    input:
    tuple val(proband), path(fastqDir), file(fastq1), file(fastq2) from pairedFastQ_ch
    file bed from bed_ch
    file slurmConfig from inputSlurm_file

    output:
    stdout into terminalPipeChannel
    path("sarek-$proband/VariantCalling/$proband/HaplotypeCaller/HaplotypeCaller_${proband}.vcf.gz") into sarekDir_ch
    path("sarek-$proband/VariantCalling/$proband/HaplotypeCaller/HaplotypeCaller_${proband}.vcf.gz.tbi") into tbi_ch
    file("sarek-$proband/*") into pub_ch
    """
    /efs/sam/bin/nextflow run nf-core/sarek \
    -profile docker,slurm \
    -c $slurmConfig \
    --input $fastqDir/$proband \
    --genome GRCh37 \
    --tools haplotypecaller,snpEff \
    --target_bed $bed \
    --outdir sarek-$proband \
    --max_cpus 8 --cpus 4 --max_memory '30.GB' 
    """

}

/* STAGE 2
 * Running Exomiser
 */

// Build analysis YAML file ?? Error when trying to overwrite params with flags ??
// Produce string of HPO data for use by Exomiser


process produceHpoString {
    input:
    file(hpo) from inputHPO_files

    output:
    stdout into hpoString_ch

    """
    tr '\n' ',' < $hpo > output.hpo
    cat output.hpo
    """
}

// Concat vcfs together
process produceConcatVCF {
    input:
    path(vcfSarek) from sarekDir_ch.collect()
    path(index) from tbi_ch.collect()

    val(hpo) from hpoString_ch

    output:
    tuple path('merged.vcf.gz'), val(hpo) into vcf_ch

    """
    bcftools merge -o merged.vcf.gz --output-type z $vcfSarek
    """

}

// Produces Exomiser YAML file
process produceExomiserYAML {
    input:
    tuple path(vcf), val(hpo) from vcf_ch

    output:
    path("final-${tProband}.yml") into exomixerYAML_ch
    path(vcf) into vcfPass_ch

    script:
    text = writeExomiserYaml(params.outdirExomiser, "$inputPed_path/${tProband}.ped", tProband, "\$PWD/$vcf", hpo) 
    """
    echo "$text" > final-${tProband}.yml
    """
}
/*
process produceExomiserBatch {
    input:
    file(yaml) from exomixerYAML_ch.collect()
    path(vcf) from vcfPass_ch

    output:
    tuple file('batch.txt'), file(yaml), path(vcf) into exomiserBatch_ch
    """
    echo "${yaml.join("\n")}" > batch.txt
    """
}
*/
// Runs Exomiser
process runExomiser {
    publishDir "$params.outdir/$params.outdirExomiser", mode: 'copy'
    errorStrategy 'finish'

    input:
    file(yaml) from exomixerYAML_ch.collect()
    path(vcf) from vcfPass_ch

    output:
    stdout into terminalExomiserChannel
    file("*") into exFile_ch

    script:

    """
    java -Xms2g -Xmx18g -jar /efs/sam/bin/exomiser-cli-12.1.0/exomiser-cli-12.1.0.jar \
    --analysis $yaml \
    --spring.config.location=/efs/sam/bin/exomiser-cli-12.1.0/
    """

}

bed_ch.subscribe onComplete: {  println "cleanBed: Bed Cleaned"  }
terminalPipeChannel.subscribe {  println "runPipe: $it"  }
terminalExomiserChannel.subscribe {  println "runExomiser: $it"  }


/* HELPER FUNCTIONS
 *
 */

def writeExomiserYaml(outdirExomiser, pedFile, proband, vcf, hpo) {
    fileContents = """
analysis:
  analysisMode: PASS_ONLY
  frequencySources:
  - THOUSAND_GENOMES
  - TOPMED
  - UK10K
  - ESP_AFRICAN_AMERICAN
  - ESP_EUROPEAN_AMERICAN
  - ESP_ALL
  - EXAC_AFRICAN_INC_AFRICAN_AMERICAN
  - EXAC_AMERICAN
  - EXAC_SOUTH_ASIAN
  - EXAC_EAST_ASIAN
  - EXAC_FINNISH
  - EXAC_NON_FINNISH_EUROPEAN
  - EXAC_OTHER
  - GNOMAD_E_AFR
  - GNOMAD_E_AMR
  - GNOMAD_E_EAS
  - GNOMAD_E_FIN
  - GNOMAD_E_NFE
  - GNOMAD_E_OTH
  - GNOMAD_E_SAS
  - GNOMAD_G_AFR
  - GNOMAD_G_AMR
  - GNOMAD_G_EAS
  - GNOMAD_G_FIN
  - GNOMAD_G_NFE
  - GNOMAD_G_OTH
  - GNOMAD_G_SAS
  genomeAssembly: hg19
  hpoIds: ${hpo.split(',').collect{ "'" + it + "'" }}
  inheritanceModes:
    AUTOSOMAL_DOMINANT: 0.1
    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0
    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1
    MITOCHONDRIAL: 0.2
    X_DOMINANT: 0.1
    X_RECESSIVE_COMP_HET: 2.0
    X_RECESSIVE_HOM_ALT: 0.1
  pathogenicitySources:
  - POLYPHEN
  - MUTATION_TASTER
  - SIFT
  ped: $pedFile
  proband: '${proband}'
  steps:
  - variantEffectFilter:
      remove:
      - UPSTREAM_GENE_VARIANT
      - INTERGENIC_VARIANT
      - REGULATORY_REGION_VARIANT
      - CODING_TRANSCRIPT_INTRON_VARIANT
      - NON_CODING_TRANSCRIPT_INTRON_VARIANT
      - SYNONYMOUS_VARIANT
      - DOWNSTREAM_GENE_VARIANT
      - SPLICE_REGION_VARIANT
  - frequencyFilter:
      maxFrequency: 2.0
  - pathogenicityFilter:
      keepNonPathogenic: true
  - inheritanceFilter: {}
  - omimPrioritiser: {}
  - hiPhivePrioritiser: {}
  vcf: ${vcf}

outputOptions:
  numGenes: 50
  outputFormats:
  - TSV-GENE
  - TSV-VARIANT
  - VCF
  - HTML
  - JSON
  outputPassVariantsOnly: false
  outputPrefix: ./${proband}

    """
    return fileContents
}