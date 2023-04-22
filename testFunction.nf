
process test {

    output:
    file 'test.yaml' into final_ch

    script:
    text = writeExomiserYaml("outdirExomiser", "pedFile", "proband", "sarekDir", "dfd,fhpo")
    """
    echo "$text" > test.yaml
    """
}

final_ch.subscribe {  println "it: $it"  }

def writeExomiserYaml(outdirExomiser, pedFile, proband, sarekDir, hpo) {
    fileContents = """
    analysis:
        genomeAssembly: hg19
        vcf: ${sarekDir}/VariantCalling/${proband}/HaplotypeCaller/HaplotypeCaller_${proband}.vcf.gz
        ped: ${pedFile}
        proband: ${proband}
        hpoIds: ${hpo.split(',').collect{ '"' + it + '"'}}
        inheritanceModes: {
                AUTOSOMAL_DOMINANT: 0.1,
                AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
                AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
                X_DOMINANT: 0.1,
                X_RECESSIVE_HOM_ALT: 0.1,
                X_RECESSIVE_COMP_HET: 2.0,
                MITOCHONDRIAL: 0.2
        }
        analysisMode: PASS_ONLY
        frequencySources: [
            THOUSAND_GENOMES,
            TOPMED,
            UK10K,

            ESP_AFRICAN_AMERICAN, ESP_EUROPEAN_AMERICAN, ESP_ALL,

            EXAC_AFRICAN_INC_AFRICAN_AMERICAN, EXAC_AMERICAN,
            EXAC_SOUTH_ASIAN, EXAC_EAST_ASIAN,
            EXAC_FINNISH, EXAC_NON_FINNISH_EUROPEAN,
            EXAC_OTHER,

            GNOMAD_E_AFR,
            GNOMAD_E_AMR,
            GNOMAD_E_EAS,
            GNOMAD_E_FIN,
            GNOMAD_E_NFE,
            GNOMAD_E_OTH,
            GNOMAD_E_SAS,

            GNOMAD_G_AFR,
            GNOMAD_G_AMR,
            GNOMAD_G_EAS,
            GNOMAD_G_FIN,
            GNOMAD_G_NFE,
            GNOMAD_G_OTH,
            GNOMAD_G_SAS
        ]
        
    pathogenicitySources: [POLYPHEN, MUTATION_TASTER, SIFT]
    steps: [ 
        variantEffectFilter: {
            remove: [
                FIVE_PRIME_UTR_EXON_VARIANT,
                FIVE_PRIME_UTR_INTRON_VARIANT,
                THREE_PRIME_UTR_EXON_VARIANT,
                THREE_PRIME_UTR_INTRON_VARIANT,
                NON_CODING_TRANSCRIPT_EXON_VARIANT,
                UPSTREAM_GENE_VARIANT,
                INTERGENIC_VARIANT,
                REGULATORY_REGION_VARIANT,
                CODING_TRANSCRIPT_INTRON_VARIANT,
                NON_CODING_TRANSCRIPT_INTRON_VARIANT,
                DOWNSTREAM_GENE_VARIANT
            ]
        },
        frequencyFilter: {maxFrequency: 2.0},
        pathogenicityFilter: {keepNonPathogenic: true},
        inheritanceFilter: {},
        omimPrioritiser: {},
        hiPhivePrioritiser: {runParams: 'human'},
        ]
    outputOptions:
        outputContributingVariantsOnly: false
        numGenes: 0
    
        outputPrefix: ${outdirExomiser}/${proband}
        outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
    """
    return fileContents
}