process CRISPRCLEANR_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'
    container "docker://biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1"

    input:
    tuple val(meta), path(count_file), val(library_file)
    val(min_reads)
    val(min_targeted_genes)

    output:
    tuple val(meta), path("*_norm_table.tsv"), emit: norm_count_file
    tuple val(meta), path("*.RData"),          emit: counts_rdata
    path "versions.yml",                       emit: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(CRISPRcleanR)
    library(dplyr)
    data('${library_file}')
    count_file <- read.delim('${count_file}',header=T,sep = "\t")
    count_file_to_normalize <- count_file  %>% dplyr::left_join(get('${library_file}'), by=c("sgRNA"="Target.Context.Sequence"),multiple = "all")

    count_file_to_normalize <- count_file_to_normalize %>% 
        dplyr::select(colnames(count_file),CODE,-sgRNA)

    names(count_file_to_normalize)[names(count_file_to_normalize) == 'Gene'] <- 'gene'
    names(count_file_to_normalize)[names(count_file_to_normalize) == 'CODE'] <- 'sgRNA'
    count_file_to_normalize <- count_file_to_normalize %>% dplyr::select(sgRNA, gene, everything())

    #crisprcleanr function
    normANDfcs <- ccr.NormfoldChanges(Dframe=count_file_to_normalize,saveToFig = FALSE,min_reads=${min_reads},EXPname="${prefix}", libraryAnnotation=get('${library_file}'),display=FALSE)
    gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs[["logFCs"]],get('${library_file}'))
    correctedFCs <- ccr.GWclean(gwSortedFCs,display=FALSE,label='${meta}')
    correctedCounts <- ccr.correctCounts('${meta}',
                            normANDfcs[["norm_counts"]],
                            correctedFCs,
                            get('${library_file}'),
                            minTargetedGenes=${min_targeted_genes},
                            OutDir='./')

    write.table(correctedCounts, file=paste0("${prefix}","_norm_table.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

    #version
    version_file_path <- "versions.yml"
    version_crisprcleanr <- paste(unlist(packageVersion("CRISPRcleanR")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    crisprcleanr: ", f, sep = "")
    writeLines(version_crisprcleanr, f)
    close(f)
    """
}
process BAGEL_FC {
    label "bagel2"
    errorStrategy { task.attempt > 1 ? 'ignore': 'retry' }
    maxRetries 3
    module "/software/team113/modules/modulefiles/bagel2/f9eedca"
    publishDir "${params.outdir}/fold_changes", mode: params.publish_dir_mode

    input: 
        tuple val(meta), path(input)
    output: 
        tuple val(meta), path("*.foldchange")
    script: 
        """
        BAGEL.py fc \
        -i $input \
        -c 1 \
        -o $meta.id
        """

}

process BAGEL_BF {
    label "bagel2"
    errorStrategy { task.attempt > 1 ? 'ignore': 'retry' }
    maxRetries 3
    module "/software/team113/modules/modulefiles/bagel2/f9eedca"
    publishDir "${params.outdir}/bayes_factors", mode: params.publish_dir_mode

    input: 
        tuple val(meta), path(input), val(columns)
        path(neg_list)
        path(ceg_list)
    output: 
        tuple val(meta), path("*bf.txt"), emit: bayes_factors
    script:
    """
    BAGEL.py bf \
    -i ${input} \
    -c ${columns} \
    -n ${neg_list} \
    -e ${ceg_list} \
    -o ${meta.id}_bf.txt
    """


}

process BAGEL_PR {
    label "bagel2"
    errorStrategy { task.attempt > 1 ? 'ignore': 'retry' }
    maxRetries 3
    module "/software/team113/modules/modulefiles/bagel2/f9eedca"
    publishDir "${params.outdir}/pr_curves", mode: params.publish_dir_mode

    input: 
        tuple val(meta), path(input)
        path(neg_list)
        path(ceg_list)

    output: 
        tuple val(meta), path("*pr.txt"), emit: pr_curve
    script:
    """
    BAGEL.py pr \
    -i ${input} \
    -n ${neg_list} \
    -e ${ceg_list} \
    -o ${meta.id}_pr.txt
    """

}

workflow {
    sgRNA_ch = Channel.fromPath("corrected_counts/*_sgRNA_count.tsv")
    .map {file -> tuple([id: file.baseName.split("_sgRNA_count")[0]],file)}
    BAGEL_FC(sgRNA_ch)

    fc_input = BAGEL_FC.out.map { meta, inputFile ->
    // Read the first line of the input file using Groovy. Note: 'text' reads the file content.
    def firstLine = inputFile.text.readLines()[0]
    def numCols = firstLine.tokenize().size() - 2
    def columns = (1..numCols).join(',')
    // Return a tuple that includes the computed columns.
    return [meta, inputFile, columns]}

    neg_list = file("inputs/NEGv1.txt")
    ceg_list = file("inputs/CEGv2.txt")
    fc_input.view()
    BAGEL_BF(fc_input, neg_list, ceg_list)
    BAGEL_PR(BAGEL_BF.out.bayes_factors, neg_list, ceg_list)
    
}