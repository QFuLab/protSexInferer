#!/usr/bin/env nextflow

params.autoBLAST = false
params.minimumBLASTLength = 8
params.withoutBLASTdb = false 
params.findrAMELYThreshold = false
params.femaleMaxrAMELY = null
params.maleMinrAMELY = null

params.BLASTdbPath = "${workflow.projectDir}/BLAST_Database"
params.databaseSearchingSoftware = "PEAKS"

def installBlastAndDatabaseFunc() {

    println "Install BLAST using Conda"
    println "Download the BLAST database ..."
    
    def command = """
    mkdir ${workflow.projectDir}/BLAST_Database
    cd ${workflow.projectDir}/BLAST_Database
    perl ${workflow.projectDir}/script/toolkit/update_blastdb.pl --decompress nr > download_log 2>&1
    rm -f *.gz *.gz.md5
    """

    println "Running ${command}"
    def process = ["sh", "-c", command].execute()
    process.waitFor()
    
    if (process.exitValue() == 0) {
        println "Successfully install the BLAST and its database"
        return 0
    } else {
        return -1
    }
}

process MaxQuant2Peaks {
    publishDir "${params.outDir}/PeaksFormatInput", mode: "copy"
    
    input:
    tuple val(ind_name_input), path(raw_results)

    output:
    tuple val(ind_name), path("*.protein-peptides.csv")

    script:
    ind_name = ind_name_input
    """
    python3 ${workflow.projectDir}/script/FormatConvert/2PEAKSFormat.py ${raw_results} MaxQuant ${ind_name_input}
    """
}

process DiaNN2Peaks {
    publishDir "${params.outDir}/PeaksFormatInput", mode: "copy"
    
    input:
    tuple val(ind_name_input), path(raw_results)

    output:
    tuple val(ind_name), path("*.protein-peptides.csv")

    script:
    ind_name = ind_name_input
    """
    python3 ${workflow.projectDir}/script/FormatConvert/2PEAKSFormat.py ${raw_results} DIA-NN ${ind_name_input}
    """
}

process pFind2Peaks {
    publishDir "${params.outDir}/PeaksFormatInput", mode: "copy"
    
    input:
    tuple val(ind_name_input), path(raw_results)

    output:
    tuple val(ind_name), path("*.protein-peptides.csv")

    script:
    ind_name = ind_name_input
    """
    python3 ${workflow.projectDir}/script/FormatConvert/2PEAKSFormat.py ${raw_results} pFind ${ind_name_input}
    """
}

process autoBLAST {
    conda "blast"

    publishDir "${params.outDir}/classified_peptides", pattern: "*-classified_peptides_afterBLAST.csv", mode: "copy"
    publishDir "${params.outDir}/BLAST_results", pattern: "*-BLAST.report", mode: "copy"
    publishDir "${params.outDir}/BLAST_results", pattern: "*-removed_peptides.list", mode: "copy"
    publishDir "${params.outDir}/BLAST_results", pattern: "*-peptides-need-BLAST.fasta", mode: "copy"


    input:
    tuple val(ind_name_input), path(classified_peptides_result)
    val minimumBLASTLength

    output:
    tuple val(ind_name), path("*-classified_peptides_afterBLAST.csv"), emit: peptides_afterBLAST
    path("*-BLAST.report"), optional: true
    path("*-removed_peptides.list"), optional: true
    path("*-peptides-need-BLAST.fasta"), optional: true

    script:
    ind_name = ind_name_input
    """
    export BLASTDB=$params.BLASTdbPath

    python3 ${workflow.projectDir}/script/autoBLAST.py ${ind_name_input} ${classified_peptides_result} --minimumBLASTLength ${minimumBLASTLength}
    """
}

process classfyPeptides {
    publishDir "${params.outDir}/classified_peptides", mode: "copy"
    
    input:
    tuple val(ind_name_input), path(peptides_result)
    path RefDatabase

    output:
    tuple val(ind_name), path("*-classified_peptides.csv")

    script:
    ind_name = ind_name_input

    """
    python3 ${workflow.projectDir}/script/peptide_classifier.py ${ind_name_input} ${peptides_result} ${RefDatabase}
    """
}

process cal_rAMELY {    
    publishDir "${params.outDir}/rAMELY_ratios", mode: "copy"
    
    input:
    tuple val(ind_name_input), path(classified_peptides_result)

    output:
    path("*-rAMELY_summary_statistics.csv"), emit: peptides_csv

    script:
    """
    python3 ${workflow.projectDir}/script/rAMELY_computer.py ${ind_name_input} ${classified_peptides_result}
    """
}

process concatDf {
    cache false
    publishDir "${params.outDir}", mode: "copy", enabled: params.findrAMELYThreshold

    input:
        path inputFiles
    output:
        path "rAMELY_summary_statistics.csv", emit: rAMELY_out

    script:
    """
    python3 ${workflow.projectDir}/script/toolkit/concat.py ${inputFiles}
    """
}

process findThreshold {
    publishDir "${params.outDir}", mode: "copy"
    debug true

    input:
        val input_fn
        path rAMELY_results
        val outDir
    output:
        path "rAMELY_threshold.txt"
        path "rAMELY_distribution.pdf"

    script:
    """
    python3 ${workflow.projectDir}/script/findrAMELYThreshold.py ${input_fn} ${rAMELY_results}

    echo "All work completed! Please check the rAMELY_summary_statistics.csv, rAMELY_threshold.txt and rAMELY_distribution.pdf in the output directory." 
    echo "Finished time: \$(date)" 
    """
}

process plotrAMELY {
    publishDir "${params.outDir}", mode: "copy"
    debug true

    input:
        val input_fn
        val femaleMaxrAMELY
        val maleMinrAMELY
        val outDir

    output:
        path "rAMELY_SexDeterminer.pdf"
        path "Sex_assessment_report.csv"

    script:
    """
    python3 ${workflow.projectDir}/script/plotrAMELY.py ${input_fn} ${femaleMaxrAMELY} ${maleMinrAMELY}
    
    echo "All work completed! Please check the Sex_assessment_report.csv and rAMELY_SexDeterminer.pdf in the output directory." 
    echo "The classified AMELY/AMELX-specifc peptides and the rAMELY_ratios for each individual are also saved in their corresponding folders."
    echo "Finished time: \$(date)" 
    """
}

workflow {
    // Validate required parameters
    assert params.inputFile != null &&  params.outDir != null:
        "Error: You must provide input file (--inputFile) and output directory path (--outDir)!"
    assert params.databaseSearchingSoftware in ["PEAKS", "MaxQuant", "DIA-NN", "pFind"] :  
        "Error: The value of params.databaseSearchingSoftware should be one of 'PEAKS', 'pFind', 'MaxQuant', or 'DIA-NN'"

    // Set default thresholds if not provided
    def femaleMaxrAMELYDefault = ["PEAKS": 0.05525587, "MaxQuant": 0.01608375, "DIA-NN": 0.07319711, "pFind": 0.03290402]
    def maleMinrAMELYDefault = ["PEAKS": 0.08809166, "MaxQuant": 0.03103858, "DIA-NN": 0.08635405, "pFind": 0.05770555]

    if (params.femaleMaxrAMELY == null) {
        UserfemaleMaxrAMELY = femaleMaxrAMELYDefault[params.databaseSearchingSoftware]
    } else {
        UserfemaleMaxrAMELY = params.femaleMaxrAMELY
    }
    if (params.maleMinrAMELY == null) {
        UsermaleMinrAMELY = maleMinrAMELYDefault[params.databaseSearchingSoftware]
    } else {
        UsermaleMinrAMELY = params.maleMinrAMELY
    }

    //hello_message(params.inputFile)
    def currentDate = new java.util.Date()
    println "Sex Determiner (3 Dec 2025)"
    println "© 2025 Fan BAI, Zhongyou WU, and Qiaomei FU"
    println ""
    println "Options:"
    println "  --inputFile ${params.inputFile}"
    println "  --RefDatabase ${params.RefDatabase}"
    println "  --databaseSearchingSoftware ${params.databaseSearchingSoftware}"
    println "  --findrAMELYThreshold ${params.findrAMELYThreshold}"
    println "  --femaleMaxrAMELY ${UserfemaleMaxrAMELY}"
    println "  --maleMinrAMELY ${UsermaleMinrAMELY}"
    println "  --autoBLAST ${params.autoBLAST}"
    println "  --withoutBLASTdb ${params.withoutBLASTdb}"
    println "  --minimumBLASTLength ${params.minimumBLASTLength}"
    println "  --outDir ${params.outDir}"
    println ""
    println "Nextflow project directory is ${workflow.projectDir}"
    println "Start time: ${currentDate}"

    // Install BLAST and database
    if (params.autoBLAST && params.withoutBLASTdb) {
        println "Install BLAST and download the database"
        installBlastAndDatabaseFunc()
    } 

    // Parse input file and convert to PEAKS format if needed    
    raw_row = channel.fromPath(params.inputFile)
                                 .splitCsv()
                                 .map{ row -> [row[0], file(row[1])] }
    if (params.databaseSearchingSoftware == "PEAKS") {
        input_row = raw_row
    } else if (params.databaseSearchingSoftware == "MaxQuant") {
        input_row = MaxQuant2Peaks(raw_row).map{ output -> [output[0], file(output[1])] }
    } else if (params.databaseSearchingSoftware == "DIA-NN") {
        println "Convert DIA-NN results to PEAKS format ..."
        input_row = DiaNN2Peaks(raw_row).map{ output -> [output[0], file(output[1])] }
    } else if (params.databaseSearchingSoftware == "pFind") {
        println "Convert pFind results to PEAKS format ..."
        input_row = pFind2Peaks(raw_row).map{ output -> [output[0], file(output[1])] }
    }

    // Calculating the rAMELY values with or without automatic BLAST
    if (params.autoBLAST) {
        println "Run BLAST for classified peptides"
        println "Option in BLAST:"
        println "  --minimumBLASTLength ${params.minimumBLASTLength}"

        classfyPeptides(input_row, channel.value(file(params.RefDatabase)))
        autoBLAST(classfyPeptides.out, params.minimumBLASTLength)
        cal_rAMELY(autoBLAST.out.peptides_afterBLAST)
    } else {
        println "Without automatic BLAST for classified peptides"
        classfyPeptides(input_row, channel.value(file(params.RefDatabase)))
        cal_rAMELY(classfyPeptides.out)
    }

    // Concatenate all rAMELY results
    concatDf(cal_rAMELY.out.peptides_csv.collect())

    if (params.findrAMELYThreshold) {
        // Find the rAMELY threshold from the input data
        findThreshold(channel.fromPath(params.inputFile), concatDf.out.rAMELY_out, params.outDir)
    } else {
        // Plot the rAMELY distribution with user-defined thresholds
        plotrAMELY(concatDf.out.rAMELY_out, UserfemaleMaxrAMELY, UsermaleMinrAMELY, params.outDir)
    }
}
