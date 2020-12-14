configfile: "config/config.yaml"
configfile: "config/samples.yaml"

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

rule all:
  input:
  	expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=samples.sample_id),
  	expand("results/readDepth/{sample}.bin{binSize}.wig", sample=samples.sample_id, binSize=str(config["binSize"]))

rule bwa_mem:
    input:
        reads=["data/FASTQ/{sample}_R1.fastq.gz", "data/FASTQ/{sample}_R2.fastq.gz"]
        # lambda wildcards: config["samples"][wildcards.samples]
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=config["bwa_index"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.64.0/bio/bwa/mem"

rule samtools_index:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.64.0/bio/samtools/index"

rule read_counter:
	input:
		bam = "mapped/{sample}.bam",
		bam_index= "mapped/{sample}.bam.bai"
	output:
		"results/readDepth/{sample}.bin{binSize}.wig"
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=config["ncores"]
	log:
		"logs/readDepth/{sample}.bin{binSize}.log"
	shell:
		"{params.readCounter} {input.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
	input:
		tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
		#norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		#corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/{tumor}.RData"
	params:
		outDir="results/ichorCNA/{tumor}/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=config["ncores"]
	log:
		"logs/ichorCNA/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

