version 1.0

workflow indelValidate {
	input {
		File normalbam
		File normalbai
		File tumorbam
		File tumorbai
		String chromosome
		String indelid
		String startPos
		String endPos
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		normalbam: "normal input .bai file"
		tumorbam: "tumor input .bam file"
		tumorbam: "tumor input .bai file"
		chromosome: "chromosome on which the indel is"
		indelid: "the number of the indel in the gene, in format of $gene_$num"
		startPos: "start of indel position"
		endPos: "end of indel position"
	}

	call svaba {
		input: 
			normalbam = normalbam,
			tumorbam = tumorbam,
			normalbai = normalbai, 
			tumorbai = tumorbai, 
			chromosome = chromosome,
			indelid = indelid,
			startPos = startPos,
			endPos = endPos
	}

	call imsindel {
		input: 
			tumorbam = tumorbam,
			tumorbai = tumorbai,
			indelid = indelid,
			chromosome = chromosome,
			startPos = startPos,
			endPos = endPos
	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Validate indels for the CGI clinical pipeline"
		dependencies: 
		[
			{
				name: "svaba/1.1.0",
				url: "https://github.com/walaj/svaba"
			},
			{
				name: "imsindel/1.0.2",
				url: "https://github.com/NCGG-MGC/IMSindel"
			},
			{
				name: "samtools/1.15",
				url: "http://www.htslib.org/"
			}
		]
		output_meta: {
			imsindelOut: "List of validated indels in the window",
			svabaIndelVCF: "Validated INDELS in .vcf format",
		}
	}
	output {
		File svabaIndelVCF = "validate.svaba.somatic.indel.vcf"
		File imsindelOut = "~{chromosome}.out"
	}
}

task svaba {
	input {
		File normalbam 
		File tumorbam
		File normalbai 
		File tumorbai 
		String modules = "svaba/1.1.0 hg38-bwa-index/0.7.17"
		String basename = basename("~{tumorbam}", ".bam")
		String reference = "$HG38_BWA_INDEX_ROOT/hg38_random.fa"
		String chromosome
		String indelid
		Int startPos
		Int endPos
		Int startPosAdj = (startPos-10)
		Int endPosAdj = (endPos+10)
		Int jobMemory = 64
		Int threads = 1
		Int timeout = 10
	}

	parameter_meta {
		normalbam: "normal input .bam file"
		normalbai: "normal input .bai file"
		tumorbam: "tumor input .bam file"
		tumorbai: "tumor input .bai file"
		basename: "Base name"
		modules: "Required environment modules"
		reference: "genome reference"
		chromosome: "chromosome on which the indel is"
		indelid: "the number of the indel in the gene, in format of $gene_$num"
		startPos: "start of indel position"
		endPos: "end of indel position"
		startPosAdj: "start of indel position, adjusted for window size"
		endPosAdj: "end of indel position, adjusted for window size"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		echo -e "~{chromosome}\t~{startPosAdj}\t~{endPosAdj}" > ~{tumorbam}.~{indelid}.bed

		svaba run -a validate \
			-p ~{threads} \
			-G ~{reference} \
			-t ~{tumorbam} \
			-n ~{normalbam}  \
			-k ~{tumorbam}.~{indelid}.bed 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File svabaIndelVCF = "validate.svaba.somatic.indel.vcf"
	}

	meta {
		output_meta: {
			svabaIndelVCF: "Validated INDELS in .vcf format",
		}
	}
}

task imsindel {
	input {
		File tumorbam
		File tumorbai 
		String basename = basename("~{tumorbam}", ".bam")
		String modules = "samtools/1.15 imsindel/1.0.2 hg38-bwa-index/0.7.17"
		String reference = "$HG38_BWA_INDEX_ROOT/hg38_random.fa"
		String chromosome
		String indelid
		Int startPos
		Int endPos
		Int startPosAdj = (startPos-5000)
		Int endPosAdj = (endPos+5000)
		Int jobMemory = 64
		Int threads = 1
		Int timeout = 10
	}

	parameter_meta {
		tumorbam: "tumor input .bam file"
		tumorbai: "tumor input .bai file"
		basename: "Base name"
		modules: "Required environment modules"
		reference: "genome reference"
		chromosome: "chromosome on which the indel is"
		indelid: "the number of the indel in the gene, in format of $gene_$num"
		startPos: "start of indel position"
		endPos: "end of indel position"
		startPosAdj: "start of indel position, adjusted for window size"
		endPosAdj: "end of indel position, adjusted for window size"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail



		samtools view -b ~{tumorbam} ~{chromosome}:~{startPosAdj}-~{endPosAdj} > ~{basename}.~{indelid}.bam
		samtools index ~{basename}.~{indelid}.bam

		mkdir imsindelOutput
		mkdir temp
		imsindel --bam ~{basename}.~{indelid}.bam \
			--chr ~{chromosome} \
			--indelsize 10000 \
			--reffa ~{reference} \
			--temp temp/ \
			--outd imsindelOutput/

		cp imsindelOutput/~{chromosome}.out ~{chromosome}.out

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File imsindelOut = "~{chromosome}.out"
	}

	meta {
		output_meta: {
			imsindelOut: "List of validated indels in the window"
		}
	}
}