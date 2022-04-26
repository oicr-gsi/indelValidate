version 1.0

workflow msisensor {
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
				name: "vep/1.6.17",
				url: "https://useast.ensembl.org/info/docs/tools/vep/index.html"
			},
			{
				name: "imsindel/1.0.2",
				url: "https://github.com/NCGG-MGC/IMSindel"
			},
			{
				name: "samtools/0.1.19",
				url: "http://www.htslib.org/"
			}
		]
		output_meta: {
			imsindelOut: "List of validated indels in the window"
			svabaIndelVCF: "Validated INDELS in .vcf format",
			svabaIndelMAF: "Validated INDELS, annotated by VEP, in .maf format"
		}
	}
	output {
		File svabaIndelVCF = "validate.svaba.somatic.indel.vcf"
		File svabaIndelMAF = "~{tumorbam}.~{indelid}.maf"
		File imsindelOut = "~{chromosome}.out"
	}
}

task svaba {
	input {
		File normalbam 
		File tumorbam
		File normalbai 
		File tumorbai 
		String modules = "svaba/1.1.0 hg38-bwa-index/0.7.17 vcf2maf/1.6.17 hg38/p12 vep-hg38-cache/92 vep-hg38-exac/0.3"
		String basename = basename("~{tumorbam}", ".bam")
		String reference = ~{HG38_BWA_INDEX_ROOT}/hg38_random.fa
		String vepPath = ~{VEP_ROOT}/bin
		String vepData = ~{VEP_HG38_CACHE_ROOT}/.vep
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
		vepPath: "path to VEP"
		vepData: "path to VEP data"
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

		~{VCF2MAF_ROOT}/bin/vcf2maf --species homo_sapiens  \
			--ncbi-build GRCh38 --ref-fasta ~{reference} \
			--vep-path ~{vepPath}  --vep-data ~{vepData}  \
			--input-vcf validate.svaba.somatic.indel.vcf  \
			--output-maf ~{tumorbam}.~{indelid}.maf 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File svabaIndelVCF = "validate.svaba.somatic.indel.vcf"
		File svabaIndelMAF = "~{tumorbam}.~{indelid}.maf"
	}

	meta {
		output_meta: {
			svabaIndelVCF: "Validated INDELS in .vcf format",
			svabaIndelMAF: "Validated INDELS, annotated by VEP, in .maf format"
		}
	}
}

task imsindel {
	input {
		File tumorbam
		File tumorbai 
		String basename = basename("~{tumorbam}", ".bam")
		String modules = "samtools/0.1.19 imsindel/1.0.2 hg38-bwa-index/0.7.17"
		String reference = ~{HG38_BWA_INDEX_ROOT}/hg38_random.fa
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

		samtools view -b tumorbam ~{chrom}:~{startPosAdj}-~{endPosAdj} > ~{tumorbam}.~{indelid}.bam
		samtools index ~{tumorbam}.~{indelid}.bam

		imsindel --bam ~{tumorbam}.~{indelid}.bam \
			--chr ~{chromosome} \
			--indelsize 10000 \
			--reffa ~{reference}

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