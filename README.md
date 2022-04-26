# indelValidate

Validate indels for the CGI clinical pipeline

## Overview

## Dependencies

* [svaba 1.1.0](https://github.com/walaj/svaba)
* [imsindel 1.0.2](https://github.com/NCGG-MGC/IMSindel)
* [samtools 1.15](http://www.htslib.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run indelValidate.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`normalbam`|File|normal input .bam file
`normalbai`|File|normal input .bai file
`tumorbam`|File|tumor input .bam file
`tumorbai`|File|tumor input .bai file
`chromosome`|String|chromosome on which the indel is
`indelid`|String|the number of the indel in the gene, in format of $gene_$num
`startPos`|String|start of indel position
`endPos`|String|end of indel position


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`svaba.modules`|String|"svaba/1.1.0 hg38-bwa-index/0.7.17"|Required environment modules
`svaba.basename`|String|basename("~{tumorbam}",".bam")|Base name
`svaba.reference`|String|"$HG38_BWA_INDEX_ROOT/hg38_random.fa"|genome reference
`svaba.startPosAdj`|Int|startPos - 10|start of indel position, adjusted for window size
`svaba.endPosAdj`|Int|endPos + 10|end of indel position, adjusted for window size
`svaba.jobMemory`|Int|64|Memory allocated for this job (GB)
`svaba.threads`|Int|1|Requested CPU threads
`svaba.timeout`|Int|10|Hours before task timeout
`imsindel.basename`|String|basename("~{tumorbam}",".bam")|Base name
`imsindel.modules`|String|"samtools/1.15 imsindel/1.0.2 hg38-bwa-index/0.7.17"|Required environment modules
`imsindel.reference`|String|"$HG38_BWA_INDEX_ROOT/hg38_random.fa"|genome reference
`imsindel.startPosAdj`|Int|startPos - 5000|start of indel position, adjusted for window size
`imsindel.endPosAdj`|Int|endPos + 5000|end of indel position, adjusted for window size
`imsindel.jobMemory`|Int|64|Memory allocated for this job (GB)
`imsindel.threads`|Int|1|Requested CPU threads
`imsindel.timeout`|Int|10|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`svabaIndelVCF`|File|Validated INDELS in .vcf format
`imsindelOut`|File|List of validated indels in the window


## Commands
 This section lists the two parallel commands run by the indelValidate workflow: svaba and imsindel
 
 * Running svaba
 
 Svaba requires a .bed file (which we make by joining chromosome number with indel start and end positions (adjusted by 10bp).  Validation then requires tumor and normal bams, the bed file (for local reassembly) and the genome reference.
 
 
 
 		echo -e "~{chromosome}\t~{startPosAdj}\t~{endPosAdj}" > ~{tumorbam}.~{indelid}.bed
 
 		svaba run -a validate \
 			-p ~{threads} \
 			-G ~{reference} \
 			-t ~{tumorbam} \
 			-n ~{normalbam}  \
 			-k ~{tumorbam}.~{indelid}.bed 
 
  * Running imsindel
  
 imsindel requires a .bam of just reads around the site, so we make that with samtools to include chromosome as well as start and end positions (adjusted by 5000bp) and index it. We then run imsindel with the tumorbam and the genome reference.   
 
 
 		samtools view -b ~{tumorbam} ~{chromosome}:~{startPosAdj}-~{endPosAdj} > ~{basename}.~{indelid}.bam
 		samtools index ~{basename}.~{indelid}.bam
 
 		imsindel --bam ~{basename}.~{indelid}.bam \
 			--chr ~{chromosome} \
 			--indelsize 10000 \
 			--reffa ~{reference} \
 			--temp temp/ \
 			--outd imsindelOutput/
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
