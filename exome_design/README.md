# Exome design files

The first batch is on SureSelect V4 and V5 arrays. Use V5 design for all.

### Create Probe and Target bed files. Also create anti-target bed file.
Download from here: `https://earray.chem.agilent.com/suredesign/search.htm`

`cut -f1-3 /Users/joshuaschmidt/Downloads/S04380110_hs_hg38/S04380110_Covered.bed | sed -e '/browser/d' -e '/track/d' | sort -k1,1 -k2,2n > Agilent_V5_probes.bed`

`grep ENS /Users/joshuaschmidt/Downloads/S04380110_hs_hg38/S04380110_Targets.txt > Agilent_V5_ENST_targets.txt`

`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz`

`bedtools intersect -a gencode.v38.annotation.gff3.gz -b Agilent_V5_probes.bed |\
awk -F'[;\t =]' 'BEGIN {OFS="\t";} {chr=$1; start=$4-1; end=$5; annotation=$3;gene_id=$20; transcript_id=$16; gene_type=$18; transcript_type=$22; transcript_name=$24;exon_number=$26; } {if(gene_type=="protein_coding" && transcript_type=="protein_coding" && annotation=="exon") print chr, start, end, gene_id, exon_number, transcript_id;} {if(gene_type!="protein_coding" && annotation=="exon") print chr, start, end, gene_id, exon_number, transcript_id;}' | sed -e 's/"//g' -e 's/\.[[:digit:]]//g' >  Agilent_V5_exons.bed


`# R
library(data.table)
transcripts <- fread('Agilent_V5_ENST_targets.txt',col.names=c("transcript_id"))
exons <- fread("Agilent_V5_exons.bed",col.names=c("chr","start","end","gene_name","exon","transcript_id"))
exons <- exons[transcript_id %in% transcripts$transcript_id]
fwrite(exons[,.(chr,start,end,gene_name)],file="Agilent_V5_exons_target_transcripts.bed",sep="\t",row.names=F,col.names=F,quote=F)
`

sort -k1,1 -k2,2n Agilent_V5_exons_target_transcripts.bed | bedtools merge -i stdin -d 50 -c 4 -o distinct > Agilent_V5_targets.tmp

`#R
library(data.table)
targets <- fread('Agilent_V5_targets.tmp')
targets[,number:=1:.N,by=V4]
targets[,V4:=paste(V4,number,sep="_")]
targets[,number:=NULL]
fwrite(targets,file="Agilent_V5_targets.bed",sep="\t",row.names=F,col.names=F,quote=F)
`
