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

## Add covariates like mappability and GC content to target bed....
`MAP="/share/ScratchGeneral/jossch/mappability/hg38"`
`REF="/share/ScratchGeneral/jossch/reference/gatk/hg38"`
`TARGETS="/share/ScratchGeneral/jossch/WEScnv/exome_design"`

`
module load bedops
module load bedtools
#### MAP
bedmap --echo --mean --delim '\t' "$TARGETS"/Agilent_V5_targets.bed "$MAP"/tmp.map.bed \
> "$TARGETS"/Agilent_V5_targets-MAP100.bed

#### GC

bedtools nuc -bed "$TARGETS"/Agilent_V5_targets-MAP100.bed \
-fi "$REF"/Homo_sapiens_assembly38.fasta \
| cut -f1-5,7 | sed '1d' >  "$TARGETS"/Agilent_V5_targets-MAP100-GC.bed

bedtools nuc -bed <(awk '{OFS="\t"}; {$2=$2-500; $3=$3+500; print}' "$TARGETS"/Agilent_V5_targets-MAP100-GC.bed) \
-fi "$REF"/Homo_sapiens_assembly38.fasta \
| cut -f1-6,8 | sed '1d' | awk '{OFS="\t"}; {$2=$2+500; $3=$3-500; print}' >  "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500.bed

#### WIDTH

awk 'BEGIN {OFS="\t";} {width=sprintf("%.3f",($3-$2)/1000);} {print $0, width}' "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500.bed > "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500-WD.bed

#### MERGED
awk 'BEGIN {OFS="\t";} {merged=0;} {if($4~/,/) merged=1;} {print $0, merged} ' "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500-WD.bed > "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500-WD-M.bed

#### TXT with header

awk 'BEGIN {OFS="\t"; print "CHR", "START", "END", "TARGET", "MAP100", "GC", "GC500", "WIDTH", "MERGED";} {print $0}' "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500-WD-M.bed > "$TARGETS"/Agilent_V5_targets-MAP100-GC-GC500-WD-M.txt
`



# NEWER MANE SELECT
bgzip -dc MANE.GRCh38.v1.0.refseq_genomic.gff.gz | awk -F'[\t;=]' 'BEGIN {OFS="\t";} {chr=$1;ann=$3;start_pos=$4-1;end_pos=$5;name=$18;} {if(ann=="exon") {print chr, start_pos, end_pos, name;}}' | sort -k1,1 -k2,2n | bedtools merge -c 4 -o distinct -d -50 > MANE.GRCh38.v1.0.exons.merged.bed
