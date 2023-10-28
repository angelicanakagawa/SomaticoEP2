# SomaticoEP2
Trabalho de chamada de variantes somático 

Este repositório tem como objetivo reproduzir a execução de um pipeline para chamada de variantes somáticas, utilizando o GATK.

## Criação do diretório
Para este exercício, os dados ficarão em um diretório chamado "EP2". Para isso, utilize o comando abaixo:
```bash
 mkdir EP2
```
Entre neste diretório recém criado:
```bash
 cd EP2
```

## Download dos FASTqs
A primeira coisa que faremos será instalar o programa *sratoolkits*: 
```bash
brew install sratoolkit
pip install parallel-fastq-dump
wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
tar -zxvf sratoolkit.3.0.0-ubuntu64.tar.gz
export PATH=$PATH://workspace/somaticoEP2/EP2/sratoolkit.3.0.0-ubuntu64/bin/
echo "Aexyo" | sratoolkit.3.0.0-ubuntu64/bin/vdb-config -i
```

Agora, nós vamos baixar os arquivos. E para ser mais rápido, eles serão separados em blocos e o download será realizado em paralelo.
```bash
time parallel-fastq-dump --sra-id SRR8856724 \
--threads 10 \
--outdir ./ \
--split-files \
--gzip
```

## Instalação dos programas necessários (bwa, samtools, bedtools, gatk, vcftools)
Com os comandos abaixo, nós iremos instalar os programas que utilizaremos neste pipeline.
```bash
brew install bwa 
```
```bash
brew install samtools
```
```bash
brew install bedtools
```
```bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
unzip gatk-4.2.2.0.zip 
```
```bash
brew install vcftools
```

## Download do chr9 e indexação
Antes de fazer o mapeamento das reads com o chr9, precisamos ter exatamente o chr9 no genoma de referência hg19. Para isso, execute:
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
gunzip chr9.fa.gz
```

E para indexar a referência do chr9:
```bash
bwa index chr9.fa
samtools faidx chr9.fa
```

## Mapeamento contra o genoma de referência (somente o chr9)
Para realizar o mapeamento contra o chr9 no genoma de referência, devemos utilizar o bwa mem.
```bash
NOME=WP312; Biblioteca=Nextera; Plataforma=illumina;
bwa mem -t 10 -M -R "@RG\tID:$NOME\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" chr9.fa \
SRR8856724_1.fastq.gz SRR8856724_2.fastq.gz | samtools view -F4 -Sbu -@2 - | samtools sort -m4G -@2 -o WP312_sorted.bam
```
```bash
time samtools index WP312_sorted.bam
```

## Remoção das duplicatas

Para remover a duplicata de PCR:
```bash
time samtools rmdup WP312_sorted.bam WP312_sorted_rmdup.bam
time samtools index WP312_sorted_rmdup.bam
```

## Download dos arquivos de referência no hg19

Os arquivos de Referência: Panel of Normal (PoN), Gnomad AF, se encontram em:
https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37?project=broad-dsde-outreach

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf
```
```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx
```
```bash
wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf
```
```bash
wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx
```

## Cobertura
Para gerar o bed de cobertura, devemos:
```bash
bedtools bamtobed -i WP312_sorted_rmdup.bam > WP312_sorted_rmdup.bed
bedtools merge -i WP312_sorted_rmdup.bed > WP312_sorted_rmdup_merged.bed
bedtools sort -i WP312_sorted_rmdup_merged.bed > WP312_sorted_rmdup_merged_sorted.bed
```
```bash
bedtools coverage -a WP312_sorted_rmdup_merged_sorted.bed \
-b WP312_sorted_rmdup.bam -mean \
> WP312_coverageBed.bed
```
Gerando o bed de cobertura 20X
```bash
cat WP312_coverageBed.bed | \
awk -F "\t" '{if($4>=20){print}}' \
> WP312_coverageBed20x.bed
```

## Adicionando o chr nos arquivos de referência
```bash
grep "\#" af-only-gnomad.raw.sites.vcf > af-only-gnomad.raw.sites.chr.vcf
grep  "^9" af-only-gnomad.raw.sites.vcf |  awk '{print("chr"$0)}' >> af-only-gnomad.raw.sites.chr.vcf
bgzip af-only-gnomad.raw.sites.chr.vcf
tabix -p vcf af-only-gnomad.raw.sites.chr.vcf.gz
```
```bash
grep "\#" Mutect2-WGS-panel-b37.vcf > Mutect2-WGS-panel-b37.chr.vcf 
grep  "^9" Mutect2-WGS-panel-b37.vcf |  awk '{print("chr"$0)}' >> Mutect2-WGS-panel-b37.chr.vcf 
bgzip Mutect2-WGS-panel-b37.chr.vcf 
tabix -p vcf Mutect2-WGS-panel-b37.chr.vcf.gz
```

## Calculando a contaminação
Antes, vamos gerar o arquivo .dict
```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R chr9.fa -O chr9.dict
```
```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R chr9.fa -O chr9.interval_list -OT ACGT
```
```bash
./gatk-4.2.2.0/gatk BedToIntervalList -I WP312_coverageBed20x.bed \
-O WP312_coverageBed20x.interval_list -SD chr9.dict
```
Agora gerando a tabela de contaminação:
```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I WP312_sorted_rmdup.bam  \
	-V af-only-gnomad.raw.sites.chr.vcf.gz  \
	-L WP312_coverageBed20x.interval_list \
	-O WP312.table
```

```bash
 ./gatk-4.2.2.0/gatk CalculateContamination \
-I WP312.table \
-O WP312.contamination.table
```

## GATK Mutect2
Para realizar a chamada de variantes, utilize o comando:
```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R chr9.fa \
  -I WP312_sorted_rmdup.bam \
  --germline-resource af-only-gnomad.raw.sites.chr.vcf.gz  \
  --panel-of-normals Mutect2-WGS-panel-b37.chr.vcf.gz \
  --disable-sequence-dictionary-validation \
  -L WP312_coverageBed20x.interval_list \
  -O WP312.somatic.pon.vcf.gz
```
```bash
  ./gatk-4.2.2.0/gatk FilterMutectCalls \
-R chr9.fa \
-V WP312.somatic.pon.vcf.gz \
--contamination-table WP312.contamination.table \
-O WP312.filtered.pon.vcf.gz
```

## Comparando os vcfs
Baixe os vcfs do endereço do drive abaixo (da amostra WP312):
  https://drive.google.com/drive/folders/1m2qmd0ca2Nwb7qcK58ER0zC8-1_9uAiE

É necessário alterar o arquivo para adicionar a palavra chr:
```bash
zgrep "\#" WP312.filtered.vcf.gz > header.txt
zgrep -v "\#" WP312.filtered.vcf.gz | awk '{print("chr"$0)}' > variants.txt
cat header.txt variants.txt > WP312.filtered.chr.vcf
bgzip WP312.filtered.chr.vcf
tabix WP312.filtered.chr.vcf.gz
```

Filtrando o chr9
```bash
zgrep "^\#\|chr9" WP312.filtered.chr.vcf.gz > WP312.filtered.chr9.vcf
bgzip  WP312.filtered.chr9.vcf
tabix WP312.filtered.chr9.vcf.gz
```

E por fim, realizar a comparação:
```bash
vcf-compare WP312.filtered.pon.vcf.gz WP312.filtered.chr9.vcf.gz
```

O resultado final deve ser esse:
```bash
# This file was generated by vcf-compare.
# The command line was: vcf-compare(v0.1.14-12-gcdb80b8) WP312.filtered.pon.vcf.gz WP312.filtered.chr9.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are: 
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      169     WP312.filtered.chr9.vcf.gz (25.8%)      WP312.filtered.pon.vcf.gz (0.2%)
VN      487     WP312.filtered.chr9.vcf.gz (74.2%)
VN      77879   WP312.filtered.pon.vcf.gz (99.8%)
#SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      Number of REF matches:  168
SN      Number of ALT matches:  166
SN      Number of REF mismatches:       1
SN      Number of ALT mismatches:       2
SN      Number of samples in GT comparison:     0
```

