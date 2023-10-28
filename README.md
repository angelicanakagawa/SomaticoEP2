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
```bash
