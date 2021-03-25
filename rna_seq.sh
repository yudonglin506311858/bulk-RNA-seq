#质量控制（quality control）
fastqc -t 30 -o ./fastqc/ *.gz


#过滤接头（adapter removal）
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_rmadp_1.fastq.gz -p ${i}_rmadp_2.fastq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 /data/yudonglin/rnaseqdata/${i}_rmadp_1.fastq.gz /data/yudonglin/rnaseqdata/${i}_rmadp_2.fastq.gz -baseout /data/yudonglin/rnaseqdata/${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
done

#去接头之后再进行一次质量控制（quality control again）
fastqc -t 30 -o ./fastqc/ *.gz


#比对（alignment）
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4,B654-pro-19-12-24_FKDL202555457-1a-A5,B654-Baso-19-12-24_FKDL202555457-1a-A6,B654-poly-19-12-24_FKDL202555457-1a-A7,B654-ortho-19-12-24_FKDL202555457-1a-A8};
do hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 /data/yudonglin/RNA-seq2/${i}_fliter_1P.fq.gz -2 /data/yudonglin/RNA-seq2/${i}_fliter_2P.fq.gz -S /data/yudonglin/RNA-seq2/${i}.sam;
done

#sam文件转换为bam文件（sam to bam file）
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4,B654-pro-19-12-24_FKDL202555457-1a-A5,B654-Baso-19-12-24_FKDL202555457-1a-A6,B654-poly-19-12-24_FKDL202555457-1a-A7,B654-ortho-19-12-24_FKDL202555457-1a-A8};
do samtools view -S ${i}.sam -b > ${i}.bam
rm ${i}.sam
samtools sort ${i}.bam -o ${i}_sorted.bam
samtools index ${i}_sorted.bam;
done


#基因表达定量（gene expression matrix）
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4,B654-pro-19-12-24_FKDL202555457-1a-A5,B654-Baso-19-12-24_FKDL202555457-1a-A6,B654-poly-19-12-24_FKDL202555457-1a-A7,B654-ortho-19-12-24_FKDL202555457-1a-A8};
do featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/gencode.vM16.basic.annotation.gtf -o ${i}.count ${i}.bam >>~/count.txt 2>&1
rm ${i}.bam;
done
