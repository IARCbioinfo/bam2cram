.DELETE_ON_ERROR:

${PREF}.R1.fq.gz:
	wgsim ref.fa ${PREF}.R1.fq ${PREF}.R2.fq 2>${PREF}.wgsim.err >${PREF}.wgsim.log
	gzip ${PREF}.R2.fq
	gzip ${PREF}.R1.fq
$(PREF).bam:${PREF}.R1.fq.gz
	bwa mem -R '@RG\tID:${PREF}\tSM:${PREF}' ref.fa ${PREF}.R1.fq.gz ${PREF}.R2.fq.gz 2>${PREF}.bwa.err | samtools view -b - | samtools sort -o $@ -
	samtools index $@
all:$(PREF).bam

