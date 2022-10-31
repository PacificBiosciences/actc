  $ ${ACTC} ${TESTDIR}"/../data/tiny.clr.bam" "${TESTDIR}"/../data/tiny.ccs.bam tiny.actc.bam
  $ samtools view -c tiny.actc.bam
  68

  $ samtools view tiny.actc.bam > tiny.actc.sam
  $ diff ${TESTDIR}"/../data/tiny.actc_expected.sam" tiny.actc.sam
