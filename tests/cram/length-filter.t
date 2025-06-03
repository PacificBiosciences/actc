Input ccs bam contains some reads below the min-ccs-length
  $ ${ACTC} --min-ccs-length=10000 ${TESTDIR}"/../data/tiny.clr.bam" "${TESTDIR}"/../data/tiny.ccs.bam tiny.actc.bam  --log-level WARN
  $ samtools view -c tiny.actc.bam
  50
