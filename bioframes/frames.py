from bioframes import makeframe, load_vcf
import sequenceframes
import samframe
from func import compose
from operator import methodcaller
load_fastq = compose(makeframe, sequenceframes.fqframe)

issam = methodcaller('endswith', 'sam')
isfastq = methodcaller('endswith', ('fastq', 'fq'))
isfasta = methodcaller('endswith', ('fasta', 'fa'))
isvcf = methodcaller('endswith', 'vcf')
ispileup = methodcaller('endswith', 'pileup')

def load_frame(fh):
    return load_by_ext(fh.name)(fh)

def load_by_ext(fn):
    if issam(fn): return samframe.load_sam
    if isfastq(fn): return load_fastq
    if isvcf: return load_vcf
    if ispileup: return samframe.load_pileup



