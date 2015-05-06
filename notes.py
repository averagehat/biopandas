from bioframes.bioframes import *
from bioframes.samframe import *
from func import *
import pandas as pd
import numpy as np

df = pandas.DataFrame(np.random.normal(size=(75,5)))


fqdf.seq.apply(len).plot(kind='bar')
fqdf['mean_error'] = fqdf.error.apply(np.ndarray.mean)
fqdf[fqdf.error.apply(min) > 0.00001]
samdf.plot(kind='scatter', x='POS', y='cigar_score')


refactor useful functions out
ensure that VCF files work
try/design vcf filters (port tests)
replace map, filter, pmap, etc. with their lazy equivalents
port/figure out code for non-unique merges/joins
create demonstrative notebook
the duplicated thing in dataframes

load fasta files
could store sequence length
handle sam files with headers
allow for loading multiple fastq files at once (or really any file except maybe vcf)    using itertools.chain

varrecon
