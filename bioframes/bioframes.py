from __future__ import print_function
from functools import partial
from Bio import SeqIO
import numpy as np
from operator import attrgetter as attr
from operator import add, div, itemgetter, methodcaller
from func import partial2, compose, pmap, pifilter,  \
    psplit,  _id, compose_all, ilen, merge_dicts, starcompose, dictzip, \
    apply_to_object, dictmap, kstarcompose, cmp2, pjoin, nameddict, \
    iter_until_stop, flatten_list
from itertools import starmap, repeat
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import pandas as pd
import vcf
#TODO:

def check_np_type(_type_str):
    def inner(array):
        return array.dtype == np.dtype(_type_str)
    return inner

'''
Join fastq and SAM (merge on QNAME [and SEQ])
Join VCF and SAM (merge on POS)
Pileup
Join VCF and Pileup
'''

get_fastq = partial(SeqIO.parse, format='fastq')
get_fasta = partial(SeqIO.parse, format='fasta')
to_np_int = partial(np.array, dtype=int)
gccontent = compose(ilen, pifilter('GC'.__contains__))

minus33 = partial(add, -33)
qual_int_sanger = compose(minus33, ord)

''' Error = 10^-(Phred/10) '''
qual_to_phreds = compose(to_np_int, pmap(qual_int_sanger))
error = compose(partial(pow, 10), partial2(div, -10.0))
#TODO: handle non-sanger version
sanger_qual_str_to_error = compose(error, qual_to_phreds)
#SANGER_OFFSET = 33

def flatten_vcf(record):

    '''
    :param VCFRecord record: VCFRecord object
    :return dict: all fields as a flat dictionary, including those fields in rec.INFO
    '''
    _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID'] #'INFO'
    fields = [getattr(record, _id) for _id in _ids]
    d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
    d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
    return d

#TODO: properly handle [None], '-', etc.
def load_vcf(vcf_records):
    '''
    Convert a list of vcf Records to a pandas DataFrame.
    '''
    return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)
load_vcf  = compose(load_vcf, vcf.Reader)

return_frame = partial(compose, pd.DataFrame.__getitem__)

def col_compare(df, col, value, comp):
    half = partial(comp, value)
    boolean = compose(half, df.__getitem__)
    return boolean(col)

simple_compare = return_frame(col_compare)

def ambiguous(df):
    return df[~df.CB.isin(list('AGCT'))]

def vcalls(df):
    return df[df.REF != df.CB]

def makeframe(biodata):
    dicts = iter_until_stop(_create_dict, **biodata)
    return pd.DataFrame(dicts)

def obj_to_dict(obj, names_getters_map):
    apply_to_obj = partial2(apply_to_object, obj)
    return dictmap(apply_to_obj, names_getters_map)

def _create_dict(obj_func, columns, getters, validator=None, dictgetters=None):
    obj = apply(obj_func)
    pre_dict = obj_to_dict(obj, dictzip(columns, getters))
    #optional intermediate schema here
    if dictgetters:
        extra_dict = obj_to_dict(pre_dict, dictgetters)
        full_dict = merge_dicts(pre_dict, extra_dict)
    else:
        full_dict = pre_dict
    return full_dict if not validator else validator.validate(full_dict)


make_writer = lambda outpath, format: partial(SeqIO.write, format=format, handle=outpath)
def write_result_seqs(outpath, format, func, **kwargs):
    ''' write the result of a function that returns Bio.SeqRecord.SeqRecord objects. '''
    with open(outpath, 'w') as out:
        writer = make_writer(out, format)
        writer(func( **kwargs))
    return outpath

#TODO: how can we have the outfacing function accept both the output filehandle and the df/input handle
raw_make_id = '{0} {1}'.format
cbs = itemgetter('CB')
consensus = compose(pjoin(''), cbs)
chrom_groups = partial2(pd.DataFrame.groupby, 'CHROM')

def prepend_value(recs, prefix, name=None):
    ''' Helper to alter sequence ids. BEWARE: alters recs objects. '''
    altered_id = compose(partial(raw_make_id, prefix), attr(name))
    recs = tuple(recs)
    tuple(starmap(setattr, zip(recs, repeat(name),  map(altered_id, recs))))
    return recs

prepend_id = partial(prepend_value, name='id')

def make_consensus(_id, group):
    return SeqRecord(Seq(consensus(group), generic_dna), _id, description='')

def vcf(in_path):
    ''' Write consensus creates a consensus fasta file (from the called bases) for each chromosome. '''
    df = load_vcf(in_path)
    def get_consensus(prefix=None):
        groups = df.groupby('CHROM')
        seqgen = starmap(make_consensus, groups)
        return seqgen if not prefix else prepend_id(seqgen, prefix)

    write_consensus = partial(write_result_seqs, format='fasta', func=get_consensus)

    return nameddict('VCFFrame', {
        'write_consensus' : write_consensus,
        'df' : df
    })

def make_matrix(df, column):
    ''' Make a matrix from a column containing numpy arrays'''
    matrix = (df[column].ix[i] for i in range(len(df[column])))
    return pd.DataFrame(matrix, index=df.id ).fillna(0)

phredmatrix = partial(make_matrix, column='qual_ints')
errormatrix = partial(make_matrix, column='error')
qualtotals = compose(methodcaller('sum'), phredmatrix)
makepercent = lambda x: (1 - x) * float(100)

def seq_matrix(df, column):
    ''' A dataframe representing a 2D matrix of sequences or qualities. '''
    matrix = make_matrix(df, column)

    mget = partial(partial, matrix.__getitem__)
    t_plot = partial(matrix.T.plot, kind='bar')
    boxplot = partial(matrix.plot, kind='box')
    percentplot = compose(pd.Series.plot, partial(pd.Series.apply, func=makepercent))
    plot_mean_percent = compose(percentplot, matrix.mean)

    return nameddict('SeqMatrix', {
        'tplot' : t_plot,
        'boxplot' : boxplot,
        'meanplot' : plot_mean_percent,
        'percentplot' : percentplot,
        'm' : matrix
    })

