from __future__ import print_function
from functools import partial
from Bio import SeqIO
import pandas as pd
import numpy as np
#from past.builtins import map, filter
from operator import attrgetter as attr
from operator import add, div, itemgetter
from func import partial, partial2, compose, pmap, pfilter, pstrip, psplit, fzip, _id, compose_all
from itertools import repeat
from schema import Schema
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

class BioFrame(pd.DataFrame):
    def find_duplicates(frame):
        ''' not sure why this one works. '''
        ''' for example, find pcr duplicates by finding duplicate reads'''
        #df2[df2['b'] == df2['b'] & (df2['a'] == df2['a'])]
        ''' duplicated by defaults returns all but the first occurrence; take_last takes only the first occurence,
        so unioning them (via | (or)), returns all duplicates. accepts any number of keys. '''
        frame[frame.duplicated(['b'],take_last=True) | frame.duplicated(['b'])]


to_np_int = partial(np.array, dtype=int)
minus33 = partial(add, -33)
qual_int_sanger = compose(minus33, ord)
''' Error = 10^-(Phred/10) '''
qual_to_phreds = compose(to_np_int, pmap(qual_int_sanger))
error = compose(partial(pow, 10), partial2(div, -10))
qual_to_error = compose(pmap(error), qual_to_phreds)


get_fastq = partial(SeqIO.parse, format='fastq')
def init_all(fileh):
    fq = get_fastq(fileh)
    dicts = map(get_row, fq)
    return pd.DataFrame(dicts).set_index(index) #, index=index, columns=columns)


#SANGER_OFFSET = 33
columns = ('id', 'seq', 'quality', 'description', 'qual_ints', 'error')
SANGER = True
get_id = attr('id')
get_seq= compose(str, attr('seq'))
get_qual_ints = compose_all(np.array, itemgetter('phred_quality'), attr('_per_letter_annotations'))
get_description = attr('description')
get_quality = SeqIO.QualityIO._get_sanger_quality_str
#get_error = compose_all(np.array, pmap(error), get_qual_ints)
get_error = compose_all(np.vectorize(error), get_qual_ints)
''' applies directly to the seqrecord object '''
index = ['id']

#TODO: fix for 32-bit systems
final_schema =  Schema({
    'id' : str,
    'seq' : str,
    'quality' : str,
    'qual_ints' : check_np_type('int64'),
    'error' : check_np_type('float64'),
    'description' : str
})
_object = _id


def init(fileh):
    return pd.DataFrame(index=index, columns=columns)

def apply_each(funcs, arg):
    return fzip(funcs, repeat(arg))

def get_row(record):
    #record = next(fileh)
    import sys
    __module__ = sys.modules[__name__]
    clen = len(columns)
    get_getter = compose(attr, "get_{0}".format)
    _getters = map(get_getter, columns)
    self_getters = apply_each(_getters, __module__) #fzip(_getters, repeat(__module__, clen))
    results = apply_each(self_getters, record)
    final_dict = dict(zip(columns, results))
    final_schema.validate(final_dict)
    return final_dict
    #return pd.DataFrame(final_dict, index=index)

#get_fastq_row = compose(get_row, get_fastq)
load_fastq = init_all
class Seq(object):
    get_record = partial(SeqIO.parse, format='fastq')
    #SANGER_OFFSET = 33
    columns = ('id', 'seq', 'quality', 'description', 'qual_ints', 'error')
    SANGER = True
    get_id = attr('id')
    get_seq= compose(str, attr('seq'))
    get_qual_ints = compose(itemgetter('phred_quality'), attr('_per_letter_annotations'))
    get_description = attr('description')
    get_quality = SeqIO.QualityIO._get_sanger_quality_str
    get_error = compose(pmap(error), get_qual_ints)
#    phred_to_char = chr if not SANGER else compose(chr, lambda a: a - 33)
#    get_qual_chars = compose(pmap(phred_to_char), get_qual_ints)
    ''' applies directly to the seqrecord object '''
    index = ['id']

#get_final_dict = compose(self_getters, get_record)
    #TODO: fix for 32-bit systems
    final_schema =  {
        'id' : str,
        'seq' : str,
        'quality' : str,
        'qual_ints' : check_np_type('int64'),
        'error' : check_np_type('float64')
    }
    _object = _id
    '''
    assert len(quality) == len(error) == len(phred_scores)
    '''





#validate = scheme.validate
#TODO: could make these validations match samtools spec
#TODO: Could treat options/cigar string as their own class with their own parsing and validation.


def flatten_vcf(record):
    '''
    :param VCFRecord record: VCFRecord object
    :return dict: all fields as a flat dictionary, including those fields in rec.INFO
    '''
    _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID', 'INFO']
    fields = [getattr(record, _id) for _id in _ids]
    d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
    d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
    return d

def collection_as_df(lambdas, columns, collection):
    '''
    Create a pandas DataFrame by applying a series of functions to a collection.
    :param list lambdas: a list of functions which take exactly one argument (an objects in collection)
    :param list columns: (str) the column names, in order with lambdas
    :param list collection: a list of objects which are valid arguments for the lambdas.
    :return pandas.DataFrame the lambda results on the collection objects as a Matrix.
    '''
    assert len(lambdas) == len(columns), "lambdas must have same length as columns"
    '''use list here to force the evaluation of the functions. otherwise the lambda grabs the last obj evaluated from collection, as in a closure.'''
    values = (list( func(obj) for func in lambdas) for obj in collection)
    return pd.DataFrame(values, columns=columns)

def load_vcf(vcf_records):
    '''
    Convert a list of vcf Records to a pandas DataFrame.
    '''
    return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)

#TODO: properly handle [None], '-', etc.

vcalls = None
ambiguous = None
'''
df[operator.ge(df[field], df[field])]


def df_by_attrs(columns, collection):
    attr_getters = map(attr, columns)
    return collection_as_df(attr_getters, columns, collection)



#    np_type = attr('dtype')
#    pis = partial(partial2, operator.is_)
#    dtype_is = compose(pis, np.dtype)
#    np_type_check = compose(dtype_is, attr('dtype'))
#    compose(attr) foo

#class VCF(object):
    #TODO: uncommnet
#    fields = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID', 'INFO']
#    flat_vcf = flatten_vcf(_self)
#    columns =  flat_vcf.keys()
#    attrs = map(attr, fields)
#    _object = _id



