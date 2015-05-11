import re
import pandas as pd
from bioframes import to_np_int, sanger_qual_str_to_error, qual_to_phreds, makeframe
from itertools import groupby
from func import pmap, psplit, pstrip, compose, compose_all, merge_dicts, fzip, \
    partial2, dictmap, starcompose, split_list, nameddict
from operator import itemgetter
from functools import partial
import operator as op
from operator import attrgetter as attr
from schema import Schema, Use
from itertools import ifilter, imap
import samtools
#from pyparsing import Regex

#TODO: don't compute quality-ints twice.

''' NOTE: samfiles use ASCII of Phred-scaled base QUALity+33 '''
parse_array = compose_all(to_np_int, psplit(','), pstrip('[]'))
tabsplit = psplit('\t')

basic_scheme={
    'QNAME' : str,
    'FLAG' : int,
    'RNAME' : str,
    'POS' : int,
    'MAPQ' : int,
    'CIGAR' : str,
    'RNEXT' : str,
    'PNEXT' : int,
    'TLEN' : int,
    'SEQ' : str,
    'QUAL' : str,
}

basic_schema = Schema(dictmap(Use, basic_scheme))

options_scheme = {
    'A' : chr,
    'i' : int,
    'f' : float,
    'Z' : str,
    'H' : int, # hex
    'B' : parse_array
}

def parse_option(option_str):
#    _tag = Regex(r'[A-Za-z][A-Za-z0-9]')
    tag, _type, raw_val = psplit(':')(option_str)
    val = options_scheme[_type](raw_val)
    return tag, val


#TODO: handle empty cases (unmapped reads, *)

def parse_cigar(cigar_str):
    #? makes the regex not be too greedy
    cigar_regex = r'(?:([0-9]+)([MIDNSHPX=]))+?'
    reg = re.compile(cigar_regex)
    tups = reg.findall(cigar_str)
    key, value = itemgetter(1), itemgetter(0)
    groups = groupby(sorted(tups, key=key), key)
    get_counts = pmap(compose(int, value))
    sum_counts = compose(sum, get_counts)
    s = "cigar_{0}".format
    cigar_dict = dict( (s(name), sum_counts(nums)) for name, nums in groups)
    mismatches = sum(num for k, num in cigar_dict.items() if k not in ['cigar_M', 'cigar_='])
    return merge_dicts(cigar_dict, {'cigar_score': mismatches})

''' assert sum(itemgetter('M', 'I', 'S', '=', 'X')) == len(seq) == len(quality), \
    "cigar string M/I/S/=/X should sum to the length of the query sequence." '''

index = ['QNAME', 'POS', 'RNAME']
#TODO: filter unmapped reads
''' POS starts at 1 . . . but if POS is set as 0 for an unmapped read without coordinate. If POS is 0, no assumptions can be made about RNAME and CIGAR
Bit 0x4 is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no
assumptions can be made about RNAME , POS , CIGAR , MAPQ , and bits 0x2, 0x100, and 0x800 .'''


flag_meanings = {
0x1  : "template having multiple segments in sequencing",
0x2  : "each segment properly aligned according to the aligner",
0x4  : "segment unmapped",
0x8  : "next segment in the template unmapped",
0x10 : "SEQ being reverse complemented",
0x20 : "SEQ of the next segment in the template being reversed",
0x40 : "the rst segment in the template",
0x80 : "the last segment in the template",
0x100: "secondary alignment",
0x200: "not passing quality controls",
0x400: "PCR or optical duplicate",
0x800: "supplementary alignment"
}


eval_flag = compose(bool, op.and_)

def flag_dict(flag):
    return dict((meaning, eval_flag(bit, flag)) for bit, meaning in flag_meanings.items())

sam_columns = ("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL") #optiosn


parse_options = compose(dict, pmap(parse_option)) #, tabsplit)
#readfields = compose(tabsplit, next)
line_to_dict = compose_all(dict, partial(zip, sam_columns)) #, tabsplit)
validated_dict = compose(basic_schema.validate, line_to_dict)
fields_and_options = compose(partial2(split_list, len(sam_columns)), tabsplit)
parsers = partial(fzip, [validated_dict, parse_options])
parse_fields_and_options = compose(parsers, fields_and_options)
all_but_cigar_dict = starcompose(merge_dicts, parse_fields_and_options)
get_cigar_dict = compose(parse_cigar, itemgetter('CIGAR'))
get_flag_dict = compose(flag_dict, itemgetter('FLAG'))
get_error = compose(sanger_qual_str_to_error, itemgetter('QUAL'))

def load_sam(fh):
    dicts = map(get_row, ifilter(bool, fh.read().split('\n')))
    return pd.DataFrame(dicts)
# .set_index(index) #, index=index, columns=columns)


def get_row(row):
    result = all_but_cigar_dict(row)
    return merge_dicts(result, get_cigar_dict(result), get_flag_dict(result), {'error' :  get_error(result)})


def pileframe(fh):
    obj_attributes = ('ref','pos','refbase','depth','_bases','_bquals')# + ('error', 'qual_ints')#,'_mquals = parts
    columns =        ('chrom', 'pos', 'refbase', 'depth', 'bases', 'quals')# + ('error', 'qual_ints') #optionall bquals
    getters = map(attr, obj_attributes)
    iter_rows = imap(samtools.MPileupColumn, fh)
    dictgetters =  {
        'error' : compose(sanger_qual_str_to_error, itemgetter('quals')),
        'qual_ints' : compose(qual_to_phreds, itemgetter('quals'))
    }


  #  def initialize():
  #      return collection_as_df(list(getters), columns, list(iter_rows))

    #df = initialize()
    return {
        'obj_func' : partial(next, iter_rows),
        'columns' : columns,
        'getters' : getters,
        'validator' : None,
        'dictgetters' : dictgetters
    }

load_pileup = compose(makeframe, pileframe)
    #return nameddict('PileupFrame',  { 'df' : df })

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
















