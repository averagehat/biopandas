from functools import partial
import pandas as pd
from past.builtins import map , filter
from operator import attrgetter as attr
from itertools import starmap

def _id(x): return x

def apply_key_func(k, v, funcdict):
    return funcdict.get(k, _id)(v)

def compose_all(*funcs):
    return reduce(compose, funcs)

def k_compose(outer, **okwargs):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **ikwargs):
        _kwargs = dict( (k, apply_key_func(k, v, okwargs)) for k, v in ikwargs.items())
        return outer(*args, **_kwargs)
    return newfunc

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc

def fzip(funcs, args):
    for func, arg in izip(funcs, args):
        yield func(arg)

class Seq(object)
    #SANGER_OFFSET = 33
    SANGER = True
    get_seq_str = compose(str, attr('seq'))
    get_id = attr('id')
    get_qual_ints = compose(itemgetter('phred_quality'), attr('_per_letter_annotations'))
    get_description = attr('description')
    phred_to_char = chr if not SANGER else compose(chr, lambda a: a - 33)
    get_qual_chars = compose(pmap(phred_to_char), get_qual_ints)
    index = 'id'
    _object = _id



class VCF(object):
#    fields = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID', 'INFO']
    flat_vcf = flatten_vcf(_self)
    columns =  flat_vcf.keys()
    attrs = map(attr, fields)
    _object = _id



sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]

{
    'QNAME' : str,
    'FLAG' : int,
    'RNAME' : str,
    'POS' : int,
    'MAPQ' : int,
    'CIGAR' : str,
    'MRNM' : '*='.__contains__,
    'MPOS' : int,
    'ISIZE' : int,
    'SEQ' : str,
    'QUAL' : str,
    'OPT' : str
}




def samview_to_df(rawtext):
    samtext = '\n'.join( fixline(row) for row in rawtext.split('\n') )
    as_bytes = BytesIO(samtext)
    #Write test to ensure handles varaibles #columns
    return pd.read_csv(as_bytes, names=sam_columns, usecols=sam_columns, delimiter='\t', header=None, squeeze=True)

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



def df_by_attrs(columns, collection):
    attr_getters = map(attr, columns)
    return collection_as_df(attr_getters, columns, collection)

def reverse(collection): return collection[::-1]

# Parse options
from pyparsing import Regex
#_join = partial(reduce, lambda a, b: a+':'+b)
tag = Regex(r'[A-Za-z][A-Za-z0-9]')
_type = Regex(r'[AifZHB]')
value = Regex('[^\t]')
full = tag + ':' + _type + ':' + value
#reduce(operator.add, [tag, _type, value], ':')
#cigar_regex = r'\*|([0-9]+[MIDNSHPX=])+'
#? makes the regex not be too greedy

full.parseString('AS:i:213')
#full = _join( [tag, _type, value ] )
#compose3 = partial(reduce, compose)
compose_list = partial(reduce, compose)
compose_all = compose(compose_list, lambda *a : a)
pmap = partial(partial, map)
pstrip = partial(partial, string.strip)
pfilter = partial(partial, filter)
pstrip = lambda x: partial(string.split, chars=x)
psplit = lambda x: partial(string.split, sep=x)
#parse_array = re.compile(r'[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+').match
parse_array = compose_all(pmap(int), psplit(','), pstrip('[]')) #re.compile(r'[^\[\]]').match)
#m = re.compile(r'[^\[\]]+').match
{
    'A' : chr,
    'i' : int,
    'f' : float,
    'Z' : str,
    'H' : int, # hex
    'B' : parse_array
}

#parse cigar string
cigar_regex = r'(?:([0-9]+)([MIDNSHPX=]))+?'
reg = re.compile(cigar_regex)
tups = reg.findall('15S213M23S')
sum(starmap(to_cigar, tups))
#dict(map(reverse, tups))
def to_cigar(num, char):
    return cigar_int(char) * int(num)
def cigar_int(char):
    if char == '*': return -1
    if char == '=' : return 0
    else: return 1


#TODO: parse quality
