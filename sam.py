import re
from bioframes import to_np_int
from itertools import groupby
from func import pmap, psplit, pstrip, compose
from operator import itemgetter
from functools import partial
import operator
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
#parse_array = re.compile(r'[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+').match
#parse_array = compose_all(pmap(int), psplit(','), pstrip('[]')) #re.compile(r'[^\[\]]').match)
parse_array = compose_all(to_np_int, psplit(','), pstrip('[]')) #re.compile(r'[^\[\]]').match)
#m = re.compile(r'[^\[\]]+').match
#TODO: ASCII of Phred-scaled base QUALity+33
''' NOTE: samfiles use ASCII of Phred-scaled base QUALity+33 '''
#qual_int = ord

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
key,value = itemgetter(1), itemgetter(0)
groups = groupby(sorted(tups, key=key), key)
get_counts = pmap(compose(int, itemgetter(0)))
sum_counts = compose(sum, get_counts)
cigar_dict = dict( (name, sum_counts(nums)) for name, nums in groups)
mismatches = sum(num for key, num in cigar_dict.items() if key not in 'M=')

#dictmap(compose(sum, get_counts), dict(groups))
#sum(starmap(to_cigar, tups))

#dict(map(reverse, tups))
''' assert sum(itemgetter('M', 'I', 'S', '=', 'X')) == len(seq) == len(quality), \
    "cigar string M/I/S/=/X should sum to the length of the query sequence." '''

#TODO: parse flag
#TODO: handle empty cases (unmapped reads, *)

index = ['QNAME', 'POS', 'REF']
def to_cigar(num, char):
    return cigar_int(char) * int(num)
def cigar_int(char):
    if char == '*': return -1
    if char in ['=', 'M'] : return 0
    else: return 1
#TODO:
'''
POS starts at 1 . . . but if POS is set as 0 for an unmapped read without coordinate. If POS is 0, no assumptions can be made about RNAME and CIGAR
Bit 0x4 is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no
assumptions can be made about RNAME , POS , CIGAR , MAPQ , and bits 0x2, 0x100, and 0x800 .'''

eval_flag = compose(bool, operator.and_)
flag_results = dict((meaning, eval_flag(bit, flag)) for bit, meaning in flag_meanings.items())

#TODO: parse quality
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


sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPTIONS"]
#TODO: get_record function takes a filehandle and returns a single record via SeqIO, etc.
#basics, options = fields[:11], fields[11:]
#So functions expect a dictionary I guess
readline = compose(psplit('\t'), next)
line_to_dict = compose(dict, partial(zip, sam_columns))
get_record = compose(line_to_dict, readline)
#fields = readline()
#basics, options = fields[:len(sam_columns)], fields[len(sam_columns):]

scheme={
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

