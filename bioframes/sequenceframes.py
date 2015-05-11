from __future__ import print_function
from operator import attrgetter as attr
from operator import  itemgetter
from Bio import SeqIO
from func import compose, compose_all,  imap, apply_each, _id, get_funcs
from bioframes import error, get_fastq, check_np_type
import pandas as pd
import numpy as np
from schema import Schema
from collections import namedtuple
from functools import partial
''' All the boring stuff to interact with SeqRecord API '''
#TODO: fix for 32-bit systems



class BioFrame(pd.DataFrame):
    def find_duplicates(frame):
        ''' not sure why this one works. '''
        ''' for example, find pcr duplicates by finding duplicate reads'''
        #df2[df2['b'] == df2['b'] & (df2['a'] == df2['a'])]
        ''' duplicated by defaults returns all but the first occurrence; take_last takes only the first occurence,
        so unioning them (via | (or)), returns all duplicates. accepts any number of keys. '''
        frame[frame.duplicated(['b'],take_last=True) | frame.duplicated(['b'])]



def fqframe(fileh):
    final_schema =  Schema({
        'id' : str,
        'seq' : str,
        'quality' : str,
        'qual_ints' : check_np_type('int64'),
        'error' : check_np_type('float64'),
        'description' : str
    })

    columns = ('id', 'seq', 'quality', 'description', 'qual_ints', 'error')
    SANGER = True
    get_id = attr('id')
    get_seq= compose(str, attr('seq'))
    get_qual_ints = compose_all(np.array, itemgetter('phred_quality'), attr('_per_letter_annotations'))
    get_description = attr('description')
    get_quality = SeqIO.QualityIO._get_sanger_quality_str
    get_error = compose(error, get_qual_ints)
    getters = [get_id, get_seq, get_quality, get_description, get_qual_ints, get_error]
    metadata = {'filename': fileh.name}
    iterator = get_fastq(fileh)
    get_raw_record = partial(next, iterator)

    return { 'obj_func' : get_raw_record,
        'columns' : columns,
        'getters' : getters,
        'validator' : final_schema,
        'dictgetters' : None
    }


