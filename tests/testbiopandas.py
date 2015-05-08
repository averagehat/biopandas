import mock
from operator import attrgetter as attr, methodcaller
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import unittest
import sys
import pandas as pd
from pandas.util.testing import assert_series_equal, assert_frame_equal, assert_index_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal
import bioframes
import sequenceframes
from functools import partial
#from bioframes.bioframes import makeframe
import numpy as np
from operator import itemgetter
if sys.version[0] == '2':
        import __builtin__ as builtins  # pylint:disable=import-error
else:
        import builtins  # pylint:disable=import-error

def mock_file(func, read_data, *args, **kwargs):
    with mock.patch.object(builtins, 'open', mock.mock_open(read_data=read_data)): #, create = True) as m:
        with open('_') as handle:
            return func(handle, *args, **kwargs)


class TestFastq(unittest.TestCase):
    #@mock.patch('bioframes.get_fastq')
    def setUp(self):
        self.fastq_string = '''@read1
TTTCGAATC
+
FFFFFFFFF
@read2
CTTCGATC
+
AFFDDDDD
@read3
CCGATCAA
+
FF@@@F@F
@read4
TTTCGAATC
+
FFFFFFFFF
'''
        with open('tmp.fq', 'w') as tmp: tmp.write(self.fastq_string)

        fq = sequenceframes.fqframe(open('tmp.fq'))
        #self.df = sequenceframes.load_fastq(open('tmp.fq'))
        self.df = bioframes.makeframe(fq)
        self.r =SeqRecord(Seq("ACGTA"), id="Test", letter_annotations = {"phred_quality":[50, 40, 30, 20, 10]})
        with open('tmp.fq', 'w') as tmp: SeqIO.write([self.r], tmp, 'fastq')
        fq2 = sequenceframes.fqframe(open('tmp.fq'))
        df =  bioframes.makeframe(fq2)
        self.series = df.ix[0]
        #TODO: somehow SeqIO broke when I tried to mock_open
#        with mock.patch.object(builtins, 'open', mock.mock_open(read_data=self.fastq_string)): #, create = True) as m:
#            with open('_') as handle:
#                self.df = bf.load_fastq(handle)
        #self.df = mock_file(bf.load_fastq, read_data=self.fastq_string)

    def test_sanger_quality_error(self):
        expected = np.array([.1,  .01, .001, .0001, .00001][::-1])
        assert_array_almost_equal(self.series['error'], expected)

    def test_sanger_quality_string(self):
        self.assertEquals(self.series['quality'], 'SI?5+')

    def test_data_frame_lengths(self):
        expected_len = len(self.fastq_string.split('\n')) / 4
        self.assertEquals( expected_len, len(self.df))

    def test_dupe_reads(self):
        dupe_reads = self.df[self.df.seq == 'TTTCGAATC']
        self.assertEquals(2, len(dupe_reads))

    def test_dataframe_index(self):
        #expected_index = pd.Index(['read1', 'read2', 'read3', 'read4'])
        expected_index = pd.Index([0, 1, 2, 3])
        assert_index_equal(expected_index, self.df.index)

    def gen_assert(self, left, right):
        if hasattr(left, 'dtype'):
            return assert_array_almost_equal(left, right)
        else:
            return self.assertEquals(left, right)

#TODO: fix this test
    def test_dataframe_contents(self):
        from numpy import array
        q1, q2, q3, q4 = (array([ 0.00019953,  0.00019953,  0.00019953,  0.00019953,  0.00019953, 0.00019953,  0.00019953,  0.00019953,  0.00019953]), \
                          array([ 0.00063096,  0.00019953,  0.00019953,  0.00031623,  0.00031623, 0.00031623,  0.00031623,  0.00031623]), \
                          array([ 0.00019953,  0.00019953,  0.00079433,  0.00079433,  0.00079433, 0.00019953,  0.00079433,  0.00019953]), \
                          array([ 0.00019953,  0.00019953,  0.00019953,  0.00019953,  0.00019953, 0.00019953,  0.00019953,  0.00019953,  0.00019953]))

        columns = itemgetter( 'description','seq', 'quality', 'qual_ints', 'error')
        qlen=len( 'TTTCGAATC')
        expected1 = itemgetter(0, -1, 0, -2, 2, 1)(['read1', 'TTTCGAATC', 'FFFFFFFFF', np.repeat(37, qlen), q1 ])
        expected3 = itemgetter(0, -1, 0, -2, 2, 1)(['read4', 'TTTCGAATC', 'FFFFFFFFF', np.repeat(37, qlen), q4])
        #r1, r4, r2 = map(attr('values'), (columns(self.df[self.df['id']=='read1']), columns(self.df[self.df['id'] == 'read4']), columns(self.df[self.df['id'] == 'read2'])))
        r1, r4, r2 = itemgetter(0, 3, 1)(self.df.ix)
        r1, r4, r2 = map(methodcaller('tolist'), map(attr('values'), (r1, r4, r2)))
        map(self.gen_assert,  expected1, r1)
        map(self.gen_assert,expected3, r4)
        expected_qual_ints = np.array( [32, 37, 37, 35, 35, 35, 35, 35])
        expected2 = pd.Series(  itemgetter(0, -1, 0, -2, 2, 1)(['read2', 'CTTCGATC', 'AFFDDDDD', expected_qual_ints, q2]))
        map(self.gen_assert, expected2, r2)











#    def test_join_non_unique_dataframes(self):
#        '''
#        df1 and df2 share an index with duplicates, check that it is aligned correctly
#        '''
#        rows1 = [('A', 'A1'), ('B', 'B1'), ('A', 'A2'), ('C', 'C1')]
#        rows2 = [('A', '0A', False), ('B', '0B', True), ('A', '00A', False), ('C', '00C', True)]
#        self.df1, self.df2 = map(make_df_index0, (rows1, rows2))
#        self.df1.columns = ['0', '1']
#        self.df2.columns = ['0', '1', '2']
#        self.df1, self.df2 = self.df1.set_index('0'), self.df2.set_index('0')
#        result = a2f.join_non_unique_dataframes(self.df1, self.df2)
#        expected = pd.DataFrame(
#                [('A', 0, 'A1', '0A', True), ('B', 0, 'B1', '0B', True),
#                ('A', 1, 'A2', '00A', False), ('C', 0, 'C1', '00C', True)]
#                ).set_index(0).set_index(1, append=True)
#        assert_frame_equal(result, expected)

