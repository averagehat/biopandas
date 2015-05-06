import mock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import unittest
from bioframes import bioframes as bf
import sys
import pandas as pd
from pandas.util.testing import assert_series_equal, assert_frame_equal, assert_index_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal

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
    def setUp(self):
        self.r = SeqRecord(Seq("ACGTA"), id="Test", letter_annotations = {"phred_quality":[50, 40, 30, 20, 10]})
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
        self.df = bf.load_fastq(open('tmp.fq'))
#        with mock.patch.object(builtins, 'open', mock.mock_open(read_data=self.fastq_string)): #, create = True) as m:
#            with open('_') as handle:
#                self.df = bf.load_fastq(handle)
        #self.df = mock_file(bf.load_fastq, read_data=self.fastq_string)
        r_dict = bf.get_row(self.r)
        self.series = pd.Series(r_dict)

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
        expected_index = pd.Index(['read1', 'read2', 'read3', 'read4'])
        assert_index_equal(expected_index, self.df.index)



    def test_dataframe_contents(self):
        columns = itemgetter( 'description','seq', 'quality', 'qual_ints', 'error')
        qlen=len( 'TTTCGAATC')
        expected1 = pd.Series(['read1', 'TTTCGAATC', 'FFFFFFFFF', np.array([37]*qlen), np.array( [0.0001]*qlen)])
        expected3 = pd.Series(['read4', 'TTTCGAATC', 'FFFFFFFFF', np.array([37]*qlen), np.array( [0.0001]*qlen)])
        r1, r4, r2 = map(pd.Series, [columns(self.df.ix['read1']), columns(self.df.ix['read4']), columns(self.df.ix['read2'])])
        assert_series_equal( expected1, r1)
        assert_series_equal( expected3, r4)
        expected_error = np.array( [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001])
        expected_qual_ints = np.array( [32, 37, 37, 35, 35, 35, 35, 35])
        expected2 = pd.Series(['read2', 'CTTCGATC', 'AFFDDDDD', expected_qual_ints, expected_error])
        assert_series_equal(expected2, r2)











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

