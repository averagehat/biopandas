from biotest import BioTest
from nfilter import make_filtered, write_filtered, fqs_excluding_indices, write_post_filter
import os
import mock
#
class TestNGSFilter(BioTest):

    def setUp(self):
        self.actualfn = 'tests/testoutput/filtered.1900_S118_L001_R2_001_2015_04_24.fastq'
        self.expectedfn = 'tests/expected__R2__.fastq'
        self.inputfn = 'tests/1900/1900_S118_L001_R2_001_2015_04_24.fastq'
        self.inputdir = 'tests/1900'
        self.outdir = 'tests/testoutput'

    def tearDown(self):
        if os.path.exists(self.actualfn):
            os.remove(self.actualfn)

    def test_write_filtered_unfiltered(self):
        ''' input and output fastq files should be the same if filters make everything pass '''
        write_filtered(self.inputfn, 0, False, outdir=self.outdir)
        self.seqs_equal(self.inputfn, self.actualfn, format='fastq')

    def test_filter_pool(self):
        write_post_filter(self.inputdir, 32, True, self.outdir)
        actual = open(self.actualfn)
        expected = open(self.expectedfn)
        self.assertFilesEqual(expected, actual)

    @mock.patch('nfilter.os.listdir')
    def test_fqs_excluding_indices_extensions(self, mlistdir):
        mlistdir.return_value = candidates = ['foo.fq', 'foo.sff', 'foo.fastq', 'foo.bad']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo.fq', 'D/foo.sff', 'D/foo.fastq']
        self.assertListEqual(expected, actual)

    @mock.patch('nfilter.os.listdir')
    def test_fqs_excluding_indices_excludes_index(self, mlistdir):
        mlistdir.return_value = ['foo_1R_.fq', 'foo.sff', 'foo_I2_.fastq', 'foo__I1__.fq']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo_1R_.fq', 'D/foo.sff' ]
        self.assertListEqual(expected, actual)

    def test_filter_raises_error_on_empty_filtered_result(self):
        ''' This should raise an AssertionError because no reads will be left after that quality filter.'''
        with self.assertRaises(AssertionError):
            write_filtered(self.inputfn, 65, True)


