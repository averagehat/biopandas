

class TestClassicSam(unittest.TestCase):
    def setUp(self):
        self.samtext	=	'\n'.join(['read1\t1\tchr1\t1\t60	10M	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
                                  'read2	1	chr2	1	60	10M	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
                                  'read3	1	chr1	1	60	10M	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq'])
    def test_sam_to_df(self):
        result = pc.samview_to_df(self.samtext)
        self.assertEquals(result.columns.tolist(), self.columns)
        self.assertEquals(result.ix[2]['QNAME'], 'read3')

    def test_df_from_collection_attributes(self):
        mocks = [mock.Mock() for i in range(5)]
        [[mock_do() for i in range(index)] for index, mock_do in enumerate(mocks)]
        columns = ['call_count', 'called']
        expected = pd.DataFrame( [(i, bool(i)) for i in range(5)])
        expected.columns = columns
        result = a2f.df_from_collection_attributes(columns, mocks)
        assert_frame_equal(expected, result)


class TestSamFile(unittest.TestCase):
    def setUp(self):
        #TODO: fix this so data is accurate;  i.e.:
   # ''' assert sum(itemgetter('M', 'I', 'S', '=', 'X')) == len(seq) == len(quality), \ "cigar string M/I/S/=/X should sum to the length of the query sequence." '''

        self.samtext	=	'\n'.join(['read1\t1\tchr1\t1\t60	10M 	=	1	1	TTTCGAATC	FFFFFFFFF	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq',
                                  'read2	1	chr2	1	8	3I3M4D	=	1	1	CTTCGATC	AFFDDDDD	NM:i:3	AS:i:2	XS:i:0	RG:Z:Sanger',
                                  'read3	1	chr1	1	60	2D2M2I2M2=	=	1	1	CCGATCAA	FF@@@F@F	NM:i:3	AS:i:231	XS:i:0	RG:Z:MiSeq'])

    def test_cigar_scores(self):
        #TODO: with mock_open
        e_strings = pd.Series(['10M' '3I3M4D' '2D2M2I2M2='])
        e_m, e_i, e_d, e_total = pd.Series([10, 3, 4]), pd.Series([0, 3, 2]), pd.Series([0, 4, 2]), pd.Series([0, 7, 4])
        self.assertEquals(result.cigar_i, e_i)
        self.assertEquals(result.cigar_m, e_m)
        self.assertEquals(result.cigar_d, e_d)
        self.assertEquals(result.cigar_total, e_total)
        self.assertEquals(result.cigar_string, e_strings)

    def test_options(self):
       e_nm, e_as, e_xs, e_rg = map(pd.Series, [[3, 3, 3], [231, 2, 231], [0, 0, 0], ['MiSeq', 'Sanger', 'MiSeq']])
       results = [df.NM, df.AS, df.XS, df.RG]
       map(self.assertEquals, [e_nm, e_as, e_xs, e_rg], results)

    def test_df(self):
        expected = ['read3',	1,	'chr1',	1,	60,	'2D2M2I2M2=',	'=',	1,	1,	'CCGATCAA',	'FF@@@F@F',	3, 231, 0, 'MiSeq']#, 4, 2, 2, 4
        self.assertEquals(df.ix['read3', 1, 'chr1'][:len(expected)], expected)

    def test_flags(self):
        middle_e = [False, False, False, True]
        outer_e = [False, False, True, True]
        names = ("template having multiple segments in sequencing",
        "each segment properly aligned according to the aligner",
        "segment unmapped",
        "next segment in the template unmapped")
        inner_actual = df['read2'][names]
        outer_actual1= df['read1'][names]
        outer_actual3= df['read3'][names]
        self.assertEquals(middle_e, inner_actual)
        self.assertEquals(outer_e, outer_actual1)
        self.assertEquals(outer_e, outer_actual3)



