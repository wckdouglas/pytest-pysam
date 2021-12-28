from collections import OrderedDict

import pysam
from unittest.mock import MagicMock, patch
from bam_filter import filter_short_alignments


class PysamFakeBam:
    def __init__(self, header, reads):
        """
        a mock object that mimics the pysam.AlignmentFile object
        :param pysam.AlignmentHeader header: header of the mock sam file
        :param List[pysam.AlignedSegment] reads: reads of the mock sam file
        """
        self.header = header
        self.reads = reads

    def __iter__(self):
        return iter(self.reads)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


def mock_bam_header(contig_list):
    """
    making a mock pysam.AlignmentHeader object
    Example::
        contigs = [("chr1", 10), ("chr2", 20)]
        mock_header = mock_bam_header(contigs)
    :param List[Tuple[str, int]] contig_list: a list of tuples of (contig name, contig length)
    :return: a pysam.AlignmentHeader object
    :rtype: pysam.AlignmentHeader
    """
    header_dict = OrderedDict(
        [
            ("SQ", [dict(SN=contig[0], LN=contig[1]) for contig in contig_list]),
        ]
    )
    return pysam.AlignmentHeader.from_dict(header_dict)


def mock_alignment(
    header,
    reference_name,
    query_name,
    query_sequence,
    reference_start,
    cigar,
    flag,
    mapping_quality,
    next_reference_name=None,
    next_reference_start=None,
):
    """
    making a mock pysam.AlignedSegment object
    :param pysam.AlignmentHeader header_dict: a pysam alignment header object (can be created by mock_bam_header)
    :param str reference_name: reference name
    :param str query_name: query name
    :param str query_sequence: query sequence
    :param int reference_start: reference start
    :param list cigar: cigar
    :param int flag: flag
    :param int mapping_quality: mapping quality
    :param str next_reference_name: reference name for the paired end alignment mapped
    :param int next_reference_start: reference start of the paired end alignment
    """
    alignment = pysam.AlignedSegment(header)
    alignment.reference_name = reference_name
    alignment.query_name = query_name
    alignment.query_sequence = query_sequence
    alignment.reference_start = reference_start
    alignment.cigar = cigar
    alignment.flag = flag
    alignment.mapping_quality = mapping_quality
    if next_reference_name is not None and next_reference_start is not None and next_reference_start > 0:
        alignment.next_reference_name = next_reference_name
        alignment.next_reference_start = next_reference_start
    return alignment


def test_filter_short_alignments():
    header = mock_bam_header([("chr1", 100)])  # mock a 100 bp chr1 contig
    in_alignment = mock_alignment(
        header=header,
        reference_name="chr1",
        query_name="aln1",
        query_sequence="NNNNNNNNNNN",  # a 11 bp alignment, shouldn't be filtered out
        reference_start=10,
        cigar=[(0, 11)],
        flag=0,
        mapping_quality=30,
    )
    with patch("bam_filter.pysam.AlignmentFile") as pysam_bam:
        mock_in_bam = PysamFakeBam(header, [in_alignment])  # the mock in bam iterarotor will return our mock alignment
        mock_out_bam = MagicMock()
        pysam_bam.return_value.__enter__.side_effect = [
            mock_in_bam,
            mock_out_bam,
        ]  # first call of the pysam.AlignmentFile will return mock_in_bam, second call will be mock_out_bam
        filter_short_alignments(
            "/path/to/inbam", "/path/to/outbam"
        )  # these files are not real, because we are mocking the return of the call anyways

        mock_out_bam.write.assert_called_once_with(
            in_alignment
        )  # because the filter function wouldn't touch alignments with >10 bases
