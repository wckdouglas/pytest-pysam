from pathlib import Path

import pysam


def filter_short_alignments(in_bam_file: Path, out_bam_file: Path):
    """
    reading the input bam file and write alignments with >10 bases to output bam file

    :param str in_bam_file: file path to input bam file
    :param str out_bam_file: file path to output bam file
    """
    with pysam.AlignmentFile(in_bam_file) as inbam:
        with pysam.AlignmentFile(out_bam_file, "wb", template=inbam) as outbam:
            for aln in inbam:
                if len(aln.query_sequence) > 10:
                    outbam.write(aln)
