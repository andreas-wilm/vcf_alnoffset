#!/usr/bin/env python
"""A Fasta parser
"""

# Copyright (C) 2011-2016 Genome Institute of Singapore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.



#--- standard library imports
#
from itertools import groupby

#--- third-party imports
#
# /

#--- project specific imports
#
#/

class SeqRecord(object):
    """Fake, poor man's version of a simple Biopython's seqrecord for fasta files
    """

    def __init__(self, header=None, seq=None):
        """
        """

        assert not header.startswith('>')
        self.id = header.split(' ')[0]
        self.name = header.split(' ')[0]
        self.description = header
        self.seq = seq


def to_dict(seqrecs):
    """Similar to Biopython's SeqIO.to_dict
    """
    return dict([(s.id, s) for s in seqrecs])


def fasta_iter(fasta_name):
    """
    Given a fasta file, yield a fake SeqRecord

    Original author: Brent Pedersen
    http://www.biostars.org/p/710/
    """

    with open(fasta_name) as fh:

        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            #yield header, seq
            yield SeqRecord(header, seq)



if __name__ == "__main__":
    import doctest
    doctest.testmod()
