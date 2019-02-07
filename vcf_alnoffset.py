#!/usr/bin/env python
"""Convert vcf positions to alignment positions or unaligned positions
of a target sequence

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
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
import difflib

#--- third-party imports
#
# /

#--- project specific imports
#
import posmap
import fasta_parser
import vcf


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"


# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      dest="verbose",
                      action="store_true",
                      help=SUPPRESS_HELP) #"be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true",
                      help=SUPPRESS_HELP) #"debugging")
    parser.add_option("-a", "--msa",
                      dest="msa_file",
                      help="Multiple sequence alignment")
    parser.add_option("-i", "--vcf-in",
                      help="VCF input file")
    parser.add_option("-o", "--vcf-out",
                      default='-',
                      help="VCF output file ('- for stdout = default)'")
    parser.add_option("-s", "--seqid",
                      dest="seq_id",
                      help="Sequence id of interest (for which VCF"
                      " file was produced). If empty aligned positions are assumed!")
    parser.add_option("-m", "--map-to-id",
                      dest="map_to_id",
                      help="Convert SNP positions to unaligned"
                      " positions of this seq instead of alignment"
                      " coordinates")
    parser.add_option("", "--force-overwrite",
                      action="store_true",
                      dest="force_overwrite",
                      help="Force overwriting of output file")
    default = "fasta"
    parser.add_option("", "--fmt",
                      dest="aln_fmt",
                      default=default,
                      help="Alignment format (must be supported by"
                      " Biopython; default is '%s')" % default)
    return parser



def map_pos_of_vars(var_list_in, seq_id, map_to_id, pos_map):
    """The actual main function which will convert snp positions of
    snps in var_list_in coming from seq_id to map_to_id (or alignment
    if None) based on PosMap pos_map. If seq_id is empty positions are
    assumed to be aligned
    """

    assert seq_id or map_to_id, ("Need src or dst id")

    conv_pos_map = pos_map.convert(seq_id, map_to_id)

    var_list_offset = []
    for s in var_list_in:
        try:
            LOG.debug("Changing pos from %d to %d for %s",
                      s.POS, conv_pos_map[s.POS], s)
            s.INFO['orig-pos'] = s.POS
            #import pdb; pdb.set_trace()
            s.POS = conv_pos_map[s.POS]
        except KeyError:
            if map_to_id:
                LOG.warn("Position %d in %s has no equivalent in %s.",
                         s.POS, seq_id, map_to_id)
                continue
            else:
                LOG.fatal("INTERNAL ERROR: Position %d in %s seems invalid.",
                          s.POS, seq_id)
                raise

        if s.POS == -1:
            LOG.warn("No alignment match for this SNP. Skipping: %s", s)
            continue

        #LOG.warn("Need to append to INFO")
        if map_to_id:
            s.INFO['offset'] = map_to_id
        else:
            s.INFO['offset'] = "aligned"
        var_list_offset.append(s)

        #print "old", s
        #print "new", s
        #print

    return var_list_offset



def main():
    """
    The main function
    """


    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)


    if len(args) != 0:
        parser.error("Unrecognized args found")
        sys.exit(1)

    # file check
    for (filename, descr, direction, mandatory) in [
            (opts.vcf_in, "VCF input file", 'in', True),
            (opts.vcf_out, "VCF output file", 'out', True),
            (opts.msa_file, "MSA file", 'in', True)]:

        if not mandatory and not filename:
            continue

        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if filename == '-':
            continue

        if direction == 'in' and not os.path.exists(filename):
            LOG.fatal(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)
        if direction == 'out' and os.path.exists(filename) \
          and not opts.force_overwrite:
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)

    msa_seqs = fasta_parser.to_dict(
        fasta_parser.fasta_iter(opts.msa_file))
    if len(msa_seqs) < 2:
        LOG.fatal("%s contains only one sequence\n" % (opts.msa_file))
        sys.exit(1)



    seq_ids_to_check = []
    if opts.seq_id:
        seq_ids_to_check.append(opts.seq_id)
    if opts.map_to_id:
        seq_ids_to_check.append(opts.map_to_id)
    for seq_id in seq_ids_to_check:
        if not seq_id in msa_seqs.keys():
            bestmatch = None
            matches = difflib.get_close_matches(
                seq_id, msa_seqs.keys(), 1, 0.5)
            if matches:
                bestmatch = matches[0]
            LOG.fatal("Couldn't find seq '%s' in MSA file"
                      " (best match was: %s)" % (
                          seq_id, bestmatch))
            sys.exit(1)

    pos_map = posmap.PosMap(msa_seqs.values())
    if opts.debug:
        pos_map.output()
    vcf_reader = vcf.VCFReader(filename=opts.vcf_in)
    var_list = [v for v in vcf_reader]
    if not all([v.is_snp for v in var_list]):
        LOG.warn("Found non-SNPs in vcf file (e.g. indels). Will ignore them")
        var_list = [v for v in var_list if v.is_snp]
           
    var_list_offset = map_pos_of_vars(
        var_list, opts.seq_id, opts.map_to_id, pos_map)

    if opts.vcf_out == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(opts.vcf_out, 'w')
    vcf_writer = vcf.VCFWriter(fh_out, vcf_reader)
    for var in  var_list_offset:
        vcf_writer.write_record(var)
    if fh_out != sys.stdout:
        fh_out.close()


if __name__ == "__main__":
    main()
    LOG.debug("Successful exit")
    # FIXME in all other scripts LOG.info works but here it's always
    # printed. Couldn't find out why
    #LOG.critical("v=%s d=%s" % (opts.verbose, opts.debug))
    #LOG.critical("LOG.isEnabledFor(logging.INFO)=%s" % LOG.isEnabledFor(logging.INFO))
    #LOG.critical("LOG.getEffectiveLevel()=%s" % LOG.getEffectiveLevel())
