import os
import hts
import sequtils
import strutils
import tables
import algorithm
import docopt
import kexpr

import ./version

proc bam_filter(argv: var seq[string]): int =
  let args = docopt("""
bam-filter

  Usage: bam-filter [options] <expression> <BAM-or-CRAM>

  -t --threads <threads>       number of BAM decompression threads [default: 0]
  -f --fasta <fasta>           fasta file for use with CRAM files [default: $env_fasta].

valid expressions may contain:

  > mapq/start/pos/end/flag/insert_size (where pos is the 1-based start)

  > is_aligned is_read1 is_read2 is_supplementary is_secondary is_dup is_qcfail
  > is_reverse is_mate_reverse is_pair is_proper_pair is_mate_unmapped is_unmapped

to use aux tags, indicate them prefixed with 'tag_', e.g.:

  tag_NM < 2. Any tag present in the bam can be used in this manner.

example:

    bam-filter "tag_NM == 2 && tag_RG == 'SRR741410' && is_proper_pair" tests/HG02002.bam
  
  """, version=version(), argv=argv)
  #-O --format <output_format>  BAM/CRAM. [default: SAM]
  var
    ex = expression($args["<expression>"])
    bam:Bam
    obam:Bam
    threads = parse_int($args["--threads"])
    fasta: cstring
    format = ""
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  if ex.ke == nil:
    stderr.write_line("[bam-filter] error parsing expression:" & $args["<expression>"])
    stderr.write_line($ex.error())
    quit(2)

  open(bam, $args["<BAM-or-CRAM>"], threads=threads, fai=fasta)
  open(obam, "-", threads=threads, mode="w" & format, fai=fasta)
  obam.write_header(bam.hdr)

  var tags = new_seq[string]()
  ## get the list of tags we need to pull from each record
  if ($args["<expression>"]).contains("tag_"):
    var et = ($args["<expression>"]).split("tag_")
    echo et
    for i in countup(1, len(et) - 1, 1):
      var k: int
      while k < et[i].len and et[i][k].isAlphaAscii:
        k += 1
      assert k == 2
      tags.add(et[i][0..<k])

  for aln in bam:
    ex.clear()
    discard ke_set_int(ex.ke, "mapq", aln.mapping_quality.cint)
    discard ke_set_int(ex.ke, "start", aln.start.cint)
    discard ke_set_int(ex.ke, "pos", (aln.start + 1).cint)
    discard ke_set_int(ex.ke, "flag", (aln.flag).cint)
    discard ke_set_int(ex.ke, "end", aln.stop.cint)
    discard ke_set_int(ex.ke, "insert_size", aln.isize.cint)
    var f = aln.flag
    discard ke_set_int(ex.ke, "is_aligned", if f.unmapped: 0 else: 1)
    discard ke_set_int(ex.ke, "is_unmapped", if f.unmapped: 1 else: 0)
    discard ke_set_int(ex.ke, "is_mate_unmapped", if f.mate_unmapped: 1 else: 0)
    discard ke_set_int(ex.ke, "is_read1", if f.read1: 1 else: 0)
    discard ke_set_int(ex.ke, "is_read2", if f.read2: 1 else: 0)
    discard ke_set_int(ex.ke, "is_supplementary", if f.supplementary: 1 else: 0)
    discard ke_set_int(ex.ke, "is_secondary", if f.secondary: 1 else: 0)
    discard ke_set_int(ex.ke, "is_dup", if f.dup: 1 else: 0)
    discard ke_set_int(ex.ke, "is_reverse", if f.reverse: 1 else: 0)
    discard ke_set_int(ex.ke, "is_mate_reverse", if f.mate_reverse: 1 else: 0)
    discard ke_set_int(ex.ke, "is_pair", if f.pair: 1 else: 0)
    discard ke_set_int(ex.ke, "is_proper_pair", if f.proper_pair: 1 else: 0)

    for itag in tags:
      var io = tag[int](aln, itag)
      if io.isSome:
        discard ke_set_int(ex.ke, "tag_" & itag, io.get)
        continue
      var fo = tag[float](aln, itag)
      if fo.isSome:
        discard ke_set_real(ex.ke, "tag_" & itag, fo.get)
        continue

    if ex.error != 0:
      stderr.write_line "bad"
      quit(2)

    if not (ex.bool):
      if ex.error != 0 and len(tags) == 0:
        stderr.write_line "[bam-filter] expression error:", ex.error
        quit(2)
      continue
    if ex.error != 0 and len(tags) == 0:
        stderr.write_line "[bam-filter] expression error:", ex.error
        quit(2)
    obam.write(aln)

  obam.close()
  return 0
