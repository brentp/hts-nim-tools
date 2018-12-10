import os
import random
import hts
import docopt
import lapper
import strutils
import tables
import algorithm

import ./version

type
  cnv = ref object
    chrom: string
    start: int
    stop: int
    prob: float64

proc start(r: cnv): int {.inline.} = return r.start
proc stop(r: cnv): int {.inline.} = return r.stop

proc bed_line_to_cnv(line: string): cnv =
  var
   cse = line.strip().split('\t', 7)

  if len(cse) < 3:
    stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
    return nil
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    prob:float64

  if len(cse) >= 4:
    prob = parse_float(cse[3])
  else:
    prob = 0.5

  return cnv(chrom: cse[0], start: s, stop: e, prob:prob)

proc bed_to_cnv_table(bed: string): TableRef[string, seq[cnv]] =
  var bed_regions = newTable[string, seq[cnv]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr = kstring_t(l:0, m:0, s:nil)
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_cnv($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[cnv]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: cnv): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc overlap_p(a: cnv, b:Record): float =
  var o = min(a.stop, b.stop) - max(a.start, b.start)
  if o < 0: return 0
  return o.float / float(b.stop - a.start)

proc internal_sampler(ibam:Bam, obam:var Bam, regions:TableRef[string, seq[cnv]]) =
  var last_chrom: string
  var lap:Lapper[cnv]
  var res = new_seq[cnv](20)
  randomize()
  for record in ibam:
    if not regions.contains(record.chrom):
      obam.write(record)
      continue
    if record.chrom != last_chrom:
      lap = lapify(regions[record.chrom])
      last_chrom = record.chrom

    discard lap.find(record.start, record.stop, res)
    if len(res) == 0:
        obam.write(record)
        continue

    var po = overlap_p(res[0], record)

    if res[0].prob < 1 and random(1.0) < res[0].prob:
      if po == 1 or po < 3.0 * random(1.0):
        obam.write(record)
      continue

    # TODO: allow writing single-end reads.
    #if res[0].prob > 1:
    #  obam.write_record(record)
    #  if random(1.0) > 1 / prob:
    #    obam.write_record(record)


proc copy_number_sampler(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  $version

  Usage: copy-number-sampler [options] <BED> <BAM-or-CRAM>

BED format looks like:

  chr1\t1123\t345\t0.5

where the final column indicates the sampling probability.

Arguments:                                                                                                                                                 

  <BED>          the bed file containing regions in which to sample reads; 4th column is a float indicating sampling probability.
  <BAM-or-CRAM>  the alignment file for which to calculate depth.

Options:

  -f --fasta <fasta>          fasta file for use with CRAM files [default: $env_fasta].
  -h --help                     show help
  """ % ["version", version(), "env_fasta", env_fasta])

  let args = docopt(doc, version=version(), argv=argv)

  var fasta: cstring 
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var bam:Bam
  open(bam, $args["<BAM-or-CRAM>"], fai=fasta, threads=1)

  var obam:Bam
  open(obam, "sampled.bam", fai=fasta, mode="wb", threads=1)

  obam.write_header(bam.hdr)

  var regions = bed_to_cnv_table($args["<BED>"])
  internal_sampler(bam, obam, regions)
  obam.close()
  return 0
