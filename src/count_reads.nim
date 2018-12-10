import os
import hts
import docopt
import lapper
import strutils
import tables
import algorithm

import ./version

type
  region_t = ref object
    chrom: string
    start: int
    stop: int
    name: string
    count: int

proc inc_count(r:region_t) = inc(r.count)
proc start(r: region_t): int {.inline.} = return r.start
proc stop(r: region_t): int {.inline.} = return r.stop
proc tostring(r: region_t, s:var string) {.inline.} =
  s.set_len(0)
  s.add(r.chrom & "\t" & $r.start & "\t" & $r.stop & "\t")
  if r.name != "":
    s.add(r.name & "\t")
  s.add($r.count)

proc bed_line_to_region(line: string): region_t =
  var
   cse = line.strip().split('\t', 5)

  if len(cse) < 3:
    stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
    return nil
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
  if len(cse) > 3:
   reg.name = cse[3]
  return reg

proc bed_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc internal_count(bam:Bam, mapq:uint8, eflag:uint16, regions:TableRef[string, seq[region_t]]) =
  for chrom in regions.keys():
    if not regions.contains(chrom) or regions[chrom].len == 0:
      continue
    var lap:Lapper[region_t] = lapify(regions[chrom])
    for aln in bam.querys(chrom):
      if aln.mapping_quality < mapq: continue
      if (aln.flag and eflag) != 0: continue

      lap.each_seek(aln.start, aln.stop, inc_count)
    var s = new_string_of_cap(1000)
    for region in regions[chrom]:
      region.tostring(s)
      echo s

proc count_reads(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  $version

  Usage: count-reads [options] <BED> <BAM-or-CRAM>

Arguments:                                                                                                                                                 

  <BED>          the bed file containing regions in which to count reads.
  <BAM-or-CRAM>  the alignment file for which to calculate depth.

Options:

  -t --threads <threads>      number of BAM decompression threads [default: 0]
  -f --fasta <fasta>          fasta file for use with CRAM files [default: $env_fasta].
  -F --flag <FLAG>            exclude reads with any of the bits in FLAG set [default: 1796]
  -Q --mapq <mapq>            mapping quality threshold [default: 0]
  -h --help                     show help
  """ % ["version", version(), "env_fasta", env_fasta])

  let args = docopt(doc, version=version(), argv=argv)
  let mapq = parse_int($args["--mapq"])

  var fasta: cstring 
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam

  open(bam, $args["<BAM-or-CRAM>"], threads=threads, index=true, fai=fasta)
  if bam.idx == nil:
    stderr.write_line ("count-reads: requires bam/cram index")

  var regions = bed_to_table($args["<BED>"])
  internal_count(bam, uint8(mapq), eflag, regions)
  return 0
