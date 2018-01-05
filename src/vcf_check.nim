import os
import hts
import docopt
import strutils

import ./version

proc which_bin(pos:int, chunk_size:int): int {.inline.} =
  return int(pos / chunk_size)

proc extend(L:var seq[uint32], newIdx:int) {.inline.} =
  var newLen = newIdx + 1
  if newLen <= L.len: return
  var o = L.len
  L.set_len(newLen)
  for i in o..<newLen:
    L[i] = 0

proc count(v:VCF, chrom:string, maf:float32, chunk_size:int): seq[uint32] =
  ## count the query for this chrom
  var counts = new_seq_of_cap[uint32](1000)
  var afs = new_seq[float32](3)
  #stderr.write_line("[vcf-check] querying chrom:" & chrom)
  var n = 0
  var nm = 0
  for variant in v.query(chrom):
    n += 1
    if variant.info.floats("AF", afs) != Status.OK: continue
    if max(afs) < maf: continue
    nm += 1
    var bin = which_bin(variant.start, chunk_size)
    counts.extend(bin)
    inc(counts[bin])
  #stderr.write_line("n:" & $n & " nm:" & $nm)
  return counts

proc write_counts(chrom:string, q_counts: var seq[uint32], db_counts: var seq[uint32], chunk_size:int) =
  if len(q_counts) < len(db_counts):
    q_counts.extend(db_counts.len - 1)
  if len(db_counts) < len(q_counts):
    db_counts.extend(q_counts.len - 1)

  for i in 0..<q_counts.len:
    if q_counts[i] == 0 and db_counts[i] == 0: continue
    var start = i * chunk_size
    #var stop = start + chunk_size
    echo chrom & "\t" & $start & "\t" & $q_counts[i] & "\t" & $db_counts[i]

proc check_vcf(q:string, db:string, chunk_size:int, maf:float32) =
  var last_chrom:string
  var afs = new_seq[float32](3)
  var dbv:VCF
  var q_counts = new_seq_of_cap[uint32](1000)
  discard open(dbv, db)
  var qv:VCF
  discard open(qv, q)

  for variant in qv:
    if last_chrom == nil:
      last_chrom = $variant.CHROM
    if last_chrom != $variant.CHROM:
      stderr.write_line("[vcf-check] chrom:" & last_chrom)
      var db_counts = count(dbv, last_chrom, maf, chunk_size)
      write_counts(last_chrom, q_counts, db_counts, chunk_size)
      q_counts = new_seq_of_cap[uint32](1000)

      last_chrom = $variant.CHROM

    if variant.info.floats("AF", afs) != Status.OK: continue
    if max(afs) < maf: continue

    var bin = which_bin(variant.start, chunk_size)
    q_counts.extend(bin)
    inc(q_counts[bin])

  if last_chrom != nil:
    stderr.write_line("[vcf-check] chrom:" & last_chrom)
    var db_counts = count(dbv, last_chrom, maf, chunk_size)
    write_counts(last_chrom, q_counts, db_counts, chunk_size)

proc vcf_check(argv: var seq[string]): int =
  let doc = format("""

Check a VCF against a background to make sure that there
are no large missing chunks.

  $version

  Usage: vcf-check [options] <BACKGROUND_VCF> <VCF>


Arguments:                                                                                                                                                 
  <BACKGROUND_VCF>        population VCF/BCF with expected sites
  <VCF>                   query VCF/BCF to check

Options:

  -c --chunk <INT>        chunk size for genome [default: 100000]
  -m --maf <FLOAT>        allele frequency  cutoff [default: 0.1]
  """ % ["version", version()])

  let args = docopt(doc, version=version(), argv=argv)
  let chunk = parse_int($args["--chunk"])
  let maf = parse_float($args["--maf"])
  check_vcf($args["<VCF>"], $args["<BACKGROUND_VCF>"], chunk, maf.float32)
  return 0
