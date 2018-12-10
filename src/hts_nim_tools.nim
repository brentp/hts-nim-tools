import os
import hts
import sequtils
import strutils
import tables
import algorithm
import docopt
import kexpr

import ./version

include ./bam_filter
include ./count_reads
include ./vcf_check
include ./copy_number_sampler

var progs = {"bam-filter": bam_filter,
             "count-reads": count_reads,
             "vcf-check": vcf_check,
             "copy-number-sampler": copy_number_sampler
             }.toTable
var helps = {"bam-filter": "filter BAM/CRAM/SAM files with a simple expression language",
             "count-reads": "count BAM/CRAM reads in regions given in a BED file",
             "vcf-check": "check regions of a VCF against a background for missing chunks",
             "copy-number-sampler": "sample BAM regions with a probability given in a BED file (for creating CN truth-sets)"
             }.toTable

proc main() =

  var args = commandLineParams()
  if len(args) < 1 or not progs.contains(args[0]):
    var hkeys = toSeq(keys(helps))
    sort(hkeys, proc(a, b: string): int =
      if a < b: return -1
      else: return 1
      )
    echo format("\nhts-nim utility programs.\nversion: $#\n", version())

    for k in hkeys:
      echo format("	• $1: $2", k & repeat(" ", 20 - len(k)), helps[k])
    echo ""
  else:
    var p = args[0]; args.delete(0)
    quit(progs[p](args))

when isMainModule:
  main()
