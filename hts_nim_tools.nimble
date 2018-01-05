# Package
version       = "0.1.5"
author        = "Brent Pedersen"
description   = "hts-nim command-line tools"
license       = "MIT"

# Dependencies
requires "nim >= 0.17.2", "c2nim >= 0.9.10", "docopt", "lapper", "hts", "kexpr"

srcDir = "src"
bin = @["hts_nim_tools"]

task named_build, "custom build":
  mkdir "bin"
  exec "nimble c --out:bin/hts-nim-tools src/hts_nim_tools"
