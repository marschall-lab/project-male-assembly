## Build instructions

Build the following containers sequentially. The final container is a working
`Verkko` installation with all necessary dependencies.

1. base_env.def
2. jemalloc_mummer.def
3. ga_mbg.def
4. verkko.def
    - current version in def file: `#1a00b60`

Note that the auto discovery of `Verkko` to detect its setup folder and 
associated tools is a bit inscrutable. Hence, it's safest to run `Verkko`
as `/repos/verkko/bin/verkko`, i.e.

`./verkko.sif /repos/verkko/bin/verkko [parameters]`
