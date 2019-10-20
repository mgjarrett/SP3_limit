#!/bin/bash

#gauss_cases="gauss128_chebyshev64 gauss64_chebyshev64 gauss8_chebyshev128 gauss8_chebyshev16 gauss8_chebyshev64
#gauss2_chebyshev64 gauss4_chebyshev4 gauss4_chebyshev64 gauss4_chebyshev8 gauss8_chebyshev8"
gauss_cases="gauss128_chebyshev64 gauss64_chebyshev64 gauss8_chebyshev128 gauss8_chebyshev16 gauss8_chebyshev64"
yamamoto_cases="yamamoto2_chebyshev16 yamamoto_4_chebyshev16 yamamoto6_chebyshev16"

for casename in $gauss_cases
do
  echo ${casename}
  python read_logfiles.py $casename
done

#for casename in $yamamoto_cases
#do
#  echo ${casename}
#  python read_logfiles.py $casename
#done

