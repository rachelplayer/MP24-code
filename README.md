A Central Limit Approach for Ring-LWE Noise Analysis
==========================================================

This repository contains the code used to generate Tables 1--4 in [MP24].

Experimental data in Tables 1 and 2 was obtained using the provided HElib code [HElib] (with HElib version 2.2.1). Experimental data in Tables 3 and 4 was obtained using the provided Microsoft SEAL code [SEAL] (with SEAL version 4.0). Our noise growth estimates in Tables 1--4 were obtained using a python script.


Running the code
----------------

**Heuristics** 
The python script `generate_bgv_heuristics_tables.py` generates the noise growth estimates reported in the columns `[CLP20]` and `Ours` in Tables 1--4. It is best run using SageMath [SAGE]. 

Within Sage:
`load("generate_bgv_heuristics_tables.py")`

**HElib**
The HElib files `BGV_clp20.cpp` (for Table 1) and `BGV_deep.cpp` (for Table 2) were developed to run with HElib (version 2.2.1). With that version of HElib installed, add the folders `BGV_CLP20` and `BGV_deep` to the folder HElib/examples/. These files can then be compiled and run as for the other HElib examples. 

In /HElib/examples:
`cmake .`
`make`

Then in /HElib/examples/bin:
`./BGV_deep`

Note that the files `BGV_clp20.cpp` and `BGV_deep.cpp` require a slight modification to the Ctxt class, namely that the `Ctxt::tensorProduct()` function is made public.

**SEAL**
The provided files `4_bgv_basics_CLP20.cpp` (for Table 3) and `4_bgv_basics_bgv_deep.cpp` (for Table 4) were developed to run with SEAL (version 4.0). With that version of SEAL installed, they can be swapped in for the file `4_bgv_basics.cpp` in the SEAL examples (SEAL/native/examples) and compiled and run as for the original SEAL examples.

In SEAL/
`cmake -S . -B build -DSEAL_BUILD_EXAMPLES=ON`
`cmake --build build`

In /SEAL/build/bin
`./sealexamples`


Bibliography
------------
[CLP20] Anamaria Costache, Kim Laine, Rachel Player. Evaluating the effective- ness of heuristic worst-case noise analysis in FHE. In ESORICS 2020. Preprint available at: https://eprint.iacr.org/2019/493

[HElib] HElib. https://github.com/homenc/HElib

[MP24] Sean Murphy, Rachel Player. IACR Communications in Cryptology. Volume 1, No. 2, 2024. Preprint available at: https://eprint.iacr.org/2019/452.pdf 

[SAGE] https://www.sagemath.org

[SEAL] Microsoft SEAL (release 4.0). Microsoft Research, Redmond, WA. March, 2022. https://github.com/microsoft/SEAL