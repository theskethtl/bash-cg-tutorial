# tuttlelab
Some useful TuttleLab scripts. WX.pdb, make_box.py, count_water.py, eq_parameter_file.mdp, min_parameter_file.mdp all written by Travis Hesketh, TuttleLab, 2019. LICENSE applies to these files ONLY. See respective authors for other files.

top_all36_prot.rtf part of the CHARMM parameter file distribution, (C) MacKerrell Lab, University of Maryland, 2018
http://mackerell.umaryland.edu/charmm_ff.shtml

martinize.py (http://cgmartini.nl/images/tools/martinize/martinize-2.6/martinize.py), martini_v2.2.itp (http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp), and ions.itp (http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp) all (C) University of Groningen, found on http://cgmartini.nl.

Modifications to martinize.py:
  ```sed -i "1837d; 1836d" martinize.py```

  Modifications to ions.itp (renamed from martini_v2.0_ions.itp)
    ```sed -i "s/NA+/NA/g, s/CL-/CL/g" martini_v2.0_ions.itp```
