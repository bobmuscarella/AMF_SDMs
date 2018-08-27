# Automated species distribution modeling for AM Fungi

Bob Muscarella < <bob.muscarella@gmail.com> >

*(updated August 27, 2018)*

The file `Batch_run_AMF_SDMs.R` was used to automatically build SDMs for a list of AM Fungi at the global scale.  For each taxon, it builds models using climate + resource variables, climate-only variables, and resource-only variables.  It builds models with two methods of selecting background records: *target* and *random*.  The *target* method includes all occurrence records for all non-focal species as bacground points (i.e., the *target group background*).  The *random* method uses 10,000 random points in the study area.

The file `Summary_Plots.R` can be used to generate some basic summary plots of which models were selected, and which variables were most important.
