Script to compare the results from 2pt_pipeline and WLpipe on DES Y1 3x2pt analysis. 

To run, edit the path for the WLpipe and 2pt_pipeline results:

```
wlpipe_outpath= <your_path>
desy1_outpath = <your_path>
```

Then, run:

```bash
python plot_output.py
```
and 

```bash
python compare_output.py
```
