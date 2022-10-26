# Pulsar
Pulsar scripts

- pazi
  - process_cal.py processes a folded pulsar calibration file with pazi.
  - zapped_chan.py outputs the pazi command to zap the channels that have been zapped in the input file.
  
- plotting
  - freq_plot.py generates pulse phase vs frequency plots with pav. It first generates .ps files, converts them to pdfs (requires ps2pdf), combines all pdfs into one file, then deletes all .ps files.
  - pol_plot.py generates pulse phase vs flux density & polarisation with pav. 
  
- single_pulse
  - Contains scripts for single pulse searches.
  - run_split_v1/2.py splits the band into multiple sub-bands.
  - run_pipeline.py/run_pipeline_prepsubband.py runs the single pulse searching. The _prepsubband variant automatically determines the DM steps.
  - cal_samp0.py generates the samples and widths required for grab_pulse.py
  - grab_pulse.py extracts pulses from search mode files.
  - de-disperse.py de-disperses the extracted pulses.
