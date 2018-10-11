# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 08:41:33 2018

@author: Brendan
"""

"""
#######################################################
### set date/wavelength to use in yaml config file ###
######################################################
"""

import yaml
import sys

count = int(sys.argv[1])

fname = 'specFit_config.yaml'

with open(fname) as f:
    cfg = yaml.load(f)

date = cfg['dates'][count]
wave = cfg['wavelengths'][count]
print(date, wave, count, flush=True)
    
cfg['date'] = date
cfg['wavelength'] = wave

with open(fname, 'w') as f:
    yaml.dump(cfg, f, default_flow_style=False)