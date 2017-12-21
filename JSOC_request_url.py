# -*- coding: utf-8 -*-
"""
Created on Fri Nov 03 14:35:37 2017

@author: Brendan
"""

from __future__ import absolute_import, division, print_function
#import os
import drms
import numpy as np

"""
Creates a JSOC export request.  Saves the resulting URL in file that is loaded
by MPI download script.
"""

import sys

directory = sys.argv[1]
date = sys.argv[2]
wavelength = sys.argv[3]
t_start = sys.argv[4]
duration = sys.argv[5]

# email registered at JSOC (or set as variable) 
# If you have not registered your email yet, you can do this on the
# JSOC website at: http://jsoc.stanford.edu/ajax/register_email.html
email = 'bgallag6@gmu.edu'

# Use 'as-is' instead of 'fits', if record keywords are not needed in the
# FITS header. This greatly reduces the server load!
export_protocol = 'fits'
#export_protocol = 'as-is'

#out_dir = 'S:/FITS/20170505/1600/' 

# Create DRMS client, use debug=True to see the query URLs.
c = drms.Client(verbose=True)

q_date = '%s.%s.%s' % (date[0:4], date[4:6], date[6:8])

# Data export query string
#qstr = 'aia.lev1_uv_24s[2011.09.02_00:00/1h][1600]'
#qstr = 'aia.lev1_euv_12s[%s_%s/%sh][%s]' % (q_date, t_start, duration, wavelength)
if wavelength == 'UV':
    qstr = 'aia.lev1_uv_24s[%s_%s/%sh]' % (q_date, t_start, duration)
if wavelength == '1600' or wavelength == '1700':
    qstr = 'aia.lev1_uv_24s[%s_%s/%sh][%s]' % (q_date, t_start, duration, wavelength) 

print('Data export query:\n  %s\n' % qstr)

# Submit export request using the 'fits' protocol
print('Submitting export request...')
r = c.export(qstr, method='url', protocol=export_protocol, email=email)

r.wait(timeout=2500, sleep=60)
# Print request URL.
print('\nRequest URL: %s' % r.request_url)
print('%d file(s) available for download.\n' % len(r.urls))

r_url = str(r.request_url)

np.save('%s/FITS/%s/%s_request_url.npy' % (directory, date, date), r_url)