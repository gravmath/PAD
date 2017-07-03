# -*- coding: utf-8 -*-
import numpy as np
import os
import re

###############################################################################
###############################################################################
#  File Scraping Utility Functions
###############################################################################
###############################################################################
def construct_file_selector(obs,instrument,mode,descriptor,ver):
    if instrument == 'mec':
        pattern_str  = '%s_%s_%s_l2_%s_\d{8}_v%s.cdf'  % (obs,instrument,mode,descriptor,ver)
    else:
        pattern_str  = '%s_%s_%s_l2_%s_\d{14}_v%s.cdf' % (obs,instrument,mode,descriptor,ver)
    print 'encoding a search for %s!' % pattern_str
    return re.compile(pattern_str)
    
def scrape_files(regex_pattern,targ_dir):
    
    filtered_files   = []
    for root, dirs, files in os.walk(targ_dir):
        for file in files:
            if regex_pattern.match(file):
                full_filename = (root+'/'+file).replace('\\','/')
                filtered_files.append(full_filename)
                
    return np.sort(np.array(filtered_files).flatten())

def file_epoch(filename):
    basename = os.path.basename(filename)
    return basename.split('_')[5]


