import numpy as np
import os

###############################################################################
###############################################################################
#  File Scraping
###############################################################################
###############################################################################
def scrape_files(regex_pattern,targ_dir):
    
    filtered_files   = []
    for root, dirs, files in os.walk(targ_dir):
        for file in files:
            if regex_pattern.match(file):
                full_filename = (root+'/'+file).replace('\\','/')
                filtered_files.append(full_filename)
                
    return np.array(filtered_files).flatten()
