# -*- coding: utf-8 -*-
import numpy as np
import os
import re

###############################################################################
###############################################################################
#  File Scraping Utility Functions
###############################################################################
###############################################################################
def construct_file_selector(cursor,obs,instrument,mode,descriptor):
    command = 'SELECT timestamp from Instrument_Description\
                    WHERE instrument = "%s" and mode = "%s" and descriptor="%s";' \
               % (instrument,mode,descriptor)
    cursor.execute(command)
    timestamp = cursor.fetchall()[0][0]
    if descriptor == '':
        pattern_str  = '%s_%s_%s_l2_\d{%d}_v+\d.+\d.+\d.cdf'    % (obs,instrument,mode,timestamp)
    else:
        pattern_str  = '%s_%s_%s_l2_%s_\d{%d}_v+\d.+\d.+\d.cdf' % (obs,instrument,mode,descriptor,timestamp)
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


def get_info(filename):
    basename   = os.path.basename(filename)
    components = basename.split('_')
    
    obs            = components[0]
    instr          = components[1]
    mode           = components[2]
    level          = components[3]
    if len(components) == 7: #has a descriptor
        descriptor = components[4]
        timestamp  = components[5]
        ver_string = components[6]
    if len(components) == 6: #no descriptor
        descriptor = ''
        timestamp  = components[4]
        ver_string = components[5]

    ver = tuple([int(v) for v in ver_string[1:-4].split('.')])

    return (obs, instr, mode, level, descriptor, timestamp), ver
    
def construct_file_dict(master_list):
    file_dict = {}
    for item in master_list:
        index, ver = get_info(item)
        if index in file_dict:
            curr_ver = file_dict[index]['version']
            if ver > curr_ver:
                file_dict[index]['version']  = ver
                file_dict[index]['filename'] = item
        else:
            file_dict[index] = {}
            file_dict[index]['version']  = ver
            file_dict[index]['filename'] = item
    return file_dict

def populate_db(cursor,table_name,file_dict):
    cursor.execute('DROP TABLE %s;',table_name)
    for k in sorted(file_dict.keys()):
        command = 'INSERT INTO   %s (obs, instrument, mode, level, descriptor, timestamp) \
                          VALUES    (%s, %s, %s, %s, %s, %s, %s);' % k
        cursor.execute(command)
    
    cursor.commit()

    
