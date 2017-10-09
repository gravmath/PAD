# -*- coding: utf-8 -*-
import numpy as np
import os
import re

from DB import *   #get the database connection

###############################################################################
###############################################################################
#  File Scraping Utility Functions
###############################################################################
###############################################################################


###############################################################################
def get_instrument_timestring(instrument,mode,descriptor):
    my_command = 'select timestring from instruments where instrument = "%s"\
                  and mode = "%s" and descriptor = "%s";'\
                 % (instrument,mode,descriptor)
    instrument_cursor.execute(my_command)
    return int(instrument_cursor.fetchall()[0][0])


###############################################################################
def construct_file_selector(obs,instrument,mode,descriptor):
    timestring = get_instrument_timestring(instrument,mode,descriptor)
    if descriptor == '':
        pattern_str  = '%s_%s_%s_l2_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,timestring)
    else:
        pattern_str  = '%s_%s_%s_l2_%s_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,descriptor,timestring)
    print 'Encoding a search for %s!' % pattern_str
    return re.compile(pattern_str)

   
###############################################################################
def scrape_files(regex_pattern,targ_dir):
    filtered_files = []
    counter        = 0
    for root, dirs, files in os.walk(targ_dir):
        for file in files:
            if regex_pattern.match(file):
                full_filename = (root+'/'+file).replace('\\','/')
                filtered_files.append(full_filename)
            counter += 1
            if counter % 1000 == 0:
               print '*',
                
    return np.sort(np.array(filtered_files).flatten())


###############################################################################
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

    ver                 = tuple([int(v) for v in ver_string[1:-4].split('.')])
    formatted_timestamp = format_timestamp(timestamp)

    return (obs, instr, mode, level, descriptor, formatted_timestamp), ver
    
###############################################################################
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


###############################################################################
def format_timestamp(timestamp):
    timelength = len(timestamp)
    if timelength == 8:
        return timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]
    if timelength == 14:
        return timestamp[0:4] +'-'+timestamp[4:6]  +'-'+timestamp[6:8]+' '+\
               timestamp[8:10]+':'+timestamp[10:12]+':'+timestamp[12:14]

###############################################################################
def populate_db(table_name,file_dict):
    for k in sorted(file_dict.keys()):
        key_tuple  = k
        data_tuple = (file_dict[k]['filename'],file_dict[k]['version'])
        tot_tuple  = (table_name,)+key_tuple + data_tuple
        command    = 'INSERT INTO  %s VALUES\
                      ("%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s");'\
                      % tot_tuple
        instrument_cursor.execute(command)
    instrument_conn.commit()

  
###############################################################################
def clean_db(table_name):
    #note that sqlite doesn't have an explicit TRUNCATE command
    #rather one uses 'DELETE FROM table_name;'
    #see https://www.techonthenet.com/sqlite/truncate.php
    instrument_cursor.execute('DELETE FROM %s;' % table_name)
    instrument_conn.commit()

###############################################################################
def scrape_all_data_CDFs():
   #specify the table
   CDF_table = 'CDF_files'

   #setup the instrument list
   instrument_cursor.execute('Select * from instruments;')
   instrument_conn.commit()
   mon = instrument_cursor.fetchall()

   #setup the observatory list
   obs = ['mms1','mms2','mms3','mms4']

   #clean out the database
   clean_db(CDF_table)

   #scrape for everything
   for m in mon:
       for o in obs:
           instrument  = m[0]
           mode        = m[1]
           descriptor  = m[2]
           file_select = construct_file_selector(o,instrument,mode,descriptor)
           start_dir   = PAD_data_dir+o+'/'+instrument+'/'+mode+'/'+'l2'
           print 'Starting scrape in %s' % start_dir
           raw_list    = scrape_files(file_select,start_dir)
           print
           file_dict   = construct_file_dict(raw_list)
           print  'Found %d %s-%s-%s-%s files' %\
                  (len(file_dict),o,instrument,mode,descriptor)
           populate_db(CDF_table,file_dict)
           print
