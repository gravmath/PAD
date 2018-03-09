# -*- coding: utf-8 -*-
import numpy as np
import os
import re

###############################################################################
###############################################################################
#  File Scraping Utility Functions
#
#
#  last modified - 3/9/2018
###############################################################################
###############################################################################
###############################################################################
def construct_file_selector(obs,instrument,mode,level,descriptor,timestring):
    """Helper function that constructs a MMS regex file pattern
       for looking on the server for files.
       
       Arguments:
          obs:         'mms1', 'mms2', 'mms3', 'mms4'
          instrument:  'fgm', 'fpi', (others include 'edp', etc. but they 
                        aren't on our server - caveat emptor)
          mode:        'brst','fast', or 'srvy'
          level:       'l1a', 'l1b','l1c','l2'
          descriptor:  'des-cnts',  'dis-cnts', 'des-moms', 'dis-moms', 
                       'des-debug', 'dis-debug'
          timestring:  either 14 for 'yyyymmddHHMMSS' or 8 for 'yyyymmdd'

       Returns:
           compiled regex pattern used to match files
       
       Example use:  
           mms1_pattern = 
              construct_file_selector('mms1','fpi','fast','l1a','des-cnts',14)
       
       Note:  All available versions are encoded.  Other functions eliminate
              the older/earlier ones.
    """
    if descriptor == '':
        pattern_str  = '%s_%s_%s_%s_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,level,timestring)
    else:
        pattern_str  = '%s_%s_%s_%s_%s_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,level,descriptor,timestring)
    print 'Encoding a search for %s!' % pattern_str
    return re.compile(pattern_str)

    
 
###############################################################################
def scrape_files(regex_pattern,targ_dir):
    """Helper function designed to take a regex_pattern encoded by 
       construct_file_selector and return a list of all files that match
       the pattern. 
       
       Arguments:
          regex_pattern:    any valid file pattern for filtering results
          targ_dir:         directory in which to start the os.walk
          
       Returns:
          file_list:      1-D list of all matching files
           
       Example use:  
          my_file_list = scrape_files(my_pattern,my_directory)

       Notes:  
          the list has all versions of a given file.
    """
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
def construct_file_dict(master_list):
    """Helper function designed to a list of files and eliminate older versions
       of duplicate names. 
       
       Arguments:
          master_list:    a list of filenames
          
       Returns:
          file_dict:      dictionary of files with the following properties:
                          level 1:  key   = timestamp (e.g. '2018-03-02 22:00:00')
                                    value = level 2 dict
                          level 2:  
                                    keys   = ['filename', 'obs', 'version']
                                    values = fully decorated filename
                                             observatory designator (e.g. 'mms1')
                                             version tuple (e.g. (3,3,0)
                                             
      Example use:  
          my_file_dict = construct_file_dict(my_file_list)

       Notes:  
          none
    """
    file_dict = {}
    for item in master_list:
        record, ver = get_info(item)
        index       = record[5]  #formatted_timestamp
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
def get_info(filename):
    """
    A helper function that takes a filename and finds the naming parameters
    of the file.  These naming parameters are:

       Arguments:
           filename:   typical 
           targ_dir:        root directory where the os.walk starts
           
       Returns:
           filtered_files:  list of all files that match the pattern 
                            (includes path)
                            
       Example use:
           my_files = scrape_files(my_pattern,'c:/Yuggoth/'
       
       Notes:  
           1) obs        - observatory chosen from 'mms1', 'mms2', 'mms3', 'mms4'
           2) instr      - instrument (e.g. 'fgm', 'fpi', 'edp', etc.)
           3) mode       - instrument operational mode from 'brst', 'fast', 'srvy'
           4) level      - 'l1a', 'l1b', 'l1c', 'l2' or whatever else
           5) descriptor - 'des-cnts', 'dis-moms', etc.  (optional field)
           6) formatted_timestamp - 8 or 14 digits representing YYYYmmdd or YYYYmmddHHMMSS
           7) ver        - file version
           
           The first six parameters are returned as a tuple (for possible use in a database).
           The version is returned as a single object.
    """
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
def format_timestamp(timestamp):
    """
    A simple helper function that takes a SOC timestamp (14 or 8 digits) and
    returns it in a human-readable format.  For example:
    
        timestamp = 20161005
    
    returns as
    
        formatted_timestamp = 2016-10-05
        
    and
    
        timestamp = 20161005123456
        
    returns as
    
        formatted_timestamp = 2016-10-05 12:34:56
        
        
       Arguments:
          timestamp:  string of the timestamp in SOC format
        
       Returns:
          formatted_timestamp:    string of the timestamp in human format
                         
       Example use:
          formatted_timestamp = format_timestamp(timestamp)  

       Notes:  
          none        
    """
    timelength = len(timestamp)
    if timelength == 8:
        return timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]
    if timelength == 14:
        return timestamp[0:4] +'-'+timestamp[4:6]  +'-'+timestamp[6:8]+' '+\
               timestamp[8:10]+':'+timestamp[10:12]+':'+timestamp[12:14]

