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
def construct_file_selector(obs,instrument,mode,level,descriptor,timestamp):
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
          timestamp:  either 14 for 'yyyymmddHHMMSS' or 8 for 'yyyymmdd'

       Returns:
           compiled regex pattern used to match files
       
       Example use:  
           mms1_pattern = 
              construct_file_selector('mms1','fpi','fast','l1a','des-cnts',14)
              
           mms1_pattern = 
              construct_file_selector('mms1','fgm','brst','l2','',14)
           
       
       Note:  All available versions are encoded.  Other functions eliminate
              the older/earlier ones.
    """
    if descriptor == '':
        pattern_str  = '%s_%s_%s_%s_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,level,timestamp)
    else:
        pattern_str  = '%s_%s_%s_%s_%s_\d{%d}_v+\d.+\d.+\d.cdf'\
                       % (obs,instrument,mode,level,descriptor,timestamp)
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
    """Helper function designed to filer a list of files and eliminate older 
       versions of duplicate names. 
       
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
               
#############################################################################
def config_directories(basedir):
    """A helper function design to help set the appropriate search_dirs
       for the various instruments.
       
       Arguments:
          base_dir:     either '/fpiprd1/data/ftp/' if running on the server
                        or  'z:/data/ftp' is running on Windows via Samba

       Returns:
          dictionary with the following keys (<n> = 1,2,3,4 - for mms<n>):

              dsp_bpsd<n>:           search_dir for dsp for electric field
                                     for mms<n>
              dsp_epsd<n>:           search_dir for dsp for electric field
                                     for mms<n>
              edp_scpot<n>:          search_dir for fast s/c potential 
                                     for mms<n>
              fgm_brst<n>:           search_dir for fgm burst for mms<n>
              fgm_srvy<n>:           search_dir for fgm survey for mms<n>
              fpi_brst_des-dist<n>:  search_dir for fpi burst electron 
                                     distributions for mms<n>
              fpi_brst_dis-dist<n>:  search_dir for fpi burst ion distributions
                                     for mms<n>
              fpi_fast_des-dist<n>:  search_dir for fpi burst electron 
                                     distributions for mms<n>
              fpi_fast_dis-dist<n>:  search_dir for fpi burst ion distributions
                                     for mms<n>
              fpi_brst_des-moms<n>:  search_dir for fpi burst electron moments
                                     for mms<n>
              fpi_brst_dis-moms<n>:  search_dir for fpi burst ion moments
                                     for mms<n>
              fpi_fast_des-moms<n>:  search_dir for fpi burst electron momentss
                                     for mms<n>
              fpi_fast_dis-moms<n>:  search_dir for fpi burst ion moments
                                     for mms<n>
              mec_srvy_epht89d<n>:   search_dir for mec survey
     
       Example use:  
           my_search_dict = config_directories('z:/data/ftp/')
       
       Note:  Only l2 directories encoded
    """
    
    search_dirs = {}
    for obs in ['mms1','mms2','mms3','mms4']:
        obs_path       = basedir+obs+'/'
        obs_num        = obs[-1]
       
        search_dirs['dsp_fast_bpsd'+obs_num]       = obs_path+'dsp/fast/l2/bpsd/'
        search_dirs['dsp_fast_epsd'+obs_num]       = obs_path+'dsp/fast/l2/epsd/'
        search_dirs['edp_fast_scpot'+obs_num]      = obs_path+'edp_spdf/fast/l2/scpot/'
        search_dirs['fgm_brst'+obs_num]            = obs_path+'fgm/brst/l2/'
        search_dirs['fgm_srvy'+obs_num]            = obs_path+'fgm/brst/l2/srvy/l2/'    
        search_dirs['fpi_brst_des-dist'+obs_num]   = obs_path+'fpi/brst/l2/des-dist/'
        search_dirs['fpi_brst_dis-dist'+obs_num]   = obs_path+'fpi/brst/l2/dis-dist/'
        search_dirs['fpi_fast_des-dist'+obs_num]   = obs_path+'fpi/fast/l2/des-dist/'
        search_dirs['fpi_fast_dis-dist'+obs_num]   = obs_path+'fpi/fast/l2/dis-dist/'
        search_dirs['fpi_brst_des-moms'+obs_num]   = obs_path+'fpi/brst/l2/des-moms/'
        search_dirs['fpi_brst_dis-moms'+obs_num]   = obs_path+'fpi/brst/l2/dis-moms/'
        search_dirs['fpi_fast_des-moms'+obs_num]   = obs_path+'fpi/fast/l2/des-moms/'
        search_dirs['fpi_fast_dis-moms'+obs_num]   = obs_path+'fpi/fast/l2/dis-moms/'
        search_dirs['mec_srvy_epht89d'+obs_num]    = obs_path+'mec/srvy/l2/epht89d/'
 
    return search_dirs
               
               
###############################################################################
def get_my_files(obs,instrument,mode,level,descriptor,timestamp,search_dir):
    """The core function for creating a dictionary of file indices that
       meet certain parameters.
       
       Arguments:
          obs:         'mms1', 'mms2', 'mms3', 'mms4'
          instrument:  'fgm', 'fpi', (others include 'edp', etc. but they ,base
                        aren't on our server - caveat emptor)
          mode:        'brst','fast', or 'srvy'
          level:       'l1a', 'l1b','l1c','l2'
          descriptor:  'des-cnts',  'dis-cnts', 'des-moms', 'dis-moms', 
                       'des-debug', 'dis-debug'
          timestamp:   either 14 for 'yyyymmddHHMMSS' or 8 for 'yyyymmdd'

       Returns:
           a filtered file dictionary of latest versions
       
       Example use:  
           mms1_file = 
              get_my_files('mms1','fpi','fast','l1a','des-cnts',14)
              
           mms1_file = 
              get_my_files('mms1','fgm','brst','l2','',14)
           
       
       Note:  N/A
    """
    my_pattern = construct_file_selector(obs,instrument,mode,level,descriptor,timestamp)
    my_list    = scrape_files(my_pattern,search_dir)
    my_dict    = construct_file_dict(my_list)    

    return my_dict               

###############################################################################
def get_my_l2_files(obs,instrument,mode,descriptor,base_dir):
    """The core function for creating a dictionary of file indices that
       meet certain parameters.
       
       Arguments:
          obs:         'mms1', 'mms2', 'mms3', 'mms4'
          instrument:  'fgm', 'fpi', (others include 'edp', etc. but they ,base
                        aren't on our server - caveat emptor)
          mode:        'brst','fast', or 'srvy'
          descriptor:  'des-cnts',  'dis-cnts', 'des-moms', 'dis-moms', 
                       'des-debug', 'dis-debug'
          timestamp:  either 14 for 'yyyymmddHHMMSS' or 8 for 'yyyymmdd'

       Returns:
           a filtered file dictionary of latest versions
       
       Example use:  
           mms1_file = 
              get_my_l2_files('mms1','fpi','fast','des-cnts','/fpiprd1/data/ftp/')
              
           mms1_file = 
              get_my_l2_files('mms1','fgm','brst','','z:/data/ftp/')
           
       
       Note:  N/A
    """

    obs_num    = obs[-1]
    my_dirs    = config_directories(base_dir)

    if descriptor == '':
        my_key = instrument+'_'+mode+obs_num
    else:
        my_key = instrument+'_'+mode+'_'+descriptor+obs_num

    if instrument in ['dsp','bpsd','epsd','mec']:
        timestamp = 8
    if instrument in ['edp','edp_spdf','scpot','fpi']:
        timestamp = 14
    if instrument == 'fgm' and mode == 'brst':
        timestamp = 14
    if instrument == 'fgm' and mode == 'fast':
        timestamp = 8

    search_dir = my_dirs[my_key]

    my_pattern = construct_file_selector(obs,instrument,mode,'l2',descriptor,timestamp)
    my_list    = scrape_files(my_pattern,search_dir)
    print '\n'
    my_dict    = construct_file_dict(my_list)    
    print 'Found %s unique files' % (len(my_dict.keys()))

    return my_dict               

###############################################################################
def limit_time_range(start_time,stop_time,my_dict):
    """The core function for taking a file dictionary and limiting its time
       range.
       
       Arguments:
          start_time:  in 'YYYY-MM-DD hh:mm:ss' format ('2015-10-16 13:00:00')
          stop_time:   in 'YYYY-MM-DD hh:mm:ss' format ('2015-10-16 13:00:00')
          my_dict:     a output product from get_my_files, get_my_l2_files, or
                       construct_file_dict

       Returns:
           a filtered file list
       
       Example use:  
           my_fgm_files =  limit_time_range(start_time,stop_time,fpi_fast_des_moms_dict')
              
       Note:  N/A
    """
    filtered_list = []
    for time in sorted(my_dict.keys()):
        if time > start_time and time < stop_time:
            #print time, my_dict[time]['filename']
            filtered_list.append(my_dict[time]['filename'])

    return filtered_list