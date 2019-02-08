# -*- coding: utf-8 -*-
import cPickle  as pickle
import datetime as dt
import numpy    as np
import os
import re
import requests
import shutil

schiff_auth        = ('cschiff', 'es_fr_GSFC8dLASP')
pre_string         = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/file_names/science?'
start_pre          = 'start_date='
end_pre            = 'end_date='
descriptor_pre     = 'descriptor='
sc_id_pre          = 'sc_id='
data_level         = 'data_level=l2'
data_rate_mode_pre = 'data_rate_mode='
cnkt               = '&'
instrument_pre     = 'instrument_id='
download_pre       = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/download/science?file='

survey_sets = [['dsp','fast','bpsd'],
               ['dsp','fast','epsd'],
               ['edp','fast','scpot'],
               ['fgm','srvy',''],
               ['mec','srvy','epht89d']]
               
fast_sets   = [['fpi','fast','des-moms'],
               ['fpi','fast','dis-moms'],
               ['hpca','srvy','moments']]

fast_dist   = [['fpi','fast','des-dist'],
               ['fpi','fast','dis-dist'],
               ['hpca','srvy','ion']]

brst_sets   = [['edp','brst','dce'],
               ['fgm','brst',''],
               ['fpi','brst','des-moms'],
               ['fpi','brst','dis-moms'],
               ['hpca','brst','moments']]

brst_dist   = [['fpi','brst','dis-dist'],
               ['fpi','brst','des-moms'],
               ['hpca','brst','ion']]

mini_sets = [['dsp','fast','bpsd'],
             ['dsp','fast','epsd']]               
###############################################################################
###############################################################################
#  File Scraping Utility Functions
#
#
#  last modified - 11/21/2018
#
#  The function listings fall into n broad categories
#  1)  find and report back on files on the disk (scrape files to a list of names)
#      - construct_file_selector
#      - scrape_files
#      - construct_mms_file_dict
#      - get_mms_filename_info
#      - format_mms_timestamp
#      - config_mms_directories
#      - record_my_mms_files
#      - record_my_mms_l2_files
#  2)  selecting a subset of the files from the application of 1
#      - limit_mms_time_range
#      - create_mms_inventory_dict
#      - scrape_a_drive
#      - filter_for_prj
#      - load_filelist_for_prj
#  3)  download files from the MMS SDC
#      - download_file
#      - SDC_data_summary
#      - retrieve_SDC_data
###############################################################################
###############################################################################
###############################################################################
def construct_file_selector(obs,instrument,mode,level,descriptor,timestamp,fh="None"):
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
    if fh == "None":
        print "Encoding a search for %s!" % pattern_str
    else:
        fh.write("Encoding a search for %s!\n" % pattern_str)
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
def construct_mms_file_dict(master_list):
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
          my_file_dict = construct_mms_file_dict(my_file_list)

       Notes:  
          none
    """
    file_dict = {}
    for item in master_list:
        record, ver = get_mms_filename_info(item)
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
def get_mms_filename_info(filename):
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
    formatted_timestamp = format_mms_timestamp(timestamp)

    return (obs, instr, mode, level, descriptor, formatted_timestamp), ver
    
###############################################################################
def format_mms_timestamp(timestamp):
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
          formatted_timestamp = format_mms_timestamp(timestamp)  

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
def config_mms_directories(basedir,file_source):
    """A helper function design to help set the appropriate search_dirs
       for the various instruments.
       
       Arguments:
          base_dir:     either '/fpiprd1/data/ftp/' if running on the server
                        or  'z:/data/ftp' is running on Windows via Samba

       Returns:
          dictionary with the following keys (<n> = 1,2,3,4 - for mms<n>):

              dsp_fast_bpsd<n>:      search_dir for dsp for electric field
                                     for mms<n>
              dsp_fast_epsd<n>:      search_dir for dsp for electric field
                                     for mms<n>
              edp_fast_scpot<n>:     search_dir for fast s/c potential 
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
           my_search_dict = config_mms_directories('z:/data/ftp/','server')
       
       Note:  Only l2 directories encoded
    """
    
    search_dirs = {}
    for obs in ['mms1','mms2','mms3','mms4']:
        obs_path       = basedir+obs+'/'
        obs_num        = obs[-1]
       
        search_dirs['dsp_fast_bpsd'+obs_num]       = obs_path+'dsp/fast/l2/bpsd/'
        search_dirs['dsp_fast_epsd'+obs_num]       = obs_path+'dsp/fast/l2/epsd/'
        search_dirs['fgm_brst'+obs_num]            = obs_path+'fgm/brst/l2/'
        search_dirs['fgm_srvy'+obs_num]            = obs_path+'fgm/srvy/l2/'    
        search_dirs['fpi_brst_des-dist'+obs_num]   = obs_path+'fpi/brst/l2/des-dist/'
        search_dirs['fpi_brst_dis-dist'+obs_num]   = obs_path+'fpi/brst/l2/dis-dist/'
        search_dirs['fpi_fast_des-dist'+obs_num]   = obs_path+'fpi/fast/l2/des-dist/'
        search_dirs['fpi_fast_dis-dist'+obs_num]   = obs_path+'fpi/fast/l2/dis-dist/'
        search_dirs['fpi_brst_des-moms'+obs_num]   = obs_path+'fpi/brst/l2/des-moms/'
        search_dirs['fpi_brst_dis-moms'+obs_num]   = obs_path+'fpi/brst/l2/dis-moms/'
        search_dirs['fpi_fast_des-moms'+obs_num]   = obs_path+'fpi/fast/l2/des-moms/'
        search_dirs['fpi_fast_dis-moms'+obs_num]   = obs_path+'fpi/fast/l2/dis-moms/'
        search_dirs['mec_srvy_epht89d'+obs_num]    = obs_path+'mec/srvy/l2/epht89d/'
        if file_source == 'server':
            search_dirs['edp_fast_scpot'+obs_num]      = obs_path+'edp_spdf/fast/l2/scpot/'
            search_dirs['edp_brst_dce'+obs_num]        = obs_path+'edp_spdf/brst/l2/dce/'
            search_dirs['hpca_brst_ion'+obs_num]       = obs_path+'hpca_spdf/brst/l2/ion/'        
            search_dirs['hpca_brst_moments'+obs_num]   = obs_path+'hpca_spdf/brst/l2/moments/'        
            search_dirs['hpca_srvy_ion'+obs_num]       = obs_path+'hpca_spdf/srvy/l2/ion/'        
            search_dirs['hpca_srvy_moments'+obs_num]   = obs_path+'hpca_spdf/srvy/l2/moments/'        
        if file_source == 'local':
            search_dirs['edp_fast_scpot'+obs_num]      = obs_path+'edp/fast/l2/scpot/'
            search_dirs['edp_brst_dce'+obs_num]        = obs_path+'edp/brst/l2/dce/'
            search_dirs['hpca_brst_ion'+obs_num]       = obs_path+'hpca/brst/l2/ion/'        
            search_dirs['hpca_brst_moments'+obs_num]   = obs_path+'hpca/brst/l2/moments/'        
            search_dirs['hpca_srvy_ion'+obs_num]       = obs_path+'hpca/srvy/l2/ion/'        
            search_dirs['hpca_srvy_moments'+obs_num]   = obs_path+'hpca/srvy/l2/moments/'        
        
 
    return search_dirs
               
###############################################################################
def record_my_mms_files(obs,instrument,mode,level,descriptor,timestamp,search_dir,fh="None"):
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
              record_my_mms_files('mms1','fpi','fast','l1a','des-cnts',14)
              
           mms1_file = 
              record_my_mms_files('mms1','fgm','brst','l2','',14)
           
       
       Note:  N/A
    """
    my_pattern = construct_file_selector(obs,instrument,mode,level,descriptor,timestamp,fh)
    my_list    = scrape_files(my_pattern,search_dir)
    my_dict    = construct_mms_file_dict(my_list)    
    if fh == "None":
        print "Found %s unique files in %s\n" % (len(my_dict.keys()),search_dir)
    else:
       fh.write("Found %s unique files in %s\n\n" % (len(my_dict.keys()),search_dir))

    return my_dict               

###############################################################################
def record_my_mms_l2_files(obs,instrument,mode,descriptor,base_dir,file_source,fh="None"):
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
              record_my_mms_l2_files('mms1','fpi','fast','des-cnts','/fpiprd1/data/ftp/')
              
           mms1_file = 
              record_my_mms_l2_files('mms1','fgm','brst','','z:/data/ftp/')
           
       
       Note:  N/A
    """

    obs_num    = obs[-1]
    my_dirs    = config_mms_directories(base_dir,file_source)

    if descriptor == '':
        my_key = instrument+'_'+mode+obs_num
    else:
        my_key = instrument+'_'+mode+'_'+descriptor+obs_num

    if instrument in ['dsp','bpsd','epsd','mec']:
        timestamp = 8
    elif instrument in ['edp','edp_spdf','hpca','hpca_spdf','scpot','fpi']:
        timestamp = 14
    elif instrument == 'fgm' and mode == 'brst':
        timestamp = 14
    elif instrument == 'fgm' and mode == 'srvy':
        timestamp = 8
    else:
        timestamp = 8


    search_dir = my_dirs[my_key]

    my_pattern = construct_file_selector(obs,instrument,mode,'l2',descriptor,timestamp,fh)
    my_list    = scrape_files(my_pattern,search_dir)
    my_dict    = construct_mms_file_dict(my_list)   
    if fh == "None":
        print "Found %s unique files in %s\n" % (len(my_dict.keys()),search_dir)
    else:
       fh.write("Found %s unique files in %s\n\n" % (len(my_dict.keys()),search_dir))
    
    return my_dict               

###############################################################################
def limit_mms_time_range(start_time,stop_time,my_dict):
    """The core function for taking a file dictionary and limiting its time
       range.
       
       Arguments:
          start_time:  in 'YYYY-MM-DD hh:mm:ss' format ('2015-10-16 13:00:00')
          stop_time:   in 'YYYY-MM-DD hh:mm:ss' format ('2015-10-16 13:00:00')
          my_dict:     a output product from record_my_mms_files, record_my_mms_l2_files, or
                       construct_mms_file_dict

       Returns:
           a filtered file list
       
       Example use:  
           my_fgm_files =  limit_mms_time_range(start_time,stop_time,fpi_fast_des_moms_dict')
              
       Note:  N/A
    """
    filtered_list = []
    for time in sorted(my_dict.keys()):
        if len(time) == 10:  #8 digit timestamp
            compare_time = time+' 00:00:00'
        else:
            compare_time = time
        if compare_time >= start_time and compare_time <= stop_time:
            #print time, my_dict[time]['filename']
            filtered_list.append(my_dict[time]['filename'])

    return filtered_list
    
###############################################################################
def create_mms_inventory_dict(obs_name):
    """
    Helper function used to create datastructures for
    housing a list of relevant files for use by the Munger
    """
    my_file_dict = {'bpsd_f':[],
                    'epsd_f':[],
                    'dce_b':[],
                    'scpot_f':[],                    
                    'fgm_b':[],
                    'fgm_s':[],
                    'edist_b':[],
                    'emoms_b':[],
                    'idist_b':[],
                    'imoms_b':[],
                    'edist_f':[],
                    'emoms_f':[],
                    'idist_f':[],
                    'imoms_f':[],
                    'hpca_ion_b':[],
                    'hpca_moms_b':[],
                    'hpca_ion_s':[],
                    'hpca_moms_s':[],
                    'mec_s':[],
                    'ecnts_b':[]
                    }
    return my_file_dict
    
###############################################################################
def scrape_a_drive(base_dir,master_filelist_filename,file_source):
    """
    Helper function that scans the base_dir to find all instances
    used to create datastructures for
    housing a list of relevant files for use by the Munger
    
    obs                      - string with either 'mms1', 'mms2', etc.
    base_dir                 - starting (hence base) directory where the CDFs are found
    master_filelist_filename - filename of the pickle file where all the filenames are kept
    file_source              - string with either 'server' or 'local' (should be 'local' unless automating); the difference is
                               due to the soft links on the server and, as a result, the deviation from the SDC conventions

    """    
    workflow_log = open(base_dir+"workflow_log.dat",'a')
    now          = dt.datetime.now()
    workflow_log.write("\n\n******************************************************************************\n")     
    workflow_log.write("Ready to scrape a drive at %s\n" % (now,))
    workflow_log.write("Scraping together files from %s and saving to %s\n\n" % (base_dir,master_filelist_filename))
    
    try:
        scanned_files = open(master_filelist_filename, "w")
    except:
        workflow_log.write("Can't open the scraped file list!\n")

    for obs in ['mms1','mms2','mms3','mms4']:
        bpsd_fast_dict         = record_my_mms_l2_files(obs,'dsp', 'fast','bpsd',    base_dir,file_source,workflow_log)
        epsd_fast_dict         = record_my_mms_l2_files(obs,'dsp', 'fast','epsd',    base_dir,file_source,workflow_log)
        dce_brst_dict          = record_my_mms_l2_files(obs,'edp', 'brst','dce',     base_dir,file_source,workflow_log)
        scpot_fast_dict        = record_my_mms_l2_files(obs,'edp', 'fast','scpot',   base_dir,file_source,workflow_log)
        fgm_brst_dict          = record_my_mms_l2_files(obs,'fgm', 'brst','',        base_dir,file_source,workflow_log)
        fgm_srvy_dict          = record_my_mms_l2_files(obs,'fgm', 'srvy','',        base_dir,file_source,workflow_log)        
        fpi_brst_des_dist_dict = record_my_mms_l2_files(obs,'fpi', 'brst','des-dist',base_dir,file_source,workflow_log)
        fpi_brst_des_moms_dict = record_my_mms_l2_files(obs,'fpi', 'brst','des-moms',base_dir,file_source,workflow_log)
        fpi_brst_dis_dist_dict = record_my_mms_l2_files(obs,'fpi', 'brst','dis-dist',base_dir,file_source,workflow_log)
        fpi_brst_dis_moms_dict = record_my_mms_l2_files(obs,'fpi', 'brst','dis-moms',base_dir,file_source,workflow_log)
        fpi_fast_des_dist_dict = record_my_mms_l2_files(obs,'fpi', 'fast','des-dist',base_dir,file_source,workflow_log)
        fpi_fast_des_moms_dict = record_my_mms_l2_files(obs,'fpi', 'fast','des-moms',base_dir,file_source,workflow_log)
        fpi_fast_dis_dist_dict = record_my_mms_l2_files(obs,'fpi', 'fast','dis-dist',base_dir,file_source,workflow_log)
        fpi_fast_dis_moms_dict = record_my_mms_l2_files(obs,'fpi', 'fast','dis-moms',base_dir,file_source,workflow_log)
        hpca_brst_ion_dict     = record_my_mms_l2_files(obs,'hpca','brst','ion',     base_dir,file_source,workflow_log)
        hpca_brst_moments_dict = record_my_mms_l2_files(obs,'hpca','brst','moments', base_dir,file_source,workflow_log)
        hpca_srvy_ion_dict     = record_my_mms_l2_files(obs,'hpca','srvy','ion',     base_dir,file_source,workflow_log)
        hpca_srvy_moments_dict = record_my_mms_l2_files(obs,'hpca','srvy','moments', base_dir,file_source,workflow_log)
        mec_srvy_epht89d_dict  = record_my_mms_l2_files(obs,'mec', 'srvy','epht89d', base_dir,file_source,workflow_log)
        fpi_brst_des_cnts_dict = record_my_mms_files(obs,'fpi','brst','l1a','des-cnts',14,base_dir,workflow_log)
        try:
            pickle.dump(bpsd_fast_dict,scanned_files)         
            pickle.dump(epsd_fast_dict,scanned_files)         
            pickle.dump(dce_brst_dict,scanned_files)          
            pickle.dump(scpot_fast_dict,scanned_files)    
            pickle.dump(fgm_brst_dict,scanned_files)
            pickle.dump(fgm_srvy_dict,scanned_files)        
            pickle.dump(fpi_brst_des_dist_dict,scanned_files) 
            pickle.dump(fpi_brst_des_moms_dict,scanned_files) 
            pickle.dump(fpi_brst_dis_dist_dict,scanned_files) 
            pickle.dump(fpi_brst_dis_moms_dict,scanned_files) 
            pickle.dump(fpi_fast_des_dist_dict,scanned_files) 
            pickle.dump(fpi_fast_des_moms_dict,scanned_files) 
            pickle.dump(fpi_fast_dis_dist_dict,scanned_files) 
            pickle.dump(fpi_fast_dis_moms_dict,scanned_files) 
            pickle.dump(hpca_brst_ion_dict,scanned_files) 
            pickle.dump(hpca_brst_moments_dict,scanned_files)         
            pickle.dump(hpca_srvy_ion_dict,scanned_files) 
            pickle.dump(hpca_srvy_moments_dict,scanned_files)         
            pickle.dump(mec_srvy_epht89d_dict,scanned_files)
            pickle.dump(fpi_brst_des_cnts_dict,scanned_files)
        except:
            workflow_log.write("Can't save the scraped file list for obs %s!\n" % (obs,))
            return False
        
    scanned_files.close()
    now = dt.datetime.now()
    workflow_log.write("Done with scraping at %s\n" % (now,))
    workflow_log.write("******************************************************************************\n\n")     
    workflow_log.close()
    return True
    
###############################################################################
def filter_for_prj(base_dir,time_filters,project_filelist_filename,master_filelist_filename):
    """
    File to restore all the data structures for projects
    
    time_filters              - dictionary of time ranges
    project_filelist_filename - name of the pickle file with a list of filenames germane to the project
    master_filelist_filename  - the filelist pickle file (either master or project)
    """
    workflow_log = open(base_dir+"workflow_log.dat",'a')
    now          = dt.datetime.now()
    workflow_log.write("\n\n******************************************************************************\n")     
    workflow_log.write("Preparing to load and filter for a project at %s\n" % (now,))    
    workflow_log.write("Reading in scraped files on %s\n" % (master_filelist_filename,))

    start_time_brst = time_filters['start_time_brst']
    stop_time_brst  = time_filters['stop_time_brst']
    start_time_fast = time_filters['start_time_fast']
    stop_time_fast  = time_filters['stop_time_fast']
    start_time_srvy = time_filters['start_time_srvy']
    stop_time_srvy  = time_filters['stop_time_srvy']
    
    workflow_log.write("")
    workflow_log.write("Survey files limited to:      %s and %s\n" % (start_time_srvy,stop_time_srvy))    
    workflow_log.write("Fast Survey files limited to: %s and %s\n" % (start_time_fast,stop_time_fast))        
    workflow_log.write("Burst files limited to:       %s and %s\n" % (start_time_fast,stop_time_fast))        

    try:
        scanned_files          = open(master_filelist_filename,"r")
    except:
        workflow_log.write("Can't load the scraped file list!\n")
        return False

    prj_files = {}
    
    for obs in ['mms1','mms2','mms3','mms4']:
        bpsd_fast_dict         = pickle.load(scanned_files)
        epsd_fast_dict         = pickle.load(scanned_files)
        dce_brst_dict          = pickle.load(scanned_files)
        scpot_fast_dict        = pickle.load(scanned_files)
        fgm_brst_dict          = pickle.load(scanned_files)
        fgm_srvy_dict          = pickle.load(scanned_files)        
        fpi_brst_des_dist_dict = pickle.load(scanned_files)
        fpi_brst_des_moms_dict = pickle.load(scanned_files)
        fpi_brst_dis_dist_dict = pickle.load(scanned_files)
        fpi_brst_dis_moms_dict = pickle.load(scanned_files)
        fpi_fast_des_dist_dict = pickle.load(scanned_files)
        fpi_fast_des_moms_dict = pickle.load(scanned_files)
        fpi_fast_dis_dist_dict = pickle.load(scanned_files)
        fpi_fast_dis_moms_dict = pickle.load(scanned_files)
        hpca_brst_ion_dict     = pickle.load(scanned_files)
        hpca_brst_moments_dict = pickle.load(scanned_files)
        hpca_srvy_ion_dict     = pickle.load(scanned_files)
        hpca_srvy_moments_dict = pickle.load(scanned_files)
        mec_srvy_epht89d_dict  = pickle.load(scanned_files)
        fpi_brst_des_cnts_dict = pickle.load(scanned_files)


        prj_files[obs] = create_mms_inventory_dict(obs)  
        
        prj_files[obs]['bpsd_f']      = limit_mms_time_range(start_time_srvy,stop_time_srvy,bpsd_fast_dict)
        prj_files[obs]['epsd_f']      = limit_mms_time_range(start_time_srvy,stop_time_srvy,epsd_fast_dict)
        prj_files[obs]['dce_b']       = limit_mms_time_range(start_time_brst,stop_time_brst,dce_brst_dict)
        prj_files[obs]['scpot_f']     = limit_mms_time_range(start_time_srvy,stop_time_srvy,scpot_fast_dict)
        prj_files[obs]['fgm_b']       = limit_mms_time_range(start_time_brst,stop_time_brst,fgm_brst_dict)
        prj_files[obs]['fgm_s']       = limit_mms_time_range(start_time_srvy,stop_time_srvy,fgm_srvy_dict)    
        prj_files[obs]['edist_b']     = limit_mms_time_range(start_time_brst,stop_time_brst,fpi_brst_des_dist_dict)
        prj_files[obs]['emoms_b']     = limit_mms_time_range(start_time_brst,stop_time_brst,fpi_brst_des_moms_dict)
        prj_files[obs]['idist_b']     = limit_mms_time_range(start_time_brst,stop_time_brst,fpi_brst_dis_dist_dict)
        prj_files[obs]['imoms_b']     = limit_mms_time_range(start_time_brst,stop_time_brst,fpi_brst_dis_moms_dict)
        prj_files[obs]['edist_f']     = limit_mms_time_range(start_time_fast,stop_time_fast,fpi_fast_des_dist_dict)
        prj_files[obs]['emoms_f']     = limit_mms_time_range(start_time_fast,stop_time_fast,fpi_fast_des_moms_dict)
        prj_files[obs]['idist_f']     = limit_mms_time_range(start_time_fast,stop_time_fast,fpi_fast_dis_dist_dict)
        prj_files[obs]['imoms_f']     = limit_mms_time_range(start_time_fast,stop_time_fast,fpi_fast_dis_moms_dict)
        prj_files[obs]['hpca_ion_b']  = limit_mms_time_range(start_time_brst,stop_time_brst,hpca_brst_ion_dict)
        prj_files[obs]['hpca_moms_b'] = limit_mms_time_range(start_time_brst,stop_time_brst,hpca_brst_moments_dict)
        prj_files[obs]['hpca_ion_f']  = limit_mms_time_range(start_time_brst,stop_time_brst,hpca_srvy_ion_dict)
        prj_files[obs]['hpca_moms_f'] = limit_mms_time_range(start_time_brst,stop_time_brst,hpca_srvy_moments_dict)
        prj_files[obs]['mec_s']       = limit_mms_time_range(start_time_srvy,stop_time_srvy,mec_srvy_epht89d_dict)
        prj_files[obs]['ecnts_b']     = limit_mms_time_range(start_time_brst,stop_time_brst,fpi_brst_des_cnts_dict)        
    
    scanned_files.close()
    #import pdb; pdb.set_trace()      
    try:
        prj_file_list = open(project_filelist_filename,'w')
        pickle.dump(prj_files,prj_file_list)
        prj_file_list.close()        
    except:
        workflow_log.write("Can't save the local project file list!\n")
        workflow_log.write("******************************************************************************\n")            
        return False    
        
    now = dt.datetime.now()
    workflow_log.write("Done with filtering at %s!\n" % (now,))
    workflow_log.write("******************************************************************************\n")    
    workflow_log.close()
    return prj_files

###############################################################################    
def load_filelist_for_prj(project_filelist_filename):
    """
    File to restore all the data structures for projects
    
    project_filelist_filename - name of the pickle file with a list of filenames germane to the project
    """
    try:
        scanned_files          = open(project_filelist_filename,"r")
    except:
        print "Can't load the scraped file list!"
        return False

    prj_files = pickle.load(scanned_files)
    
    scanned_files.close()
    
    return prj_files   
    
###############################################################################
def inventory_mms_file(filename,base_dir,download_dir):
    """inventory_mms_file is designed to take a filename and place the 
       corresponding file into the appropriate directory 
       (creating the directory if required).
       
       It's categorization is based on the MMS science data conventions.  
       The terminal and source directories are always
       relative to variables base_dir (and curr_dir) and download_dir that
       are input (except curr_dir, which is base_dir dressed locally)

       Arguments:
          filename:  base filename
          base_dir:  directory where to find the file downloaded from the server
          my_dict:     a output product from record_my_mms_files, record_my_mms_l2_files, or
                       construct_mms_file_dict

       Returns:
           True or False
       
       Example use:  
           inventory_mms_file(file_item,inventory_base_dir,download_dir)
              
       Note:  N/A
       
       
       """
    
    #pieces - an array that contains the separate descriptors of a science file
    pieces    = filename.split('_')
    mode      = pieces[2]
    epoch_str = pieces[-2]
    year      = epoch_str[0:4]
    month     = epoch_str[4:6]
    day       = epoch_str[6:8]

    #get rid of the epoch_str and file version
    del pieces[-1]
    del pieces[-1]
    
    #put year and month into pieces
    pieces.append(year)
    pieces.append(month)
    
    if mode == 'brst':
        pieces.append(day)
    
    #now that all the pieces are available start making the directories
    curr_dir = base_dir[0:-1]
    #print pieces
    for p in pieces:
        #print curr_dir
        curr_dir = curr_dir+'/'+p
        if os.path.isdir(curr_dir):
            pass
        else:
            os.mkdir(curr_dir)
            
    #finally move the file into that directory
    try:
        shutil.move(download_dir+filename,curr_dir+'/'+filename)
        return True
    except:
        return False

###############################################################################
def download_file(request_str,local_dir):
    """
    A helper function to download a file from the SDC.  It is based
    on the code found at 
    https://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    """
    
    local_filename =local_dir+request_str.split('=')[-1]
    # NOTE the stream=True parameter
    r = requests.get(request_str, auth=schiff_auth, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                #f.flush() commented by recommendation from J.F.Sebastian
    f.close()
    return local_filename

###############################################################################
def SDC_data_summary(data_sets,start_date,end_date,fh='None'):
    """
    A helper function to summarize the number of files from the SDC.
    """
    
    for data_set in data_sets:
        instr          = data_set[0]
        data_rate_mode = data_set[1]
        descriptor     = data_set[2]
        my_sdc_request = pre_string+\
                         start_pre+start_date+cnkt+\
                         end_pre+end_date+cnkt+\
                         instrument_pre+instr+cnkt+\
                         data_rate_mode_pre+data_rate_mode+cnkt+\
                         data_level
        if len(descriptor) > 0:
            my_sdc_request = my_sdc_request+cnkt+descriptor_pre+descriptor
        csv_string = requests.get(my_sdc_request, auth=schiff_auth).content
        csv_list   = csv_string.split(',')
        file_list  = [string.split('/')[-1] for string in csv_list]
        if fh == 'None':
            print data_set
            print len(file_list)
        else:
            fh.write("\n%s\n%s\n" % (data_set,len(file_list)))

###############################################################################
def retrieve_SDC_data(data_sets,download_dir,inventory_base_dir,start_date,end_date):
    """
    A helper function to retrieve data from the SDC.
    """
    report_file = open(download_dir+'SDC_retrieval_log.dat','a')
    now         = dt.datetime.now()
    #report the preliminary information)
    report_file.write("\n\n******************************************************************************\n")
    report_file.write("SDC Retrieval Log:  %s\n" % (now,))
    report_file.write("Summary of SDC files over range from %s to %s\n" % (start_date,end_date))
    
    #import pdb; pdb.set_trace()   
    for data_set in data_sets:
        instr          = data_set[0]
        data_rate_mode = data_set[1]
        descriptor     = data_set[2]
        my_sdc_request = pre_string+\
                         start_pre+start_date+cnkt+\
                         end_pre+end_date+cnkt+\
                         instrument_pre+instr+cnkt+\
                         data_rate_mode_pre+data_rate_mode+cnkt+\
                         data_level
        if len(descriptor) > 0:
            my_sdc_request = my_sdc_request+cnkt+descriptor_pre+descriptor
        SDC_data_summary([data_set],start_date,end_date,report_file)            
        csv_string = requests.get(my_sdc_request, auth=schiff_auth).content    
        csv_list   = csv_string.rstrip().split(',')
        file_list  = [string.split('/')[-1] for string in csv_list]
        for file_item in file_list:
            my_download_request = download_pre+file_item
            report_file.write("\nDownloading %s\n" % (file_item,))
            attempt_count = 0
            while attempt_count < 10:
                try:
                    file_name = download_file(my_download_request,download_dir)
                    attempt_count = 10
                    report_file.write("Success for %s!\n" %(file_item,))
                    try:
                        inventory_mms_file(file_item,inventory_base_dir,download_dir)
                        report_file.write("Inventoried %s!\n" % (file_item,))
                    except:
                        report_file.write("Couldn't move %s.\n" % (file_item,))                    
                except:
                    attempt_count += 1
                    report_file.write("Attempt %s for %s\n" % (file_item, len(file_item)))

                
    now = dt.datetime.now()
    report_file.write("\nCompleted retrieval from the SDC at %s - bye\n" % (now,))
    report_file.write("******************************************************************************\n\n\n")    
    report_file.close()
    return True
                
###############################################################################

###############################################################################

###############################################################################

###############################################################################

###############################################################################

###############################################################################

        