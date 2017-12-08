import re
import os

def find_debug_CDFs(obs,ver,basedir,year,month):
    file_pattern = re.compile('%s_fpi_fast_l2_des-debug_\d{14}_v%s.cdf' % (obs,ver))
    
    targ_dir = (basedir+'%s/%s/') % (obs,year,str(month).zfill(2))
    print 'scraping in %s' % targ_dir
 
    filtered_files = []
    for root, dirs, files in os.walk(targ_dir):
        for file in files:
            if file_pattern.match(file):
                filtered_files.append(targ_dir+file)
    return filtered_files
    
def get_debug_CDF_epoch(filename):
    basename = os.path.basename(filename)
    return basename.split('_')[5]    
    
def get_debug_CDF_ver(filename):
    basename = os.path.basename(filename)
    return (basename.split('_')[6])[:-4]
    
def cull_list(file_list):
    #sort list and assume that files with the same epoch now appear
    #in one contiguous block
    sorted_file_list = sorted(file_list)
    
    culled_list = []
    curr_file   = sorted_file_list[0]
    curr_epoch  = get_debug_CDF_epoch(curr_file)
    curr_ver    = get_debug_CDF_ver(curr_file)
    for sf in sorted_file_list:
        sf_epoch = get_debug_CDF_epoch(sf)
        sf_ver   = get_debug_CDF_ver(sf)
        if sf_epoch == curr_epoch:
            if sf_ver > curr_ver:
                curr_epoch = sf_epoch
                curr_ver   = sf_ver
                curr_file  = sf
        if sf_epoch > curr_epoch:
            #store the latest file and report
            culled_list.append(curr_file)
              #print 'Just added %s to the database' % curr_file
              #print 'Database size is %s' % str(len(culled_list))
            #reset values
            curr_epoch = sf_epoch
            curr_ver   = sf_ver
            curr_file  = sf
            
    last_file  = sorted_file_list[-1]
    last_epoch = get_debug_CDF_epoch(last_file)
    last_ver   = get_debug_CDF_ver(last_file)
    
    if last_epoch in culled_list:
        match_file = culled_list.index(last_epoch)
        match_ver  = get_debug_CDF_ver(match_file)
        if last_ver > match_ver:
            culled_list.append(last_file)
    else:
        culled_list.append(last_file)
    
    return culled_list
