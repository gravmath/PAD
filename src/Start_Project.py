import os
import pickle
import shutil
import sys
sys.path.append('../src')
import Scraper

###############################################################################
def inventory_file(filename,base_dir,download_dir):
    """inventory_file is designed to take a filename and place it into the
       appropriate directory (creating the directory if required) based on the
       MMS science conventions.  The terminal and source directories are always
       relative to variables base_dir (and curr_dir) and download_dir that
       are input (except curr_dir, which is base_dir dressed locally)"""
    
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
    curr_dir = base_dir
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
        shutil.move(download_dir+'/'+filename,curr_dir+'/'+filename)
        return True
    except:
        return False
        
###############################################################################       
#simple, quick and dirty script to start a project

prj_settings = open(sys.argv[1],'r')
glob         = prj_settings.readlines()
prj_settings.close()

prj_dict = {}

for j in range(1,len(glob)):
    namelist = [q.strip() for q in glob[j].split('=')]
    prj_dict[namelist[0]] = namelist[1]

base_dir               = prj_dict['base_dir']
download_dir           = prj_dict['download_dir']
local_dir              = prj_dict['local_dir']
scan_or_load           = prj_dict['scan_or_load']
observatories          = ['mms1','mms2','mms3','mms4']

my_files = {'mms1':{},'mms2':{},'mms3':{},'mms4':{}}

my_files['mms1'] = {'bpsd_f':[],
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
                    'mec_s':[]
                    }

my_files['mms2'] = {'bpsd_f':[],
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
                    'mec_s':[]
                    }

my_files['mms3'] = {'bpsd_f':[],
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
                    'mec_s':[]
                    }

my_files['mms4'] = {'bpsd_f':[],
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
                    'mec_s':[]
                    }
                    
start_time_brst = prj_dict['start_time_brst']
stop_time_brst  = prj_dict['stop_time_brst']
start_time_fast = prj_dict['start_time_fast']
stop_time_fast  = prj_dict['stop_time_fast'] 
start_time_srvy = prj_dict['start_time_srvy']
stop_time_srvy  = prj_dict['stop_time_srvy']
 
for obs in observatories:
    if prj_dict['file source'] == 'server':
        scanned_files_filename = 'c:/Users/cschiff/Documents/github/scanned_files_%s.fpiprd1' % (obs,)
    else:
        scanned_files_filename = 'c:/Users/cschiff/Documents/github/scanned_files_%s.local' % (obs,)

    if prj_dict['scan_or_load'] == 'scan':
        print 'Preparing to scrape together files for %s from the %s drive and save to %s' % (obs, prj_dict['file source'],scanned_files_filename)
        bpsd_fast_dict         = Scraper.get_my_l2_files(obs,'dsp', 'fast','bpsd',    base_dir,prj_dict['file source'])
        epsd_fast_dict         = Scraper.get_my_l2_files(obs,'dsp', 'fast','epsd',    base_dir,prj_dict['file source'])
        dce_brst_dict          = Scraper.get_my_l2_files(obs,'edp', 'brst','dce',     base_dir,prj_dict['file source'])
        scpot_fast_dict        = Scraper.get_my_l2_files(obs,'edp', 'fast','scpot',   base_dir,prj_dict['file source'])
        fgm_brst_dict          = Scraper.get_my_l2_files(obs,'fgm', 'brst','',        base_dir,prj_dict['file source'])
        fgm_srvy_dict          = Scraper.get_my_l2_files(obs,'fgm', 'srvy','',        base_dir,prj_dict['file source'])        
        fpi_brst_des_dist_dict = Scraper.get_my_l2_files(obs,'fpi', 'brst','des-dist',base_dir,prj_dict['file source'])
        fpi_brst_des_moms_dict = Scraper.get_my_l2_files(obs,'fpi', 'brst','des-moms',base_dir,prj_dict['file source'])
        fpi_brst_dis_dist_dict = Scraper.get_my_l2_files(obs,'fpi', 'brst','dis-dist',base_dir,prj_dict['file source'])
        fpi_brst_dis_moms_dict = Scraper.get_my_l2_files(obs,'fpi', 'brst','dis-moms',base_dir,prj_dict['file source'])
        fpi_fast_des_dist_dict = Scraper.get_my_l2_files(obs,'fpi', 'fast','des-dist',base_dir,prj_dict['file source'])
        fpi_fast_des_moms_dict = Scraper.get_my_l2_files(obs,'fpi', 'fast','des-moms',base_dir,prj_dict['file source'])
        fpi_fast_dis_dist_dict = Scraper.get_my_l2_files(obs,'fpi', 'fast','dis-dist',base_dir,prj_dict['file source'])
        fpi_fast_dis_moms_dict = Scraper.get_my_l2_files(obs,'fpi', 'fast','dis-moms',base_dir,prj_dict['file source'])
        hpca_brst_ion_dict     = Scraper.get_my_l2_files(obs,'hpca','brst','ion',     base_dir,prj_dict['file source'])
        hpca_brst_moments_dict = Scraper.get_my_l2_files(obs,'hpca','brst','moments', base_dir,prj_dict['file source'])
        mec_srvy_epht89d_dict  = Scraper.get_my_l2_files(obs,'mec', 'srvy','epht89d', base_dir,prj_dict['file source'])
        scanned_files = open(scanned_files_filename, "w")
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
        pickle.dump(mec_srvy_epht89d_dict,scanned_files)
        scanned_files.close()
    else:
        print 'Reading in scraped files for %s on the %s drive from %s' % (obs, prj_dict['file source'],scanned_files_filename)        
        scanned_files          = open(scanned_files_filename,"r")
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
        mec_srvy_epht89d_dict  = pickle.load(scanned_files)
        scanned_files.close()

    my_files['%s' % (obs,)]['bpsd_f']      = Scraper.limit_time_range(start_time_srvy,stop_time_srvy,bpsd_fast_dict)
    my_files['%s' % (obs,)]['epsd_f']      = Scraper.limit_time_range(start_time_srvy,stop_time_srvy,epsd_fast_dict)
    my_files['%s' % (obs,)]['dce_b']       = Scraper.limit_time_range(start_time_brst,stop_time_brst,dce_brst_dict)
    my_files['%s' % (obs,)]['scpot_f']     = Scraper.limit_time_range(start_time_srvy,stop_time_srvy,scpot_fast_dict)
    my_files['%s' % (obs,)]['fgm_b']       = Scraper.limit_time_range(start_time_brst,stop_time_brst,fgm_brst_dict)
    my_files['%s' % (obs,)]['fgm_s']       = Scraper.limit_time_range(start_time_srvy,stop_time_srvy,fgm_srvy_dict)    
    my_files['%s' % (obs,)]['edist_b']     = Scraper.limit_time_range(start_time_brst,stop_time_brst,fpi_brst_des_dist_dict)
    my_files['%s' % (obs,)]['emoms_b']     = Scraper.limit_time_range(start_time_brst,stop_time_brst,fpi_brst_des_moms_dict)
    my_files['%s' % (obs,)]['idist_b']     = Scraper.limit_time_range(start_time_brst,stop_time_brst,fpi_brst_dis_dist_dict)
    my_files['%s' % (obs,)]['imoms_b']     = Scraper.limit_time_range(start_time_brst,stop_time_brst,fpi_brst_dis_moms_dict)
    my_files['%s' % (obs,)]['edist_f']     = Scraper.limit_time_range(start_time_fast,stop_time_fast,fpi_fast_des_dist_dict)
    my_files['%s' % (obs,)]['emoms_f']     = Scraper.limit_time_range(start_time_fast,stop_time_fast,fpi_fast_des_moms_dict)
    my_files['%s' % (obs,)]['idist_f']     = Scraper.limit_time_range(start_time_fast,stop_time_fast,fpi_fast_dis_dist_dict)
    my_files['%s' % (obs,)]['imoms_f']     = Scraper.limit_time_range(start_time_fast,stop_time_fast,fpi_fast_dis_moms_dict)
    my_files['%s' % (obs,)]['hpca_ion_b']  = Scraper.limit_time_range(start_time_brst,stop_time_brst,hpca_brst_ion_dict)
    my_files['%s' % (obs,)]['hpca_moms_b'] = Scraper.limit_time_range(start_time_brst,stop_time_brst,hpca_brst_moments_dict)
    my_files['%s' % (obs,)]['mec_s']       = Scraper.limit_time_range(start_time_srvy,stop_time_srvy,mec_srvy_epht89d_dict)
    
    
    if prj_dict['load_local'] == 'True':
        for k in my_files['%s' % (obs,)].keys():
            for f in my_files['%s' % (obs,)][k]:
                 bf = os.path.basename(f)
                 print bf
                 shutil.copy(f,download_dir+bf)      
                 inventory_file(bf,local_dir,download_dir)
                 
prj_file = open(prj_dict['prj_filename'],'w')
pickle.dump(my_files,prj_file)
prj_file.close()