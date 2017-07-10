import datetime          as dt
import matplotlib        as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy             as np
import sqlite3
import sys
import spacepy.pycdf     as pycdf

sys.path.append('c:/users/cschiff/Documents/GitHub/PAD/')
import Grapher
import Plasma_Plotter

sqlite_file = 'c:/users/cschiff/Documents/GitHub/PAD/instrument_data_db.sqlite'
conn = sqlite3.connect(sqlite_file)
cursor = conn.cursor()
fpi_prd1 = 'x:'

obs         = 'mms1'
year        = 2017
month       = 5
day         = 1
start_epoch = dt.datetime(year,month,day)

instrument = 'fpi'
instr_mode = 'fast'
ver        = '3.1.0'

mec_mode   = 'srvy'

for i in range(60):
    curr_epoch = start_epoch + dt.timedelta(days=i)
    curr_year  = curr_epoch.year
    curr_month = curr_epoch.month
    curr_day   = curr_epoch.day
    print curr_epoch
    
    #execute the des queries
    descriptor = 'des-debug'
    eresult = cursor.execute('Select filename from fpi_data where\
                                            obs        = "%s" and\
                                            mode       = "%s" and\
                                            descriptor = "%s" and\
                                            year       =  %s  and\
                                            month      =  %s  and\
                                            day        =  %s;' % 
                                           (obs,instr_mode,descriptor,curr_year,curr_month,curr_day))
    e = eresult.fetchall()

    descriptor = 'dis-debug'
    iresult = cursor.execute('Select filename from fpi_data where\
                                            obs        = "%s" and\
                                            mode       = "%s" and\
                                            descriptor = "%s" and\
                                            year       =  %s  and\
                                            month      =  %s  and\
                                            day        =  %s;' % 
                                            (obs,instr_mode,descriptor,curr_year,curr_month,curr_day))
    i = iresult.fetchall()
    
    descriptor = 'ephts04d'
    
    for N in range(len(e)):
        print N, len(e)
        #open the files
        efile  = pycdf.CDF(fpi_prd1+e[N][0])
        ifile  = pycdf.CDF(fpi_prd1+i[N][0])
        
        #import the data
        e_t = np.asarray(efile['Epoch'])
        i_t = np.asarray(ifile['Epoch'])
        e_E = np.asarray(efile['%s_des_energy_%s'%(obs,instr_mode)][0,:])
        i_E = np.asarray(ifile['%s_dis_energy_%s'%(obs,instr_mode)][0,:])
        e_n = np.asarray(efile['%s_des_numberdensity_%s' % (obs,instr_mode)])
        i_n = np.asarray(ifile['%s_dis_numberdensity_%s' % (obs,instr_mode)])
        e_Tperp = np.asarray(efile['%s_des_tempperp_%s' % (obs,instr_mode)])
        i_Tperp = np.asarray(ifile['%s_dis_tempperp_%s' % (obs,instr_mode)])
        e_V = np.asarray(efile['%s_des_bulkv_gse_%s' % (obs,instr_mode)])
        i_V = np.asarray(ifile['%s_dis_bulkv_gse_%s' % (obs,instr_mode)])
        e_B = np.asarray(efile['%s_des_b_gse_srvy' % obs])
        i_B = np.asarray(ifile['%s_dis_b_gse_srvy' % obs])
        e_omni = np.asarray(efile['%s_des_energyspectr_omni_%s' %(obs,instr_mode)])
        i_omni = np.asarray(ifile['%s_dis_energyspectr_omni_%s' %(obs,instr_mode)])
        e_scpot = np.asarray(efile['%s_des_scpot_max_%s' % (obs,instr_mode)])
        i_scpot = np.asarray(ifile['%s_dis_scpot_max_%s' % (obs,instr_mode)])
        #close files
        efile.close()
        ifile.close()
        
        fig,axs = plt.subplots(nrows=7,ncols=1,sharex=True,figsize=(12,20))
        Plasma_Plotter.make_density_panel (axs[0],obs,e_t,e_n,i_t,i_n)
        Plasma_Plotter.make_Tperp_panel   (axs[1],obs,e_t,e_Tperp,i_t,i_Tperp)
        Plasma_Plotter.make_eVvector_panel(axs[2],obs,e_t,e_V,)
        Plasma_Plotter.make_iVvector_panel(axs[3],obs,i_t,i_V)
        Plasma_Plotter.make_Et_panel      (fig,axs[4],obs,e_t,e_E,e_omni,e_scpot)
        Plasma_Plotter.make_Et_panel      (fig,axs[5],obs,i_t,i_E,i_omni,i_scpot)
        Plasma_Plotter.make_Bvector_panel (axs[6],obs,e_t,e_B)
        Plasma_Plotter.add_info_box(fig,'%s-%s-%s'%(day,month,year),[-0.02,0.017,0.1,0.1])
        Plasma_Plotter.pos_on_time_axis(fig,axs[6],cursor,fpi_prd1,obs,mec_mode,descriptor,year,month,day)
        fig.savefig(fpi_prd1+'/fpishare/Conrad/test/FastSurvey_plot_%s_%s_%s_%s_%s.png'%(obs,year,month,day,N))
        fig.clf()