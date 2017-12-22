#study data
obs     = 'mms1'
year    = '2017'
month   = '06'
day     = '08'
basedir = '/fpiprd1/fpishare/Conrad/Yuggoth/'
#basedir = 'c:/Yuggoth/'

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
sys.path.append('c:/Users/Conrad/Documents/GitHub/PAD/src/')
sys.path.append('/home/cschiff/PAD/src/')
import Burst_Munger as munge
import Grapher
import PAD

epoch_strings = ['20170608133933',
                 '20170608134133',
                 '20170608134403',
                 '20170608134623',
                 '20170608134853',
                 '20170608135003',
                 '20170608135133',
                 '20170608135353',
                 '20170608135623',
                 '20170608135803',
                 '20170608135943',
                 '20170608140143',
                 '20170608140353',
                 '20170608140603',
                 '20170608140813',
                 '20170608141033',
                 '20170608141243',
                 '20170608141453',
                 '20170608141703',
                 '20170608141923',
                 '20170608142133',
                 '20170608142303',
                 '20170608143333',
                 '20170608143453']

munge.config_directories(basedir,obs,year,month,day)

Be   = munge.munge_moms(epoch_strings,'des')
