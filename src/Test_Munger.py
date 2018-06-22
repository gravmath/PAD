import cPickle as pickle
import numpy as np
import scipy.interpolate as interp
import sys

import Analyzer
import Grapher
import Munger
import PAD
import Plasma_Plotter
import Scraper

obs = 'mms1'
obs_files = pickle.load(open('c:/Users/cschiff/Documents/github/Plasma_Projects/Jim_Burch_Science_Unit_Test/Jim_Burch_Science_Unit_Test.data_local','r')) 
m1f = obs_files[obs]

dce_munge   = Munger.make_munge_via_translation(obs,'dce',Munger.dce_delta,m1f['dce_b'],Munger.dce_translation)
print 
fgm_munge   = Munger.make_munge_via_translation(obs,'fgm',Munger.fgm_delta,m1f['fgm_b'],Munger.fgm_translation) 
print 
emoms_munge = Munger.make_munge_via_translation(obs,'emoms',Munger.des_delta,m1f['emoms_b'],Munger.emoms_translation)
print 
imoms_munge = Munger.make_munge_via_translation(obs,'imoms',Munger.dis_delta,m1f['imoms_b'],Munger.imoms_translation)
print 
edist_munge = Munger.make_munge_via_translation(obs,'edist',Munger.des_delta,m1f['edist_b'],Munger.edist_translation)
print 
idist_munge = Munger.make_munge_via_translation(obs,'idist',Munger.dis_delta,m1f['idist_b'],Munger.idist_translation)
print 
bpsd_munge  = Munger.make_munge_via_translation(obs,'bpsd',0,m1f['bpsd_f'],Munger.bpsd_translation)
print 
epsd_munge  = Munger.make_munge_via_translation(obs,'epsd',0,m1f['epsd_f'],Munger.epsd_translation)
print 
mec_munge   = Munger.make_munge_via_translation(obs,'mec',0,m1f['mec_s'],Munger.mec_translation)
print 
scpot_munge = Munger.make_munge_via_translation(obs,'scpot',0,m1f['scpot_f'],Munger.scpot_translation)

efgm_munge   = Munger.adapt_munge_to_munge(fgm_munge,emoms_munge)