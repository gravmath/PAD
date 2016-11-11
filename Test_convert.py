#A simple comment - that can be removed
import Convert
import numpy as np
#import pdb

epoch_dict          = {}
epoch_dict['year']  = 2015
epoch_dict['month'] = 8
epoch_dict['day']   = 20
epoch_dict['hour']  = 12
epoch_dict['min']   = 0
epoch_dict['sec']   = 0

B_GSE = np.array([25.54237747,-243.367813,637.019104])
r_GSE = np.array([-7.1,-5.5,-2.0])

#pdb.set_trace()
print Convert.convert_GSE_to_GSM(epoch_dict).dot(r_GSE)