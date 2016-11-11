import numpy  as np
import SunEph 
import equation_of_time

eclip_pole_GCI = np.array([0.0,-0.39777715575399,0.917482062146321])
dipole_ECEF    = np.array([-0.05752247,0.17291892,-0.98325491])

def convert_GCI_to_GSE(epoch_dict):
    S, V  = SunEph.CalcSun_Low(epoch_dict)
    X_hat = S/np.sqrt(S.dot(S))
    Y     = np.cross(eclip_pole_GCI,X_hat)
    Y_hat = Y/np.sqrt(Y.dot(Y))
    Z_hat = np.cross(X_hat,Y_hat)
    
    return np.vstack((X_hat,Y_hat,Z_hat))
    
def convert_ECEF_to_GCI(epoch_dict):
    gha   = equation_of_time.calculate_GMST(epoch_dict)
    cos_g = np.cos(gha)
    sin_g = np.sin(gha)
    
    return np.array([[cos_g,-sin_g,0],[sin_g,cos_g,0],[0,0,1]])
    
def convert_GSE_to_GSM(epoch_dict):
    A_ecef_to_gci = convert_ECEF_to_GCI(epoch_dict)
    A_gse_gci     = convert_GCI_to_GSE(epoch_dict)
    S, V          = SunEph.CalcSun_Low(epoch_dict)
    X_hat         = S/np.sqrt(S.dot(S))
    D             = np.dot(A_ecef_to_gci,dipole_ECEF)
    Y             = np.cross(S,D)
    Y_hat         = Y/np.sqrt(Y.dot(Y))
    Z_hat         = np.cross(X_hat,Y_hat)
    
    X_hat_GSE     = A_gse_gci[0,:]
    Y_hat_GSE     = A_gse_gci[1,:]
    Z_hat_GSE     = A_gse_gci[2,:]
    
    A00           = np.dot(X_hat_GSE,X_hat)
    A01           = np.dot(Y_hat_GSE,X_hat)
    A02           = np.dot(Z_hat_GSE,X_hat)
    A10           = np.dot(X_hat_GSE,Y_hat)
    A11           = np.dot(Y_hat_GSE,Y_hat)
    A12           = np.dot(Z_hat_GSE,Y_hat)
    A20           = np.dot(X_hat_GSE,Z_hat)
    A21           = np.dot(Y_hat_GSE,Z_hat)
    A22           = np.dot(Z_hat_GSE,Z_hat)

    return np.array([[A00,A01,A02],[A10,A11,A12],[A20,A21,A22]])

def convert_GSM_to_ECEF(epoch_dict):
    T_GCI_ECEF = convert_ECEF_to_GCI(epoch_dict)
    T_GSE_GCI  = convert_GCI_to_GSE(epoch_dict)
    T_GSM_GSE  = convert_GSE_to_GSM(epoch_dict)
    
    T_GSM_ECEF = T_GSM_GSE.dot(T_GSE_GCI.dot(T_GCI_ECEF))
    
    return T_GSM_ECEF.transpose()