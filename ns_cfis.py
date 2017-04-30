from astropy.io import fits
from matplotlib import pyplot as plt
from numpy import *
from sklearn.ensemble import RandomForestClassifier as RF
from sklearn.ensemble import IsolationForest as IF
from sklearn.metrics import confusion_matrix, roc_curve
from sklearn.model_selection import train_test_split
import pandas as pd
import seaborn as sns
import os, sys

class fits_cat:
    def __init__(self, fitsfile):
        self.fits_name = fitsfile
        self.fits_file = fits.open(self.fits_name)
        self.header = self.fits_file[1].header
        self.catalog = self.fits_file[1].data
        self.stationary0 = self.catalog[self.catalog['MATCHES'] > 0]
        self.stationary1 = self.catalog[self.catalog['MATCHES'] > 1]
        self.stationary2 = self.catalog[self.catalog['MATCHES'] > 2]
        self.ns0 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] >1)]
        self.ns1 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] >2)]
        self.real_feature = self.feature(self.stationary2)
        self.ns_feature = self.feature(self.ns1)
    def plot(self):
        #plt.xlim(165.5, 168.0)
        #plt.ylim(31.5, 32.0)
        #plt.xlim(165, 168.0)
        #plt.plot(self.stationary0['X_WORLD'], self.stationary0['Y_WORLD'], '.', label='MATCH >0')
        #plt.plot(self.stationary1['X_WORLD'], self.stationary1['Y_WORLD'], '.', label='MATCH >1')
        #plt.plot(self.ns['X_WORLD'], self.ns['Y_WORLD'], '.', label='NS, OVERLAPS >1')
        #plt.plot(self.ns1['X_WORLD'], self.ns1['Y_WORLD'], '.', label='NS, OVERLAPS >2')
        plt.plot(self.catalog['X_WORLD'], self.catalog['Y_WORLD'], marker = '.', label='ALL')
        plt.legend(loc=2)
    def feature(self, cat):
        return array([cat['MAG_ISO'], cat['MAGERR_ISO'], cat['MAG_ISOCOR'], cat['MAGERR_ISOCOR'], 
                       cat['MAG_APER'], cat['MAGERR_APER'], cat['MAG_AUTO'], cat['MAGERR_AUTO'],
                       cat['MAG_PETRO'], cat['MAGERR_PETRO'], cat['MAG_WIN'], cat['MAGERR_WIN'],
                       cat['SNR_WIN'], cat['MAG_SOMFIT'], cat['MAGERR_SOMFIT'], cat['KRON_RADIUS'],
                       cat['PETRO_RADIUS'], cat['BACKGROUND'], cat['ISOAREA_IMAGE'], cat['ISOAREAF_IMAGE'],
                       cat['A_IMAGE'], cat['B_IMAGE'], cat['THETA_IMAGE'], cat['ERRX2_IMAGE'],
                       cat['ERRY2_IMAGE'], cat['ERRXY_IMAGE'], cat['ERRCXX_IMAGE'], cat['ERRCYY_IMAGE'],
                       cat['ERRCXY_IMAGE'], cat['CXXWIN_IMAGE'], cat['CYYWIN_IMAGE'], cat['CXYWIN_IMAGE'],
                       cat['AWIN_IMAGE'], cat['BWIN_IMAGE'], cat['ERRX2WIN_IMAGE'], cat['ERRY2WIN_IMAGE'],
                       cat['ERRXYWIN_IMAGE'], cat['ERRCXXWIN_IMAGE'], cat['ERRCYYWIN_IMAGE'], cat['ERRCXYWIN_IMAGE'],
                       cat['MU_THRESHOLD'], cat['MU_MAX'], cat['ISO0'], cat['ISO1'], cat['ISO2'], cat['ISO3'],
                       cat['ISO4'], cat['ISO5'], cat['ISO6'], cat['ISO7'], cat['FWHM_IMAGE'], cat['ELONGATION'],
                       cat['ELLIPTICITY'], cat['POLAR_IMAGE'], cat['POLARWIN_IMAGE'], cat['FLUX_RADIUS'], 
                       cat['FWHMPSF_IMAGE'], cat['MAG_PSF'], cat['MAGERR_PSF'], cat['CHI2_PSF'],
                       cat['ERRCXXPSF_IMAGE'], cat['ERRCYYPSF_IMAGE'], cat['ERRCXYPSF_IMAGE'],
                       cat['ERRAPSF_IMAGE'], cat['ERRBPSF_IMAGE'], cat['ERRAPSF_WORLD'], cat['ERRBPSF_WORLD'],
                       cat['ERRTHETAPSF_SKY']]).T
        
class IF_real_bogus:
    def __init__(self, feature):
        self.feature = feature
        #self.Class = []
        #[self.Class.append('real') for i in range(len(self.feature))]
        #self.trainning, self.test, self.Class_train, self.Class_test = train_test_split(self.real_feature,\
        #                        self.real_Class, test_size=0.4, random_state=0)
        #print self.trainning     
    def train(self):
        self.IFmod = IF(n_estimators = 160)
        self.IFmod.fit(self.feature)
        
    def validation(self):
        #score = self.IFmod.decision_function(self.feature)
        result = self.IFmod.predict(self.feature)
        self.normal = (result == 1)
        self.abnormal = (result == -1)
        print "total:{0}, normal:{1}, abnormal:{2}, total/normal:{3}".format(len(self.feature), 
            	self.normal.sum(), self.abnormal.sum(), self.normal.sum()/float(len(self.feature)))
                    
def build_model(fits_tab):
	sources = fits_cat(fits_tab)
	print len(sources.stationary2), len(sources.ns1)
	real_bogus = IF_real_bogus(sources.real_feature)
	real_bogus.train()
	real_bogus.validation()
	return real_bogus.IFmod

def ns_cat(fits_tab):
	#try:
	sources = fits_cat(fits_tab)
	print fits_tab
	real_bogus = IF_real_bogus(sources.ns_feature)
	real_bogus.IFmod = RB_model
	real_bogus.validation()
	ns_cat = sources.ns1[real_bogus.normal]
	fits.writeto('{}.ns'.format(fits_tab), ns_cat, header=sources.header)
	#except ValueError:
		#print 'error: {0}'.format(fits_tab)
	
def main():
	global RB_model
	workdir = '/sciproc/disk2/cfis/catalog_20170419/'
	fits_list = filter(lambda x: x.endswith('fits'), os.listdir(workdir))
	fits_list = [workdir + i for i in fits_list]
	fits_list = ['/sciproc/disk2/cfis/catalog_20170419/HPX_02559_RA_178.6_DEC_+35.7_cat.fits']
	training_fits = 'HPX_02937_RA_160.3_DEC_+31.4_cat.fits'
	RB_model = build_model(workdir + training_fits)
	map(ns_cat, fits_list)
	
if __name__ == '__main__':
	main()
	
	
	
	
	
	    
