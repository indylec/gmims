"""
This implements the Theil-Sen linear regression estimator for 2d data points.
The jist of it is:
It returns the median all computed slope value between pairs (x_i, y_i), (x_j, y_j), (x_i > x_j)
where slope = (y_i - y_j)/(x_i - x_j)


Very robust to outliers.

"""
import numpy as np
import bottleneck #very fast searching and sorting written in Cython.
import itertools
from scipy.misc import comb
#import matplotlib.pyplot as plt

def theil_sen(x,y, sample= "auto", n_samples = 1e7):
    """
    Computes the Theil-Sen estimator for 2d data.
    parameters:
        x: 1-d np array, the control variate
        y: 1-d np.array, the ind variate.
        sample: if n>100, the performance can be worse, so we sample n_samples.
                Set to False to not sample.
        n_samples: how many points to sample.
    
    This complexity is O(n**2), which can be poor for large n. We will perform a sampling
    of data points to get an unbiased, but larger variance estimator. 
    The sampling will be done by picking two points at random, and computing the slope,
    up to n_samples times.
    
    """
    assert x.shape[0] == y.shape[0], "x and y must be the same shape."
    n = x.shape[0]
    
    if n < 100 or not sample:
        ix = np.argsort( x )
        slopes = np.empty( int(n*(n-1)*0.5) )
       # print '...calculating slopes...'
        for c, pair in enumerate(itertools.combinations( range(n),2 ) ): #it creates range(n) =( 
            i,j = ix[pair[0]], ix[pair[1]]
            slopes[c] = slope( x[i], x[j], y[i],y[j] )
        #print 'slope min and max are:',np.amin(slopes),np.amax(slopes)
        rank_up=int(0.5*int(comb(n,2))+1.96 * np.sqrt(n*(n-1.)*(2.*n-5.)/18.)+1.)
        rank_low=int(0.5*int(comb(n,2))- 1.96 * np.sqrt(n*(n-1.)*(2.*n-5.)/18.))
        #c95=np.percentile(slopes,(5,95))
        slopes_sort=np.sort(slopes)
        slope_up=slopes_sort[rank_up]
        slope_low=slopes_sort[rank_low]
    else:
        i1 = np.random.randint(0,int(n), int(n_samples))
        i2 = np.random.randint(0,int(n), int(n_samples))
        #print '...checking for unwanted zeros...'
        zero_check=np.where(np.abs((x[i1]-x[i2])) != 0)
        i1=i1[zero_check]
        i2=i2[zero_check]
        #print '...calculating slopes...'
        slopes = slope( x[i1], x[i2], y[i1], y[i2] )
        #print 'slope min and max are:',np.amin(slopes),np.amax(slopes)
        #c95=np.percentile(slopes,(5,95))
        #pdb.set_trace()
    
    slope_ = bottleneck.nanmedian( slopes )
    #np.savetxt('slope_one_sample_test.txt',slopes)
    #print '...done! Now finding intercepts...'
    #find the optimal b as the median of y_i - slope*x_i
    intercepts = np.empty( n )
    for c in xrange(n):
        intercepts[c] = y[c] - slope_*x[c]
    #np.savetxt('intercept_one_sample_test.txt',slopes)

    
   
    #c95i=np.percentile(intercepts,(5,95))
   
    intercept_ = bottleneck.median( intercepts )

    return np.array( [slope_,intercept_,slope_up, slope_low] )
        
        
        
def slope( x_1, x_2, y_1, y_2):
    return (1 - 2*(x_1>x_2) )*( (y_2 - y_1)/np.abs((x_2-x_1)) )
    
    
    
    
if __name__=="__main__":
    x = np.asarray( [ 0.0000, 0.2987, 0.4648, 0.5762, 0.8386 ] ) 
    y = np.asarray( [ 56751, 57037, 56979, 57074, 57422 ] ) 
    print theil_sen( x, y )
    
