'''
Created on Dec 10, 2012

@author: jwe

 NAME:
       MEANCLIP

 PURPOSE:
       Computes an iteratively sigma-clipped mean on a data set
 EXPLANATION:
       Clipping is done about median, but mean is returned.
       Called by SKYADJ_CUBE

 CATEGORY:
       Statistics

 CALLING SEQUENCE:
       MEANCLIP, Data, Mean, [ Sigma, SUBS =
              CLIPSIG=, MAXITER=, CONVERGE_NUM=, /VERBOSE, /DOUBLE ]

 INPUT POSITIONAL PARAMETERS:
       Data:     Input data, any numeric array
       
 OUTPUT POSITIONAL PARAMETERS:
       Mean:     N-sigma clipped mean.
       Sigma:    Standard deviation of remaining pixels.

 INPUT KEYWORD PARAMETERS:
       CLIPSIG:  Number of sigma at which to clip.  Default=3
       MAXITER:  Ceiling on number of clipping iterations.  Default=5
       CONVERGE_NUM:  If the proportion of rejected pixels is less
           than this fraction, the iterations stop.  Default=0.02, i.e.,
           iteration stops if fewer than 2% of pixels excluded.
       /VERBOSE:  Set this flag to get messages.
       /DOUBLE - if set then perform all computations in double precision.
                 Otherwise double precision is used only if the input
                 data is double
 OUTPUT KEYWORD PARAMETER:
       SUBS:     Subscript array for pixels finally used.


 MODIFICATION HISTORY:
       Written by:     RSH, RITSS, 21 Oct 98
       20 Jan 99 - Added SUBS, fixed misplaced paren on float call, 
                   improved doc.  RSH
       Nov 2005   Added /DOUBLE keyword, check if all pixels are removed  
                  by clipping W. Landsman 

'''
from numpy import median, isnan, sqrt, var, mean

def meanclip( image, clipsig=3, maxiter=5, 
              converge_num=0.02, verbose=False, subs=None):

    skpix = image[~isnan(image)]
    ct = len(skpix)
    lastct=1
    niter=0
    while (abs(ct-lastct)/lastct >= converge_num) and (niter < maxiter) and (ct > 0):
        niter = niter + 1
        lastct = ct
        medval = median(skpix)
        sig = sqrt(var(skpix))
        subs = skpix[abs(skpix-medval) < clipsig*sig]
        ct = len(subs)
        if ct > 0:
            skpix = subs         
    mean = mean(skpix)
    sigma = sqrt(var(skpix))
    if verbose:
        print clipsig,'-sigma clipped mean'
        print 'Mean computed in ',niter,' iterations'
        print 'Mean = ',mean,',  sigma = ',sigma
    return([mean,sigma])
