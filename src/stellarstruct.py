'''
Created on 28.04.2013

source:
http://www.astro.princeton.edu/~gk/A403/hms11.f

----- calculates Hydrogen Main Sequence models, uses subroutines:
      sch envel core step rhs state opact nburn solve prout
----- kept as fortran files:
sch11 envel11 core11 step11 rhs1 state1 opact1 nburn1 solve11 prout1

INPUT (from the keyboard):
X, Z  = hydrogen and heavy element content     (by mass fraction)
fmlz  = log10 (stellar mass)             (M Sun)
tel   = log10 (effective temperature)    (K)
fll   = log10 (surface luminosity)       (L Sun)
tcl   = log10 (central temperature)      (K)
rhcl  = log10 (central density)          (g/cm**3)
iprint = print control for the integrations, iprint<0 implies none
nm    = number of stellar model to be calculated
dm    = step in log10 (mass) between stellar models
NOTICE: tel, fll, tcl, rhcl are guessed values of the boundary parameters
        to be improved by the subroutine sch

OUTPUT:
summary of results onto file "hms11.re1" and onto the screen
more detailes (when iprint>0) onto a file "hms11.re2" and onto the screen

'''

def sch(x, z, starm, tel, fll, tcl, rhcl, iprint):
    """
----- calculates chemically homogeneous stellar model in a thermal
----- equilibrium, i.e. a model on a zero age main sequence
----- fitting results of integrations of envelope and core at the
----- fitting mass, fmf=0.3*M

INPUT:
  x     = hydrogen content                   (by mass fraction)
  z     = heavy element content              (by mass fraction)
  starm = total stellar mass                 (M Sun)
  tel   = log10 (effective temperature)      (K)         (a guess)
  fll   = log10 (surface luminosity)         (L Sun)     (a guess)
  tcl   = log10 (central temperature)        (K)         (a guess)
  rhcl  = log10 (central density)            (g/cm**3)   (a guess)
  iout  = number of the output unit for storing results on the disk
  iprint = output control 

OUTPUT:
  tel, fll, tcl, rhcl = accurate values obtained from numerical iterations
  rl                  = log10 (surface radius)    (R Sun)
  iter                = number of iterations needed to obtain those results
     if iprint>0 then some results of integrations are stored on the disk

 NOTICE: this subroutine uses subroutines: envel, core, solve

---------- auxillary variables:
      zn    = nitrogen content for CNO burning   (by mass fraction)
      fmf   = mass at the fitting point             (M Sun)
      difac = the maximum differences allowed at the fitting point 
      delta = perturbations of the logarithms of boundary parameters
      itmax = the maximum number of iterations
      iter  = the iteration counter
      difmax = the maximum difference calculated at the fitting point
      delfit(i) = the calculated differences at the fitting point for
                  the four variables = xe(i)/xc(i)-1 , i=1,2,3,4
      xe(i)  = results of the envelope integrations at the fitting point
      xc(i)  = results of the core integrations at the fitting point
      deriv(i,j) = partial derivatives of the differences at the fitting
                   point with respect to the boundary parameters
      delce(i) = corrections added to the boundary parameters
"""

    #global sunm,sunl,sunr,sunt,cgrrad,cdpm,sunmr3,sunml,g,c,fmod

    #dimension xe(10),xc(10),delfit(4),delce(4),deriv(4,5),xe1(10),xe2(10),xc3(10),xc4(10)
    deriv = np.zeros([4,5])
    print '  Mass, log Te, log L, log Tc, log rhc, iter ='
    print '%12.4f, %8.4f, %8.4f, %8.4f, %8.4f, %5i' % (starm,tel,fll,tcl,rhcl,iprint)

# ----- set the constants for the fitting and the iterations ----
    zn = z*0.3
    fmf = starm*0.3
    difac = 0.0001
    delta = 0.0001
    itmax = 15
    iter = 0
    
    if iprint >=-1:
        print ' iter  log Te   log L  log Tc  log rhc    max del b.p.  max non-fit'

# ----- a new iteration ----------------------------------------------------
    while (iter <= itmax):
        iter += 1
    
        rl, xe = envel(tel,fll,x,z,zn,fmf,starm,iprint)
        xc = core(tcl,rhcl,x,z,zn,fmf,iprint)

# ----- see if the envelope and the core integrations fit 
        delfit = xe/xc - 1.0
        difmax = max(abs(dif))
        if(difmax < difac):
# ---- iterations converged, stellar envelope and core fit each other
            if (iprint < 0):
                print ' log Te, log L, log Tc, log rhc, log R, iter = %10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%5i' % (tel,fll,tcl,rhcl,rl,iter)
            return
          
# ----- the differences at the fitting point are too large to be acceptable,
# ----- calculate derivatives of differences at the fitting point with respect
# ----- to the boundary parameters: tel, fll, tcl, rhcl
        rl, xe1 = envel(tel+delta,fll,x,z,zn,fmf,starm,iprint)
        rl, xe2 = envel(tel,fll+delta,x,z,zn,fmf,starm,iprint)
        xc3 = core(tcl+delta,rhcl,x,z,zn,fmf,iprint)
        xc4 = core(tcl,rhcl+delta,x,z,zn,fmf,iprint)
        for i in range(4):
            deriv[i,0] = ((xe1[i]/xc[i]-1.0)-delfit[i])/delta
            deriv[i,1] = ((xe2[i]/xc[i]-1.0)-delfit[i])/delta
            deriv[i,2] = ((xe[i]/xc3[i]-1.0)-delfit[i])/delta
            deriv[i,3] = ((xe[i]/xc4[i]-1.0)-delfit[i])/delta
            deriv[i,4] = delfit[i]
    

# ----- calculate the corrections to the boundary parameters
        delce, facdel = solve(deriv)

# ----- add the corrections to the boundary parameters:
        tel = tel + delce[0]
        fll = fll + delce[1]
        tcl = tcl + delce[2]
        rhcl = rhcl + delce[3]

        if iprint >= -1:
            print ' %5i,%8.4f,%8.4f,%8.4f,%8.4f,%12.6f,%12.6f' % (iter,tel,fll,tcl,rhcl,facdel,difmax)
  
# ----- end of an iteration ------------------------------------------------
    print ' iterations have not converged '
    iter = -iter
    return np.array([tel,fll,tcl,rhcl,rl,iter])

def envel(tel,fll,x,z,zn,fmf,starm,iprint):
    """
    -- integrates stellar structure equations from the surface to the fitting
    -- point for a chemically homogeneous star
    
     INPUT:
       tel   = log10 (effective temperature)       (K)
       fll   = log10 (surface luminosity)          (L Sun)
       x     = hydrogen content                    (by mass fraction)
       z     = heavy element content               (by mass fraction)
       zn    = nitrogen content for CNO burning    (by mass fraction)
       fmf   = mass at the fitting point           (M Sun)
       starm = total stellar mass                  (M Sun)
       iout   = the number of the output unit for subroutine prout
       iprint = print control, one out of "iprint" steps to be printed with
                subroutine prout,  if iprint<0 then no print
    
     OUTPUT:
       rl    = log10 (surface radius)               (R Sun)
       xx(i), i=1,2,3,4,5,6,7,8  at the fitting point
       xx(1) = temperature                          (K)
       xx(2) = density                              (g/cm**3)
       xx(3) = radius                               (R Sun)
       xx(4) = luminosity                           (L Sun)
       xx(5) = mass                                 (M Sun)
       xx(6) = x = hydrogen content            (by mass fraction)
       xx(7) = z = heavy element content       (by mass fraction)
       xx(8) = zn = nitrogen content for CNO burning (by mass fraction)
    
     NOTICE: subroutine envel1 uses:
       subroutine step       to make integration steps 
       subroutine prout      to output results onto the terminal and unit "iout"
    """
    #global sunm,sunl,sunr,sunt,cgrrad,cdpm,sunmr3,sunml,g,c,fmod
      #dimension xx(10)

#--- calculate the conditions at the photosphere: temperature (K),
#--- luminosity (L Sun), and radius (R Sun)
    tph = 10.0**tel
    flph = 10.0**fll
    rph = np.sqrt(flph)/(tph/sunt)**2
    rl = np.log10(rph)

#--- set all the variables at the outer boundary, at a very small density,
#--- i.e. at a very small optical depth, where temperature is lower than
#--- the photospheric by a factor 2**0.25;  assume the atmosphere is very
#--- thin geometrically, so the radius is the same as at the photosphere.
    xx=np.array([tph/2.0**0.25, 1.0e-12, rph, flph, starm, x, z, zn])

    if (iprint>0): prout(0,xx,-1)
    if (iprint>0): prout(0,xx,1)

    kmax = 10
    ip = iprint

#----- integrate the structure equations inward to the fitting point
    for k in range(kmax):
        xx = step(xx, fmf)
        ip -= 1
        if ip!=0:
            
            if(abs(xx[4]/fmf-1.0)<0.00001):
                # ------ integrations completed -----
                if iprint>0: prout(-k,xx,1)
                return (rl, xx)
        prout(k,xx,1)
        ip = iprint
    print '  too many integration steps in subroutine envel'
    prout(k,xx,1)
    exit()

#------ integrations completed -----
    return (rl, xx)

def core(tcl,rhcl,x,z,zn,fmf,iprint):
    """
-- integrates stellar structure equations from the center to the fitting
-- point for a chemically homogeneous star

 INPUT:
     rhcl     = log10 (central density)     (g/cm**3)
     tcl      = log10 (central temperature)     (K)
     x        = hydrogen content (by mass fraction)
     z        = heavy element content (by mass fraction)
     zn       = nitrogen content for CNO burning (by mass fraction)
     fmf      = mass at the fitting point     (M Sun)
     iout     = the number of the output unit for subroutine prout
     iprint   = print control, one out of "iprint" steps to be printed with
                subroutine prout,   if iprint<0 then no printout

 OUTPUT:
     xx(i), i=1,2,3,4,5,6,7,8  at the fitting point
     xx(1)    = temperature        (K)
     xx(2)    = density            (g/cm**3)
     xx(3)    = radius             (R Sun)
     xx(4)    = luminosity         (L Sun) 
     xx(5)    = mass               (M Sun)
     xx(6)    = x = hydrogen content  (by mass fraction)
     xx(7)    = z = heavy element content  (by mass fraction)
     xx(8)    = zn = nitrogen content for CNO burning  (by mass fraction)

 NOTICE: subroutine core uses:
  subroutine nburn      to calculate luminosity at the inner boundary
  subroutine step       to make integration steps
  subroutine prout      to output results onto the terminal and unit "iout"

--------------- variables from common/const/ used in this subroutine:
     sunmr3   = ( solar mass ) / ( 4*pi*(solar radius)**3 )     (c.g.s.)
     sunml    = ( solar mass ) / ( solar luminosity )           (c.g.s.)

--------------- auxillary variables:
    tc        = central temperature    (K)
    rhc       = central density        (g/cm**3)
    rc        = radius of the central sphere / solar radius
    flc       = luminosity generated in the central sphere / solar luminosity
    fmc       = mass of the central sphere / solar mass
    epsx      = energy generation rate in hydrogen burning   (erg/g/sec)
"""
    #global sunm,sunl,sunr,sunt,cgrrad,cdpm,sunmr3,sunml,g,c,fmod
    #  dimension xx(10)

#--- calculate variables at the surface of the inner sphere, i.e. inner
#--- boundary conditions  (this section may be a separate subroutine center)
    tc = 10.0**tcl
    rhc = 10.0**rhcl
    fmc = fmf*0.0001
    epsx,dxt = nburn(rhc,tc,x,zn)
    rc = (fmc/rhc*sunmr3*3)**(1.0/3.0)
    flc = sunml*fmc*epsx

    xx = np.array([tc, rhc, rc, flc, fmc, x, z, zn])
# ----- all variables at the starting point are ready

    if iprint>0: prout(0,xx,-1)
    if iprint>0: prout(0,xx,1)

    kmax = 1000
    ip = iprint

# --- integrate the equations out to the fitting point:
    for k in range(kmax):
        xx = step(xx,fmf)
        ip = ip-1
        if ip!=0:
            if(abs(yy(5)/fmf-1.0)<0.00001):
#----- integrations completed ----
                if iprint>0: prout(-k,yy,1)
                return yy
    prout(k,xx,1)
    ip = iprint
    print ' too many integration steps in the subroutine core1 '
    prout(k,xx,1)
    exit()

def step(xx, fmf):
    """
    -- makes one integration step for a stellar core or envelope with a second
    -- order Runge-Kutta method, the steps are carried up to the point xx(5) = fmf
    
     INPUT:
       xx(1)  = t                temperature   (K)
       xx(2)  = rh               density       (g/cm**3)
       xx(3)  = r                radius        (R Sun)
       xx(4)  = flr              luminosity    (L Sun)
       xx(5)  = fmr              mass          (M Sun)     (independent variable)
       fmf              mass at the fitting point   (M Sun)
    
     OUTPUT:
       xx(i), i=1,2,3,4,5   at the end of integration step
    
     NOTICE:  subroutine step uses:
      subroutine rhs     to calculate right hand sides of differential equations
    
    ------------- auxillary variables:
          acc(i), i=1,2,3,4,5   = maximum step size in all variables
          h   = integration step size
          h2  = half of the integration step
          yy(i) = d xx(i) / d xx(nv) as calculated by subroutine rhs1
          xi(i) = values of xx(i) variables at the middle of the integration step
    """
    #dimension xx(10),yy(10),xi(10),acc(5)
    xi = np.zeros(8)
    acc = np.array([0.05, 0.15, 0.05, 0.2, 0.2])
    yy = rhs(xx)
    
# ------ estimate the integration step = h -----
    hi = min(abs(yy[:5]/acc[:5]/xx[:5]))
    
    h = 1.0/hi
    if xx[4] > fmf :
        h = -h
    if ((xx[4]-fmf) * (xx[4]+h-fmf)) < 0.0:
        h = fmf-xx[4]

    h2 = h*0.5
# ----- make the first half of the integration step -------
    xi[:5] = xx[:5] + h2*yy[:5]
    xi[5:8] = xx[5:8]
    
    yy = rhs(xi)

# ----- make the whole integration step -------
    xx[:5] = xx[:5] + h*yy[:5]
    #print 'step() = ', xx
    return xx

def rhs(xx):
    """
    -- calculates right hand sides of the stellar structure equations
    INPUT:
         xx(1)  = t                  (K)           temperature
         xx(2)  = rh                 (g/cm**3)     density
         xx(3)  = r                  (R Sun)       radius
         xx(4)  = flr                (L Sun)       luminosity
         xx(5)  = fmr                (M Sun)       mass (independent variable)
         xx(6)  = x                  hydrogen content        (by mass fraction)
         xx(7)  = z                  heavy element content   (by mass fraction)
         xx(8)  = zn                 nitrogen content for CNO burning (by mass)
                       if xx(8) < 0.0   then there is no nuclear burning - this is
                                   to be used for example in envelope integrations
    OUTPUT:
         yy(i) = d xx(i) / d xx(1)   i=1,2,3,4,5
    
    --------- NOTICE: subroutine rhs uses:
     subroutine state    to calculate thermodynamic functions
     subroutine opact    to calculate opacity (total: radiative and conductive)
     subroutine nburn    to calculate energy generation rate in hydrogen burning
    
    --------------- variables from common/const/ used by this subroutine:
         cgrrad = (sunl/sunm)/(16*pi*c*G)
         cdpm   = (G*sunm**2)/(4*pi*sunr**4)
         sunmr3 = (sunm)/(4*pi*sunr**3)
         sunml  = (sunm)/(sunl)
                   sunm, sunl, sunr are the solar mass, luminosity and radius
    
    --------------- auxillary variables:
         p      = total pressure      (c.g.s)
         pt     = d ln p / d ln t     at constant density
         pr     = d ln p / d ln rh    at constant temperature
         prad   = radiation pressure  (c.g.s.)
         grad   = d ln t / d ln p     at constant entropy (adiabatic temperature
                                                                     gradient)
         grrad  = d ln t / d ln p     radiative temperature gradient
         grt    = d ln t / d ln p     temperature gradient in the star
         grrh   = d ln rh / d ln p    logarithmic density gradient in the star
         dpm    = d ln p / d ( Mr / M Sun )    pressure gradient
    """
    
    t,rh,r,flr,fmr,x,z,zn = xx[:8]
    epsx = 0.0
    
    (p, pt, pr, pgas, prad, grad, qt, qr) = state(rh, t, x, z)
    fkap = opact(rh, t, x, z)
    if(zn > 0.0): 
        epsx, dxt = nburn(rh,t,x,zn)
    
    grrad = cgrrad*fkap*p/prad*flr/fmr
    grt = grad
    if(grt > grrad): 
        grt = grrad
    grrh = (1.0-grt*pt)/pr
    dpm = - cdpm*fmr/r**4
    
    yy = np.array([grt*dpm*t/p,
          grrh*dpm*rh/p,
          sunmr3/rh/r**2,
          sunml*epsx,
          1.0])
    return yy

def state(rh, t, x, z):
    """
    -- calculates equation of state and thermodynamic quantities
    INPUT:
          rh  = density (g/cm**2)
          t   = temperature (K)
          x   = hydrogen content (by mass fraction)
          z   = heavy element content (by mass fraction; helium content: y=1-x-z)
    OUTPUT:
          p    = total pressure (radiation+ions+electrons; dyn/cm**2)
          pt   = d log p / d log t     (at constant density)
          pr   = d log p / d log rh    (at constant temperature)
          prad = radiation pressure
          pgas = gas pressure (ions+electrons)
          grad = adiabatic temperature gradient: (d log t / d log p) at constant
                                                                        entropy
          qt   = t*cv, where "cv" is specific heat at constant volume (per gram)
          qr   = qt * ( d log t / d log rh) at constant entropy
     ---------------- NOTICE: full ionization is assumed, but electrons may be
                              partly degenerate, partly relativistic;
                              all formulae are written so as to avoid overflows;
    
     ------ auxiliary variables:
         y     = helium content (by mass)
         fmi   = mean molecular weight of ions (oxygen stands for heavy elements)
         fme   = mean molecular weight per free electron (full ionization)
         pion  = ion pressure (perfect gas law)
         pedr  = pressure of relativistically degenerate electrons
         pednr = pressure of non-relativistically degenerate electrons
         ped   = combined pressure of degenerate electrons
         pend  = pressure on non-degenerate electrons (perfect gas law)
         pe    = combined electron pressure (partly relativistic and degenerate)
         f     = a factor indicating how relativistic are the electrons
         pet   = d log pe / d log t   (at constant density)
         per   = d log pe / d log rh  (at constant temperature)
    """
    fkh, fk1, fk2, arad = (0.8251e8, 9.91e12, 1.231e15, 2.523e-15)
    
    y = 1.0-x-z
    fmi = 1.0/(x+y/4+z/16)
    fme = 2.0/(1.0+x)
    
    prad = arad*t**4
    pion = fkh/fmi*rh*t
    pedr = fk2*(rh/fme)**(4.0/3.0)
    pednr = fk1*(rh/fme)**(5.0/3.0)
    ped = pednr/np.sqrt(1+(pednr/pedr)**2)
    pend = fkh/fme*rh*t
    pe = pend*np.sqrt(1+(ped/pend)**2)
    
    pgas = pion+pe
    p = prad+pgas
    
    f = 5.0/3.0*(ped/pednr)**2+4.0/3.0*(ped/pedr)**2
    pet = (pend/pe)**2
    per = pet+(1-pet)*f
    
    pt = (4*prad+pion+pet*pe)/p
    pr = (pion+per*pe)/p
    qt = (12*prad+1.5*pion+pet/(f-1)*pe)/rh
    qr = p/rh*pt
    grad = 1/(qt/qr*pr+pt)
    
    return (p,pt,pr,pgas,prad,grad,qt,qr)

def nburn(rh,t,x,zn):
    """
    - calculates hydrogen burning rates for SCH and HEN programs, Oct. 22, 1984
    INPUT:
       rh  = density (g/cm**3)
       t   = temperature (K)
       x   = hydrogen content (by mass fraction)
       zn  = nitrogen content (by mass fraction)
    OUTPUT:
           epsx  = energy generation rate due to p-p and CNO hydrogen burning
                          (erg/g/sec)
           dxt   = hydrogen depletion rate (g/g/sec)  (negative number)
    
    ------ auxillary variables:
        epscno  = energy generation rate in CNO cycle
        epspp   = energy generation rate in the pp reaction
    """
    if(t < 1.0e6):
        """
        if temperature is below 10**6 K the burning rate is fixed at 1.0e-30
        to avoid underflows
        """
        epsx=1.0e-30
    else:
        t613 = (t/1.0e6)**(1.0/3.0)
        t623 = t613*t613
        epscno = x*zn*rh*np.exp(64.24-152.313/t613)/t623
        epspp = x*x*rh*np.exp(14.54-33.81/t613)/t623
        epsx = epscno+epspp
    dxt = -1.667e-19*epsx
    return (epsx,dxt)

def opact(rh,t,x,z):
    """
    -- calculates opacity (radiative and conductive)
    INPUT:
          rh  = density (g/cm**3)
          t   = temperature (K)
          x   = hydrogen content (by mass fraction)
          z   = heavy element content (by mass fraction; helium: y=1-x-z)
    OUTPUT:
          fkap = opacity (cm**2/g)
          
    ----------- auxiliary variables:
      y    = helium conetent (by mass)
      zav  = average charge of ions for conductive opacity
      fke  = electron scattering opacity
      fkk  = "Kramers" (i.e. bound-free, free-free and bound-bound) opacity
      fkhm = negative hydrogen ion opacity
      fkm  = molecular opacity
      fkr  = total radiative opacity
      fkc  = conductive opacity
      fkap = total opacity (radiative and conductive)
    
    -------- ad hoc correcting factors: A, B, C, introduced to beter reproduce
    -------- results from more sophisticated program (with ionizations included)
    """
    A = 6.0
    B = 0.001
    C = 0.001

    y = 1-x-z
    zav = x+4*y+8*z

    fke = 0.2*(1+x)/(1+2.7e11*rh/t/t)/(1+(t/4.5e8)**0.86)
    fkk = 12.0/A*(1+x)*(0.001+z)*rh*(1.0e7/t)**3.5
    if(t > 4.0e4):
        fkhm = 1.0e10
        fkm = 1.0e-5
        
    else:
        fkhm = B*65.0*np.sqrt(z*rh)*(t/3000.0)**7.7
        fkm = C*0.1*z
    fkr = fkm+1.0/(1.0/fkhm+1.0/(fke+fkk))
    fkc = 1.0e10
    if(rh > 0.00001):
        fkc = zav*2.6e-7*(t/rh)**2*(1+(rh/2.0e6)**(2.0/3.0))
    fkap = 1.0/(1.0/fkr+1.0/fkc)

    return fkap 
   
def solve(deriv):
    """
    ----- solves n=4 linear algebraic equations
     INPUT:
        deriv(n,n+1) = array of coefficients of "n" linear algebraic equations
     OUTPUT:
        delce(n)   = n corrections, without delmax(i) they should satisfy n
                     equations: sum over j: deriv(i,j)*delce(j)+deriv(i,n+1)=0
                     however, delce(i) are reduced so as to make them not larger
                     than delmax(i)
        facdel     = reduction factor for the corrections, if no reduction is
                     applied then facdel indicates the largest correction
    
     - auxillary variables:
        delmax(n)  = maximum acceptable values of corrections
        n = 4      = number of equations = number of unknowns (i.e. corrections)
    """
    #dimension deriv(4,5),delmax(4),delce(4)
    delmax = np.array([0.03, 0.1, 0.01, 0.1])
    n = 4
    nm = n-1
    np = n+1

      #do 1 k=1,nm
    for k in range(nm):
        kp=k+1
        fac1=deriv[k,k]
        for i in range(kp,n):#do 2 i=kp,n
            fac2=deriv[i,k]
            for j in range(kp,np): #do 3 j=kp,np
                deriv[i,j]=deriv[i,j]*fac1-deriv[k,j]*fac2
# ------ the matrix of derivatives, deriv(i,j) is now triangular

    delce[n] = -deriv[n,np]/deriv[n,n]
    for i in range(n): #do 4 i=2,n
        i1 = n-i+1
        i2 = i1+1
        delce[i1]=-deriv[i1,np]
        for j in range(i2,n):#do 5 j=i2,n
            delce[i1]=delce(i1)-deriv[i1,j]*delce[j]
        delce[i1]=delce[i1]/deriv[i1,i1]
#----- the unknowns delce(i) have been found

#----- check the size of "delce":
    dm = max(abs(delce/delmax))
    
    facdel = 1.0
    if(dm > 1.0): 
        facdel = dm

    delce[i] /= facdel

    facdel = dm

    return (delce,facdel)

def prout(k,xx,ic):
    """
    -- prints results of envelope or core integrations on the terminal and
    -- on a unit number "ip"; if ic<0 then prints only the headings
     INPUT:
       k     = integration step number, k<0 at the fitting point
       xx(1) = temperature)                    (K)
       xx(2) = density)                        (g/cm**3)
       xx(3) = r                               (R Sun)
       xx(4) = Lr                              (L Sun)
       xx(5) = Mr                              (M Sun)
       xx(6) = x = hydrogen content            (by mass fraction)
       xx(7) = z = heavy element content       (by mass fraction)
       xx(8) = zn = nitrogen content for CNO burning   (by mass fraction) 
       ip    = number of the output unit for "printing"
       ic    = control variable: ic>0 print results, ic<0 print table heading
    
     ------ variables from common/const/ used by subroutine prout.for:
       starm   = total stellar mass / solar mass
    """
    #global sunm,sunl,sunr,sunt,cgrrad,cdpm,sunmr3,sunml,g,c,fmod
    #dimension xx(10)

    if ic > 0:
        tlog, rhlog, rlog, fllog, fmlog= [np.log10(x) for x in xx[:5]]
        print ' %5i%11.5f%11.5f%11.5f%11.5f%11.5f' % (k,tlog,rhlog,rlog,fllog,fmlog)
    else:
        print '  step    log T     log rho   lg r/Rsun  lg L/Lsun  lg M/Msun'
        
        
if __name__ == '__main__':
    import numpy as np
    #global sunm,sunl,sunr,sunt,cgrrad,cdpm,sunmr3,sunml,g,c,fmod
    # -- set the values of mathematical and physical constants in the common/const/
    sunm = 1.99e33
    sunl = 3.83e33
    sunr = 6.96e10
    sunt = 5800.0
    
    #pi=3.1415927 use np.pi instead
    g = 6.67e-8
    c = 2.998e10
    fmod = np.log(10.0)
    
    cgrrad = (sunl/sunm)/(16.0*np.pi*c*g)
    cdpm = g/4.0/np.pi*(sunm/sunr**2)**2
    sunmr3 = sunm/(4.0*np.pi*sunr**3)
    sunml = sunm/sunl
    # -- end of constants -------------------------------------------------------
    
    print """ program "stellarstruct.py" calculates a series of hydrogen
 main sequence models
 with hydrogen content = X, and heavy element content = Z;"""
          
    
    #     Changed status='new' to status='unknown' in open statements --JJG
    #TODO: Filehandling
    #open(1,file='hms11.re1',status='unknown')
    #open(2,file='hms11.re2',status='unknown')
    
    # -- enter parameters for the stellar models from the keyboard --------------
    x,z = 0.7,0.3 #input('          enter chemical composition: X, Z = ')
      
    if x < 0.0:
        exit()
    print ' results of program "hms11.for", X, Z = %8.5f/%8.5f' % (x,z)
    print '   X        Z      log M  log Te   log L  log Tc  log rhc  log R  iter'
    
    fmlz = np.log10(1.2) #  input(' log (Mass of the first star / solar mass) = ')
    tel =  np.log10(6200.0) #  input('  log (effective temperature in degrees K) = ')
    fll =  np.log10(2.0) #  input('       log (luminosity / solar luminosity) = ')
    tcl =  np.log10(15.6e6)#  input('    log (central temperature in degrees K) = ')
    rhcl =  np.log10(150.0)  #  input('          log (central density in g/cm**3) = ')
    iprint = 1 # input('                  output control (integer) = ')
    
    starm=10.0**fmlz
    tel,fll,tcl,rhcl,rl,iter = sch(x, z, starm, tel, fll, tcl, rhcl, iprint)
    print ' %6.3f,%10.7f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%4i' % (x,z,fml,tel,fll,tcl,rhcl,rl,iter)



