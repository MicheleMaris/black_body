__DESCRIPTION__=""" 
CGS       = an object returning the CGS units and constants
BlackBody = an object returning black body fluxes and conversion facilities
CMB       = a specialized blackbody for CMB

See: 
   import blackbody 
   help(blackbody.CGS)
   help(blackbody.BlackBody)
   help(blackbody.CMB)

The package manages numpy.array().

Version 1.0 - 2011 Set 1 - 2012 Apr 5 - M.Maris
"""

class __physical_parameters_cgs(object) :
   def __init__(self) :
      self.K = 1.380658e-16    # erg/K
      self.c = 2.99792458e10   # cm/sec
      self.h = 6.6260755e-27   # erg / sec i.e. erg Hz
      self.flux_Jansky = 1e23  # (cm2 sec cm sterad)/erg * Jansky
      self.ScaleToMJ = 1.e17   # 1.e23 Jy/erg cm2 sec * 1e-6 MJy/Jy
   def keys(self) :
      return self.__dict__.keys()
   def __call__(self,name) :
      units = {}
      units['K'] = 'erg/K'
      units['c'] = 'cm/sec'
      units['h'] = 'erg/sec'
      units['flux_Jansky'] = '(cm2 sec cm sterad)/erg * Jansky'
      units['ScaleToMJ'] = '1.e23 Jy/erg cm2 sec * 1e-6 MJy/Jy'
      units['TCMB'] = 'K'
      return "%s = %e : %s "%(name,self.__dict__[name],units[name])
CGS=__physical_parameters_cgs()

class __black_body(object) :
   def __init__(self) :
      return
   def __call__(self,FreqGHz,T) :
      import numpy as np
      
      return self.bbn_cgs(FreqGHz*1e9,T)*1e23*self.valid(FreqGHz,T)
   def valid(self,FreqGHz,T) :
      return (FreqGHz >= 0.) * (T >= 0.)
   def bbl_cgs(self,Lambda,T) :
      """
!
! BB(lambda,T) thermal radiance function cgs
!
! lambda in cm
!
! bbl_cgs = erg/(cm2 sec cm sterad)
!
"""
      import numpy as np
      FT = Lambda*Lambda*Lambda*Lambda*Lambda
      FT = 2.*CGS.h*CGS.c*CGS.c/FT;
      ET = np.exp( CGS.h*CGS.c/(Lambda*CGS.K*T) );
      ET = ET - 1.
      return FT / ET
   def bbn_cgs(self,nu,T) :
      """
!
! BB(nu,T) thermal radiance function cgs
!
! nu in Hz
!
! bbn_cgs = erg/(cm2 sec cm sterad)
!
      """
      import numpy as np
      FT = nu*nu*nu
      ET = CGS.c*CGS.c
      FT = 2.*CGS.h*FT/ET
      ET = CGS.h*nu/(CGS.K*T)
      ET = np.exp(ET)
      ET = ET - 1
      return FT / ET
   def bbl_rj_cgs(self,Lambda,T) :
      """
!
! BB thermal radiance function cgs in Rayleight-Jeans approx
!
! lambda in cm
! T in K
! rj_cgs = erg/(cm2 sec cm sterad)
!
      """
      FT = Lambda*Lambda*Lambda*Lambda
      return 2.*CGS.K*T*CGS.c/FT 
   def bbn_rj_cgs(self,nu,T) :
      """
!
! BB thermal radiance function cgs in Rayleight-Jeans approx
!
! nu in Hz
! T in K
!
! rj_cgs = erg/(cm2 sec Hz sterad)
!
      """
      FT = (nu/CGS.c)
      FT = FT*FT
      return 2.* CGS.K*T*FT
   def Krj2MJysr(self,nu_ghz,Hz=False) :
      """generates the conversion factor from K_rj to MJy/sr for a given frequency in GHz (Hz=True to pass it in Hz)"""
      if Hz :
         return self.bbn_rj_cgs(nu_ghz,1.)*CGS.ScaleToMJ
      return self.bbn_rj_cgs(nu_ghz*1e9,1.)*CGS.ScaleToMJ
   def Tb(self,nu,Bn,FreqGHz=False,MJySr=False) :
      """returns for a given CGS brightness the correspondiing temperature"""
      # B=2*h*n^3/c^2 /(exp(h*nu/K/T)-1)
      # 1/(log( 1/(B/(2*h*nu^3/c^2)+1 )/(h*nu)*K)
      import numpy as np
      if FreqGHz : 
         f=nu*1e9
      else :
         f=nu
      if MJySr :
         B=Bn/CGS.ScaleToMJ
      else :
         B=Bn
      a=CGS.h*f/(CGS.K*np.log((2*CGS.h*f**3/CGS.c**2)/B+1))
      return a
   def bbn_diff(self,nu_ghz,T,Hz=False,MJySr=True) :
      """

%
% [bbn_diff,bbn_diff_ratio]=bbn_diff_mks(nu_hz,T)
%
% bbn_diff =  derivative of the BB thermal radiance in mks
%             for a temperature change
%
% bbn_diff_ratio = bbn_diff/bbn
%
% nu_hz in Hz
% T in K
%
% bbn_diff_mks = erg/(cm2 sec Hz sterad)/K
% bbn_diff_ratio = 1/K
%
% 1 Jy/sterad = 1e-23 erg/(cm2 sec Hz sterad)
%
% the derivative is defines as:
%
%     d(Bnu(T))/dT = Bnu(T) X exp(X)/(exp(X)-1) / T
%
% the variation is defined as
%
%     DeltaBnu(T) = Bnu(T) X exp(X)/(exp(X)-1) DeltaT/ T
%
      """
      import numpy as np
      if not Hz :
         _nu=nu_ghz*1e9
      else :
         _nu=nu_ghz
      if MJySr :
         cf = CGS.ScaleToMJ
      else :
         cf = 1.
      FT = 2.*CGS.h*(_nu**3)/CGS.c**2;
      X =  CGS.h*_nu/(CGS.K*T);
      ET = np.exp(X);
      Lbbn = FT / (ET - 1);
      FACT = X*ET/(ET-1)/T;
      return {'bbn_diff':cf*Lbbn*FACT,'bbn_diff_ratio':cf*FACT,'bbn':cf*Lbbn}
BlackBody=__black_body()

class __cmb(object) :
   def __init__(self) :
      import numpy as np
      self.Tcmb=2.72548 #Fixsen, D. J. 1999, Volume 707, Issue 2, pp. 916-920 (2009) #2.725 # K 
      self.DeltaTDipole = 3.3e-3; #K
      self.sqC2 = np.sqrt(211)*1e-6/np.sqrt(2*(2+1)/(2*np.pi)); #K sqrt of Cl derived from quadrupole moment = l*(l+1)/(2pi)*Cl
      self.sqC3 = np.sqrt(1041)*1e-6/np.sqrt(3*(3+1)/(2*np.pi)); #K
      self.sqC100 = np.sqrt(3032.9299)*1e-6/np.sqrt(100*(100+1)/(2*np.pi)); #K
      self.source = "Source: WMAP"
   def __call__(self,FreqGHz,MJySr=True) :
      if MJySr :
         return BlackBody.bbn_cgs(FreqGHz*1e9,self.Tcmb)*1e23
      return BlackBody.bbn_cgs(FreqGHz*1e9,self.Tcmb)
      #return "Source: WMAP"
   def fluctuations(self,FreqGHz) :
      """Return a dictionary with the list of most important cmb components"""
      D=BlackBody.bbn_diff(FreqGHz*1e9,self.Tcmb,Hz=True);
      cmbf={}
      cmbf['FreqGHz'] = FreqGHz;
      cmbf['Inu'] = D['bbn']/1e-23/1e6;
      cmbf['dInu_dT'] = D['bbn_diff']/1e-23/1e6;
      cmbf['dlogInu_dT'] = D['bbn_diff_ratio'];
      cmbf['Dipole'] = cmbf['dInu_dT']*self.DeltaTDipole;
      cmbf['sqC2'] = cmbf['dInu_dT']*self.sqC2;
      cmbf['sqC3'] = cmbf['dInu_dT']*self.sqC3;
      cmbf['sqC100'] = cmbf['dInu_dT']*self.sqC100;
      return cmbf
   def Kcmb2MJysr(self,FreqGHz,TKCMB) :
      """Converts Kcmb to MJy/sr"""
      c=BlackBody.bbn_diff(FreqGHz,self.Tcmb)['bbn_diff']
      return c*TKCMB
   def MJysr2Kcmb(self,FreqGHz,IMJY_Sr) :
      """Converts MJy/sr to Kcmb"""
      c=BlackBody.bbn_diff(FreqGHz,self.Tcmb)['bbn_diff']
      return IMJY_Sr/c
   def etaDeltaT(self,FreqGHz) :
      """Computes EtaDeltaT(nu) (see Eq.(34) of Zacchei et al. (2011), A&A, 536, A5) 
      defined as (\partial B(\nu,T) / \partial T)_{T=Tcmb} / (2 K \nu^2/c^2)
      """
      import numpy as np
      X=CGS.h*FreqGHz*1e9/(CGS.K*self.Tcmb)
      eX=np.exp(X)
      return eX*(X/(eX-1.))**2
CMB=__cmb()
