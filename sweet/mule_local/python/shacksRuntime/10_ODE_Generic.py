from mule.JobCompileOptions import *


"""
IMPORTANT: Don't forget to call "mule.python.update_links" if updating anything
"""
class ODE_Generic:

    def __init__(self):
        self.ode = None
        
        self.ode_dahlquist_lambda1 = None
        self.ode_dahlquist_lambda2 = None
        self.ode_dahlquist_lambda3 = None

        self.ode_dahlquist_mu = None
        self.ode_dahlquist_phi = None

    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        if not 'runtime.ode_generic' in filter_list:
            if self.ode != None:
                uniqueIDStr += '_ode_'+str(self.ode)
                
            if self.ode_dahlquist_lambda1 != None:
                uniqueIDStr += '_lam1_'+str(self.ode_dahlquist_lambda1).replace(",", "_")
                
            if self.ode_dahlquist_lambda2 != None:
                uniqueIDStr += '_lam2_'+str(self.ode_dahlquist_lambda2).replace(",", "_")
                
            if self.ode_dahlquist_lambda3 != None:
                uniqueIDStr += '_lam3_'+str(self.ode_dahlquist_lambda3).replace(",", "_")
                
            if self.ode_dahlquist_mu != None:
                uniqueIDStr += '_mu_'+str(self.ode_dahlquist_mu).replace(",", "_")
                
            if self.ode_dahlquist_phi != None:
                uniqueIDStr += '_phi_'+str(self.ode_dahlquist_phi).replace(",", "_")
                
        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        def cplxToStr(x):
            return str(x.real)+","+str(x.imag)

        if self.ode != None:
            retRuntimeOptionsStr += " --ode="+str(self.ode)

        if self.ode_dahlquist_lambda1 != None:
            retRuntimeOptionsStr += " --ode-dahlquist-lambda1="+cplxToStr(self.ode_dahlquist_lambda1)
            
        if self.ode_dahlquist_lambda2 != None:
            retRuntimeOptionsStr += " --ode-dahlquist-lambda2="+cplxToStr(self.ode_dahlquist_lambda2)
            
        if self.ode_dahlquist_lambda3 != None:
            retRuntimeOptionsStr += " --ode-dahlquist-lambda3="+cplxToStr(self.ode_dahlquist_lambda3)
            
        if self.ode_dahlquist_mu != None:
            retRuntimeOptionsStr += " --ode-dahlquist-mu="+cplxToStr(self.ode_dahlquist_mu)
            
        if self.ode_dahlquist_phi != None:
            retRuntimeOptionsStr += " --ode-dahlquist-phi="+cplxToStr(self.ode_dahlquist_phi)

        return retRuntimeOptionsStr
