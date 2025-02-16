
from mule.JobCompileOptions import *

class XBraid:

    def __init__(self):
        # XBraid
        self.xbraid = 'none'
        self.xbraid_scalar = 'disable'
        self.xbraid_cart2d = 'disable'
        self.xbraid_sphere2d = 'disable'
        self.xbraid_cart2d_swe = 'disable'
        self.xbraid_cart2d_burgers = 'disable'


    def getSConsParams(self):
        retval = ''

        # XBraid
        retval += ' --xbraid='+self.xbraid
        retval += ' --xbraid-scalar='+self.xbraid_scalar
        retval += ' --xbraid-cart2d='+self.xbraid_cart2d
        retval += ' --xbraid-sphere2d='+self.xbraid_sphere2d
        retval += ' --xbraid-cart2d-swe='+self.xbraid_cart2d_swe
        retval += ' --xbraid-cart2d-burgers='+self.xbraid_cart2d_burgers


        return retval


    def sconsAddOptions(self, scons):

        scons.AddOption(    '--xbraid',
                dest='xbraid',
                type='choice',
                choices=['none', 'mpi'],
                default='none',
                help='Enable XBBraid (none,  mpi) [default: %default]\nOnly works, if XBraid is supported by the simulation'
        )
        self.xbraid = scons.GetOption('xbraid')

        scons.AddOption(    '--xbraid-scalar',
                dest='xbraid_scalar',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for scalar problems (enable, disable) [default: %default]'
        )
        self.xbraid_scalar = scons.GetOption('xbraid_scalar')

        scons.AddOption(    '--xbraid-cart2d',
                dest='xbraid_cart2d',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the cart2d (enable, disable) [default: %default]'
        )
        self.xbraid_cart2d = scons.GetOption('xbraid_cart2d')

        scons.AddOption(    '--xbraid-sphere2d',
                dest='xbraid_sphere2d',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the sphere2d (enable, disable) [default: %default]'
        )
        self.xbraid_sphere2d = scons.GetOption('xbraid_sphere2d')

        scons.AddOption(    '--xbraid-cart2d-swe',
                dest='xbraid_cart2d_swe',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for SWE on the cart2d (enable, disable) [default: %default]'
        )
        self.xbraid_cart2d_swe = scons.GetOption('xbraid_cart2d_swe')

        scons.AddOption(    '--xbraid-cart2d-burgers',
                dest='xbraid_cart2d_burgers',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for Burgers on the cart2d (enable, disable) [default: %default]'
        )
        self.xbraid_cart2d_burgers = scons.GetOption('xbraid_cart2d_burgers')


    def sconsValidateOptions(self):
        
        if not self.xbraid == 'none':
            if self.program == "parareal_ode":
                self.xbraid_scalar = 'enable';
            elif self.program == 'swe_cart2d' or self.program == 'burgers':
                self.xbraid_cart2d = 'enable';
                if self.program == 'swe_cart2d':
                    self.xbraid_cart2d_swe = 'enable';
                elif self.program == 'burgers':
                    self.xbraid_cart2d_burgers = 'enable';
            elif self.program == 'swe_sphere2d':
                self.xbraid_sphere2d = 'enable';
                
                
    def sconsAddFlags(self, env):
        if self.xbraid == 'none':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID=0'])
        elif self.xbraid == 'mpi':
            env.Append(CXXFLAGS=['-Ilocal_software/local/include/xbraid'])
            env.Append(LIBS=['braid'])
            env.Append(CXXFLAGS=['-DSWEET_XBRAID=1'])
        else:
            raise Exception("Invalid option '"+str(self.xbraid)+"' for XBraid")
        
        
        if self.xbraid_scalar == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_SCALAR=1'])
        if self.xbraid_cart2d == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_CART2D=1'])
        if self.xbraid_sphere2d == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_SPHERE2D=1'])
        if self.xbraid_cart2d_swe == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_CART2D_SWE=1'])
        if self.xbraid_cart2d_burgers == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_CART2D_BURGERS=1'])
        

