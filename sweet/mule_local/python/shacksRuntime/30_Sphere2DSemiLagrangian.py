
from mule.JobCompileOptions import *

"""
IMPORTANT: Don't forget to call "mule.python.update_links" if updating anything
"""
class Sphere2DSemiLagrangian:

    def __init__(self):
        self.semi_lagrangian_approximate_sphere2d_geometry = None



    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''

        if self.semi_lagrangian_approximate_sphere2d_geometry != None:
            uniqueIDStr += '_spap'+str(self.semi_lagrangian_approximate_sphere2d_geometry)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        if self.semi_lagrangian_approximate_sphere2d_geometry != None:
            retRuntimeOptionsStr += ' --semi-lagrangian-approximate-sphere2d-geometry='+str(self.semi_lagrangian_approximate_sphere2d_geometry)

        return retRuntimeOptionsStr

    