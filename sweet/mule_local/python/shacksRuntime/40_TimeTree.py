
from mule.JobCompileOptions import *

"""
IMPORTANT: Don't forget to call "mule.python.update_links" if updating anything
"""
class TimeTree:

    def __init__(self):

        self.timestepping_method = None
        self.timestepping_order = None
        self.timestepping_order2 = None
        
        self.timestep_size = None
        self.max_timesteps_nr = None
    
        self.max_simulation_time = None
        self.max_wallclock_time = None


        return


    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        
        if not 'runtime.timestepping' in filter_list:
            if self.timestepping_method != None:
                if not 'runtime.timestepping_method' in filter_list:
                    tsm = self.timestepping_method
                    tsm = tsm.replace(",", "_")
                    tsm = tsm.replace("(", "_")
                    tsm = tsm.replace(")", "_")
                    tsm = tsm.replace("=", "_")
                    tsm = tsm.replace("order", "o") # compactify
                    tsm = tsm.replace("__", "_")
                    uniqueIDStr += '_tsm_'+tsm

                if not 'runtime.timestepping_order' in filter_list:
                    if self.timestepping_order != None:
                        uniqueIDStr += '_tso'+str(self.timestepping_order)

                if not 'runtime.timestepping_order2' in filter_list:
                    if self.timestepping_order2 != None:
                        uniqueIDStr += '_tsob'+str(self.timestepping_order2)

                if not 'runtime.semi_lagrangian' in filter_list:
                    if self.semi_lagrangian_max_iterations != None:
                        uniqueIDStr += '_sli'+str(self.semi_lagrangian_max_iterations)

                    if self.semi_lagrangian_convergence_threshold != None:
                        uniqueIDStr += '_slc'+str("{:0.5e}".format(self.semi_lagrangian_convergence_threshold))


            if not 'runtime.timestepping_size' in filter_list:
                if self.timestep_size != None:
                    # Leading number is the total number of digits!
                    if self.timestep_size < 1:
                        uniqueIDStr += '_dt'+str("{:0.8e}".format(self.timestep_size))
                    else:
                        uniqueIDStr += '_dt'+str("{:08.2f}".format(self.timestep_size))

            if not 'runtime.max_timesteps_nr' in filter_list:
                if self.max_timesteps_nr != None:
                    uniqueIDStr += '_T'+str(self.max_timesteps_nr).zfill(3)

            if not 'runtime.max_wallclock_time' in filter_list:
                if self.max_wallclock_time != None:
                    uniqueIDStr += '_W'+str(self.max_wallclock_time).zfill(6)


        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        if self.timestepping_method != None:
            retRuntimeOptionsStr += " --timestepping-method="+str(self.timestepping_method)

        if self.timestepping_order != None:
            retRuntimeOptionsStr += " --timestepping-order="+str(self.timestepping_order)

        if self.timestepping_order2 != None:
            retRuntimeOptionsStr += " --timestepping-order2="+str(self.timestepping_order2)

        if self.timestep_size != None:
            retRuntimeOptionsStr += " --dt="+str(self.timestep_size)

        if self.max_timesteps_nr != None:
            retRuntimeOptionsStr += " -T "+str(self.max_timesteps_nr)

        if self.max_simulation_time != None:
            retRuntimeOptionsStr += " -t "+str(self.max_simulation_time)

        if self.max_wallclock_time != None:
            retRuntimeOptionsStr += " --max-wallclock-time="+str(self.max_wallclock_time)
        
        return retRuntimeOptionsStr

    