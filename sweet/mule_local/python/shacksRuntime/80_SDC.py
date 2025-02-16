from mule.JobCompileOptions import *


"""
IMPORTANT: Don't forget to call "mule.python.update_links" if updating anything
"""
class SDC:
    def __init__(self):

        #
        # ExpIntegrators method:
        #
        # '': default
        # 'file': load coefficients from file (terry, cauchy & butcher coefficients)
        # 'direct': Use direct solution, if available
        # 'terry': deprecated T-REXI method
        # 'butcher': deprecated Butcher-REXI method
        #

        # SDC parameters
        self.sdc_params: SWEETFileDict = None
        self.sdc_file: str = None
        self.sdc_unique_id: str = None
        self.sdc_parallel: int = None


    def load_from_dict(self, d):

        if 'sdc_params' in d:
            self.sdc_params = d['sdc_params']

        if 'sdc_file' in d:
            self.sdc_file = d['sdc_file']

        if 'sdc_parallel' in d:
            self.sdc_parallel = d['sdc_parallel']


    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''

        if 'runtime.sdc' in filter_list:
            return uniqueIDStr

        if self.sdc_params != None:
            if self.sdc_unique_id == None:
                self.sdc_unique_id = self.sdc_params['idString']

        if not 'runtime.sdc.params' in filter_list:
            if self.sdc_params != None:
                uniqueIDStr += f"_{self.sdc_unique_id}"

        if not 'runtime.sdc.file' in filter_list:
            if self.sdc_file != None:
                uniqueIDStr += f"_{self.sdc_file}"

        if not 'runtime.sdc.parallel' in filter_list:
            if self.sdc_parallel != None:
                uniqueIDStr += "_p"+str(self.sdc_parallel)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''


        if self.sdc_file != None:
            retRuntimeOptionsStr += ' --sdc-file='+str(self.sdc_file)

        if self.sdc_parallel != None:
            retRuntimeOptionsStr += ' --sdc-parallel='+str(self.sdc_parallel)


        return retRuntimeOptionsStr
    
    
