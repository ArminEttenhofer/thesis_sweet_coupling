from mule.JobCompileOptions import *

"""
IMPORTANT: Don't forget to call "mule.python.update_links" if updating anything
"""
class Parallelization:

    def __init__(self):
        self.sh_setup_num_threads = None


    def load_from_dict(self, d):
        if 'sh_setup_num_threads' in d:
            self.sh_setup_num_threads = int(d['sh_setup_num_threads'])


    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''

        if not 'parallelization' in filter_list:
            if self.sh_setup_num_threads != None:
                uniqueIDStr += '_nts'+str(self.sh_setup_num_threads)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        if self.sh_setup_num_threads != None:
            retRuntimeOptionsStr += ' --sh-setup-num-threads='+str(self.sh_setup_num_threads)

        return retRuntimeOptionsStr

    
