#! /usr/bin/env python3

import sys
import os
import copy

from mule.JobCompileOptions import *
from mule.InfoError import *
from mule.SWEETFileDict import SWEETFileDict
import mule.Shacks

__all__ = ['JobRuntimeOptions']


class JobRuntimeOptions(InfoError):

    def __init__(self, dummy_init = False):

        self.init_phase = True
        
        InfoError.__init__(self, "JobRuntimeOptions")
        
        #self.shacksRuntime = mule.Shacks.getShacksDict("shacksRuntime").values()
        self.shacksRuntime = mule.Shacks.Shacks("shacksRuntime")
        
        for _ in self.shacksRuntime.values():
            _.__init__(self)

        # String to job directory.
        self.p_job_dirpath = None


        self.space_res_spectral = None
        self.space_res_grid = None
        self.sh_setup_num_threads = None
        self.sh_setup_verbosity = None


        self.verbosity = 0

        self.gui = None

        # Use 14 digits per default
        self.floating_point_output_digits = 12

        self.sphere2d_radius = None

        self.cart2d_domain_size = None    # Cart2D: Domain size

        self.space_grid_use_c_staggering = None
        self.space_use_spectral_basis_diffs = None
        
        self.reuse_plans = None

        self.f_sphere2d = None
        self.semi_lagrangian_max_iterations = None
        self.semi_lagrangian_convergence_threshold = None


        self.instability_checks = None
        
        self.normal_mode_analysis = None


        # SDC parameters
        self.paramsSDC: SWEETFileDict = None
        self.sdcParallel: int = None


        self.gravitation= None
        self.h0 = None
        self.sphere2d_rotating_coriolis_omega = None
        self.viscosity = None
        self.viscosity_order = None

        # Specify benchmark name
        self.benchmark_name = None
        self.benchmark_galewsky_umax = -1
        self.benchmark_galewsky_hamp = -1
        self.benchmark_galewsky_phi2 = -1
        self.benchmark_galewsky_geostropic_setup = None
        self.benchmark_normal_modes_case = None

        self.benchmark_advection_rotation_angle = None
        self.benchmark_advection_velocity = None

        #self.uselineardiv = None
        self.use_nonlinear_only_visc = None

        self.compute_errors = None

        self.comma_separated_tags = None

        #
        # User defined parameters
        # Each new entry must set three values:
        # 
        # 'id': The id for the unique id
        # 'value': The value of the parameter
        # 'option': The program option
        #
        # example:
        # self.user_defined_parameters['test-mode'] = {'id': 'tm', 'value': 1, 'option': '--test-mode='}
        #
        self.user_defined_parameters = {}

        self.init_phase = False



    def __setattr__(self, name, value):

        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value



    def load_from_dict(self, d):
        
        for _ in self.shacksRuntime.values():
            _.load_from_dict(self, d)
        
        if 'timestep_size' in d:
            self.timestep_size = float(d['timestep_size'])

        if 'reuse_plans' in d:
            if self.reuse_plans is not None:
                self.reuse_plans = int(d['reuse_plans'])

        if 'comma_separated_tags' in d:
            self.comma_separated_tags = d['comma_separated_tags']



    def getUniqueID(self, compileOptions : JobCompileOptions, filter_list : list = []):
        uniqueIDStr = ''

        if not 'runtime.benchmark' in filter_list:
            if self.benchmark_name != None:
                uniqueIDStr += '_b'+str(self.benchmark_name)

            if not 'runtime.benchmark.galewsky_params' in filter_list:
                if self.benchmark_galewsky_umax != None:
                    uniqueIDStr += '_bgu'+str("{:.4E}".format(self.benchmark_galewsky_umax))

                if self.benchmark_galewsky_hamp != None:
                    uniqueIDStr += '_bgh'+str("{:.4E}".format(self.benchmark_galewsky_hamp))

                if self.benchmark_galewsky_phi2 != None:
                    uniqueIDStr += '_bgp'+str("{:.4E}".format(self.benchmark_galewsky_phi2))

                if self.benchmark_galewsky_geostropic_setup != None:
                    uniqueIDStr += '_bgs'+self.benchmark_galewsky_geostropic_setup

        if not 'runtime.normal_modes_params' in filter_list:
            if self.benchmark_normal_modes_case != None:
                uniqueIDStr += '_bcase'+str(self.benchmark_normal_modes_case)

        if not 'runtime.simparams' in filter_list:
            if self.gravitation!= None:
                uniqueIDStr += '_g'+str("{:05.2f}".format(self.gravitation))
            if self.h0 != None:
                uniqueIDStr += '_h'+str("{:010.3f}".format(self.h0))
            if self.sphere2d_rotating_coriolis_omega != None:
                uniqueIDStr += '_f'+str("{:e}".format(self.sphere2d_rotating_coriolis_omega))

            if compileOptions.sphere2d_spectral_space == 'enable':
                if self.sphere2d_radius != None:
                    uniqueIDStr += '_a'+str(self.sphere2d_radius)
                if self.f_sphere2d != None:
                    uniqueIDStr += '_fsph'+str(self.f_sphere2d)

            if self.viscosity != None:
                uniqueIDStr += '_u'+str(self.viscosity)
            if self.viscosity_order != None:
                uniqueIDStr += '_U'+str(self.viscosity_order)

            if self.benchmark_advection_rotation_angle != None:
                uniqueIDStr += '_ar'+str(self.benchmark_advection_rotation_angle)

            if self.benchmark_advection_velocity != None:
                uniqueIDStr += '_av'+str(self.benchmark_advection_velocity).replace(",", "_")

        if 'timestep' in filter_list:
            raise Exception("Deprecated")
        
        if not 'runtime.paramSDC' in filter_list:
            if self.paramsSDC != None:
                uniqueIDStr += f"_{self.paramsSDC['idString']}"


        for _ in self.shacksRuntime.values():
            uniqueIDStr += _.getUniqueID(self, compileOptions, filter_list)

        
        if not 'runtime.disc_space' in filter_list:
            if self.space_res_spectral != None:
                if isinstance(self.space_res_spectral, (list, tuple)):
                    uniqueIDStr += '_M'+str("x".join([str(x).zfill(4) for x in self.space_res_spectral]))
                else:
                    uniqueIDStr += '_M'+str(self.space_res_spectral).zfill(4)

            if self.space_res_grid != None:
                if isinstance(self.space_res_grid, (list, tuple)):
                    uniqueIDStr += '_N'+str("x".join([str(x).zfill(4) for x in self.space_res_grid]))
                else:
                    uniqueIDStr += '_N'+str(self.space_res_grid).zfill(4)

            if self.cart2d_domain_size != None:
                if isinstance(self.cart2d_domain_size, (list, tuple)):
                    uniqueIDStr += '_X'+str("x".join([str(x).zfill(4) for x in self.cart2d_domain_size]))
                else:
                    uniqueIDStr += '_X'+str(self.cart2d_domain_size)

            if self.space_use_spectral_basis_diffs != None:
                uniqueIDStr += '_spd'+str(self.space_use_spectral_basis_diffs)

        if not 'runtime.reuse_plans' in filter_list:
            if self.reuse_plans != None:
                uniqueIDStr += '_plans'+str(self.reuse_plans)

        if not 'runtime.comma_separated_tags' in filter_list:
            if self.comma_separated_tags != None:
                uniqueIDStr += '_tags'+str(self.comma_separated_tags)

        if uniqueIDStr != '':
            uniqueIDStr = "RT"+uniqueIDStr

        for key, param in self.user_defined_parameters.items():
            if param['id'] != '':
                uniqueIDStr += '_'+param['id']+str(param['value'])

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        self.cleanup_options()

        def addOption(name, value, option):

            if isinstance(option, list):
                # One attribute for multiple options
                exec = ''
                if hasattr(value, '__iter__'):
                    # Multiple values for multiple options
                    for opt, val in zip(option, value):
                        exec += f' {opt}{val}'
                else:
                    # One same value for multiple options
                    for opt in option:
                        exec += f' {opt}{value}'
                return exec

            if isinstance(value, SWEETFileDict):
                # Attribute is a SWEETFileDict
                filePath = os.path.join(self.p_job_dirpath, f'{name}.sweet')
                value.writeToFile(filePath)
                exec = f' {option}{filePath}'
                return exec

            # Default case
            if isinstance(value, (list, tuple)):
                value = ",".join([str(x) for x in value])
            exec = f' {option}{value}'

            return exec        

        optionList = {
            # Program options
            'verbosity': '-v ',
            'gui': '-G ',
            'floating_point_output_digits': '-d ',
            'comma_separated_tags': '--comma-separated-tags=',
            
            # Space discretization
            'space_res_spectral': '-M ',
            'space_res_grid': '-N ',
            'space_grid_use_c_staggering': '--space-grid-use-c-staggering=',
            'space_use_spectral_basis_diffs': '-S ',
            'cart2d_domain_size': ['-X ', '-Y '],

            # Benchmark options
            'benchmark_name': '--benchmark-name=',
            'benchmark_normal_modes_case': '--benchmark-normal-modes-case=',
            'benchmark_galewsky_geostropic_setup': '--benchmark-galewsky-geostropic-setup=',
            'benchmark_advection_rotation_angle': '--benchmark-advection-rotation-angle=',
            'benchmark_advection_velocity': '--benchmark-advection-velocity=',
            
            # Time-stepping options
            'instability_checks': '--instability-checks=',

            # PDE options
            'viscosity': '--pde-viscosity=',
            'viscosity_order': '--pde-viscosity-order=',
            'use_nonlinear_only_visc': '--use-nonlinear-only-visc=',
            'gravitation': '--pde-gravitation=',
            'h0': '--pde-h0=',
            'sphere2d_rotating_coriolis_omega': '-f ',
            'sphere2d_radius': '-a ',
            'f_sphere2d': '-F ',

            # Time discretization
            'semi_lagrangian_max_iterations': '--semi-lagrangian-max-iterations=',
            'semi_lagrangian_convergence_threshold': '--semi-lagrangian-convergence-threshold=',
            'paramsSDC': '--sdc-file=',
            'sdcParallel': '--sdc-parallel=',
            
            # Analysis
            'normal_mode_analysis': '--normal-mode-analysis-generation=',
            'compute_errors': '--compute-errors=',
            'reuse_plans': '--reuse-plans=',
        }

        for name, option in optionList.items():
            value = self.__getattribute__(name)
            if value is not None:
                retRuntimeOptionsStr += addOption(name, value, option)

        for _ in self.shacksRuntime.values():
            retRuntimeOptionsStr += _.getRuntimeOptions(self)

        for param in self.user_defined_parameters.values():
            retRuntimeOptionsStr += ' '+param['option']+str(param['value'])

        return retRuntimeOptionsStr



    def cleanup_options(self):

        # Cleanup old format of reuse plans to new one
        if self.reuse_plans == -1:
            self.reuse_plans = "quick"
        elif self.reuse_plans == 0:
            self.reuse_plans = "save"
        elif self.reuse_plans == 1:
            self.reuse_plans = "load"
        elif self.reuse_plans == 2:
            self.reuse_plans = "require_load"



    def get_jobscript_plan_exec_prefix(self, compile, runtime):
        content = ""

        runtime.cleanup_options()

        plan_files = []
        if compile.cart2d_spectral_space == 'enable':
            plan_files.append('sweet_fftw')

        if compile.sphere2d_spectral_space == 'enable':
            plan_files.append('shtns_cfg')
            plan_files.append('shtns_cfg_fftw')


        if runtime.reuse_plans == "quick" or runtime.reuse_plans is None:
            # Quick plan generation mode, nothing to do
            pass

        elif "load" in runtime.reuse_plans:
            content += "\n"
            # Reuse plans if available
            # No error if plans don't exist
            for i in plan_files:
                if runtime.reuse_plans == "require_load":
                    content += "cp ../"+i+" ./ || exit 1\n"
                else:
                    content += "cp ../"+i+" ./ 2>/dev/null\n"

        elif runtime.reuse_plans == "save":
            # Create plans, don't load/store them
            pass

        else:
            raise Exception("Invalid reuse_plans value '"+str(runtime.reuse_plans)+"'")

        return content



    def get_jobscript_plan_exec_suffix(self, compile, runtime):
        content = ""

        runtime.cleanup_options()

        if runtime.reuse_plans is not None:
            plan_files = []
            if compile.cart2d_spectral_space == 'enable':
                plan_files.append('sweet_fftw')
    
            if compile.sphere2d_spectral_space == 'enable':
                plan_files.append('shtns_cfg')
                plan_files.append('shtns_cfg_fftw')
    
            #
            # Reusing plans assumes them to be stored in the folder one level up in the hierarchy
            #
            if runtime.reuse_plans == "quick":
                # Quick plan generation mode, nothing to do
                pass
    
            elif "save" in runtime.reuse_plans:
                # Write back plans to main directory to reuse them
                content += "\n"
                for i in plan_files:
                    content += "cp ./"+i+" ../ 2>/dev/null\n"
                content += "\n"
                pass
    
            elif "load" in runtime.reuse_plans:
                pass
    
            else:
                raise Exception("Invalid reuse_plans value '"+str(runtime.reuse_plans)+"'")


        return content


if __name__ == "__main__":

    p = JobRuntimeOptions()
    i = p.getRuntimeOptions()
    p.info(i)

    p.info("FIN")

