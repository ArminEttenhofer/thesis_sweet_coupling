
Go to the main SWEET folder and compile your code with

```
scons --program=programs/PDE_SWESphere --mode=release --gui=enable
```

Start your simulation, e.g., with 
```
./build/programs/PDE_SWESphere_COMP_plspec_pldeal_spspec_spdeal_gui_fft_thomp_release --benchmark-name=galewsky -M 256 --dt=150 --timestepping-method=l_irk_n_erk --timestepping-order=2 -t $((24*60*60*10))
```

Activate screenshot generation and run the simulation. This will create a bunch of .bmp output files.

We now want to create a video:

* Use the scripts to convert the bmp files to png.

