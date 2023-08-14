# C-Fluid-Sim

Based on [Real-Time Fluid Dynamics for Games](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf) and [Fluid Simulation for dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html).

Part of GUI grid code from [sdl-grid](https://github.com/catsocks/sdl-grid/tree/master)

## SDL
Manually added SDL to project by adding `SDL2-2.28.2\include` to the `C/C++->Additional Include Directories`, and adding `SDL2-2.28.2\lib\x64\SDL2.lib` and `SDL2-2.28.2\lib\x64\SDL2main.lib` to `Linker->Input->Additional Dependencies`.  

`SDL.dll` has also been added to project and been set to be coppied to output during build step.
