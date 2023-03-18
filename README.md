# Etching-Sim
Cellular automata model for etching a substrate


Old python exploration when I was thinking about developing a mixed-model simulation to selectively ecth silicon.

Silicon has crystalline planes and tends to etch anisotropically along crystal plane directions.
For small optics or microfluidics, curved surfaces are useful.

In Stephen Senturia's book 'Microsystem Design' it becomes clear that the etch rate depends on the temperature and the pH.
The pH can be modified very locally with an electric field.
The temperature can be modified very locally with heat from flowing currents or from a structured light field.

This code was meant to explore how a designed MEMS structure could induce local heat gradients generating local etch gradients.
The result would be to model the ability to make smooth silicon microstructures through these applied gradients.

This code was developed as a 2D matrix approach to develop overlapping matrices that deal with one field variable at a time.
Field variables included Temperature, Material, and Voltage. 
Energy transport models would then be implemented to Etch, Heat, and relax an Electric Field Solver

The development stopped when I saw the validation using femtosecond laser 'super-etching' using dilute HF acid for silica based 3D micromachining
