# KROGER
Program for computing full and partial equilibria of point defects in solids


README

Welcome to Kroger in Matlab.  The purpose is extended calculations of defect concentrations, given thermochemistry and computed point defect formation energies.  

The main calculation engine is "defect_equilibrium_dark_with_frozen.m"

To run a scenario, we build up a "conditions" variable via a script (for now at least), then pass it and the desired defect database file to the calculation engine.  The script for now also handles formatting the output and saving it to disk.  

The example script is "defect_equilibrium_dark_with_frozen_script.m" (or something along those lines - look for "script" in the filename.  

The defect database is stored in a .mat file that loads "defects" into the workspace.  Defects is a character array so holds all the needed info about the defects and their chargestates - text names, charge, formation energy without Ef and mu terms, which sites the defect uses in the crystal, what elements need to be added or subtracted to make the defect... all that stuff.  

For better or worse, the current script keeps all sorts of things we have done in it as comments - that way we can turn things on and off by commenting/uncommenting.  


To run a scenario - best thing is to load the script into the matlab code editor, then put a debugging stop somewhere and press run.  Then step through each line to see what it is doing (including "step in" to go into the subroutines.  

The current database is "Ga2O3_Varley_all_defecs_new_072924.mat" - when you run the script it will ask you to pick the defect database - double click on this file.  

Then, it will open a dialog box asking you for the filename you want outputs to save into.  Just type a name without an extension - like "junk" and then look in the folder when it is all done.  

Some items are still in development, so ask lots of questions and/or get a demo session with Mike and his students 


Mike Scarpulla
University of Utah
mike.scarpulla@utah.edu
April 15, 2024
