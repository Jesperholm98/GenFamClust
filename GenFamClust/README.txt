INSTALLATION INSTRUCTIONS:

1) check and install a Python3.x version.
You can check if you currently have a 3.x version by typing:
'python3 --version', if you don't have a 3.x version installed you can install it by typing:
'sudo apt-get install python3.7'.

2) install numpy.
You can install the Numpy Python package by writing:
'sudo apt install python3-numpy'

3) install Python-dev for the c part to work.(#include<Python.h>)
you install this by typing:
'sudo apt-get install python3-dev'

4) try building/installing the nc script for your computer.
Go to the neighborhood_correlation-2.1 folder in GenFamClust where setup.py is located.
try running it by typing: 
'python3 setup.py build'
And after that:
'python3 setup.py install --user'

NOTE IF you are getting a "No module named 'distutils.core'" when trying to run the setup.py you need to type the following command:
'sudo apt-get install python3-distutils --reinstall'
and then
'python3 setup.py build'
And after that:
'python3 setup.py install --user' again.

Now a script called NC_standalone should be installed on your computer. If you cannot find the 'NC_standalone' script immediately, try restarting your device.


RUN INSTRUCITONS:

1) run NC module.
Go to the /Data folder inside GenFamClust and unpack the test data human_and_mouse by typing:
'gunzip human_mouse_bit_scores.dat.gz'
Now when you have some input data for NC, run NC by typing:
'NC_standalone -f human_mouse_bit_scores.dat -o nc.txt'
This should start the NC module with the c implementation and should take around 5-10 minutes for the included data with the c implementation.

2) Run GFC.
You can now run the rest of the program by first moving to the /GenFamClust folder where main.py is located and then typing:
'python3 main.py human_locs.txt mouse_locs.txt'.
The human_locs and mouse_locs are the synteny files for the respective genome. If you run the program with your own data you need to supply the relevant synteny files for both genomes.


If you want to run the program multiple times, be aware of that each module in GFC will write the result of itself to a specifically named file in the /Data folder.
To avoid loosing older results from modules you should either 1) Rename the result file of the module you may want to save or 2) Move it away from the /Data folder.
But keep in mind, the modules of GFC looks by default for input files from previous modules in the /Data folder, and writes to the result files in there as well.

If you want to run individual modules in GFC you can go to main.py file and comment away the modules you don't want to be run during the execution.
