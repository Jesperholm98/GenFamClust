Installation instructions:

1) check and install a Python3.x version.
You can check if you currently have a 3.x version by typing:
'python3 --version', if you don't have a 3.x version installed you can install it by typing:
'sudo apt-get install python3.7'.

2) install numpy.
You can install the Nunmpy Python package by writing:
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

Now a script called NC_standalone should be installed on your computer(jag behövde köra build och sedan install med setup IGEN för att skriptet skulle installeras korrekt, alternativt starta datorn och sedan kör steg 4) IGEN.)
(Första gången jag gick igenom steg 1-4 så hittades inte scriptet, men när jag testade steg 4 dagen efter så installerades scriptet och kunde hittas.)

5) run NC module(skriptet).
Go to the /Data folder inside GenFamClust and unpack the test data human_and_mouse by typing:
'gunzip human_mouse_bit_scores.dat.gz'
Now you should be able to run NC module by typing:
'NC_standalone -f human_mouse_bit_scores.dat -o nc.txt'
This should start the NC module with the c implementation and should take around 5-10 mins.

6) Run GFC.
You can now run the rest of the program by first moving to the /GenFamClust folder where main.py is located and then running:
'python3 main.py 1 2'. (1 2 placeholder)
(this will at this point just run the SyS module for testing purposes)
