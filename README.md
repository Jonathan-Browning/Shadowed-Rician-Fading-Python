# Shadowed-Rice-Fading-Matlab

The shadowed Rician fading model implemented in python. Plots the theoretical and simulated, envelope and phase porbability density functions (PDFs)

For more information this model please refer to Browning's paper (me):
"The Rician Complex Envelope under Line of Sight Shadowing".

This project uses PySimpleGUI, numpy, matplotlib and tkinter.

This project was developed on a windows OS, using Spyder IDE with Python 3.8. All the dependencies where installed by anaconda.

The input K accepts values in the range 0.001 to 50.
The input m accepts values in the range 0.001 to 50.
The input \hat{r} accepts values in the range 0.5 to 2.5.
The input \phi accepts values in the range -pi to pi.

Runing main.py to start the GUI displays:
  
![ScreenShot](https://raw.github.com/Jonathan-Browning/Shadowed-Rician-Fading-Python/main/docs/window.png)

Entering values for the Rician K factor, m the shadowing severity, the root mean sqaure of the signal \hat{r}, and \phi the phase parameter:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Shadowed-Rician-Fading-Python/main/docs/inputs.png)

The theoretical evenlope PDF and phase PDFs are plotted to compare with the simulation and gives the execution time for the theoretical calculations and simulations together:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Shadowed-Rician-Fading-Python/main/docs/results.png)