This folder contains some sample data for microrheology calculations.
The synthetic data (simulation) refers to a bead in an optical trap diffusing in 2D in a Newtonian fluid.
The real data refer to a 1.5um silica bead in an optical trap diffusing in water. Position data is obtained via video-particle tracking.
Both data series have been generated/acquired at 1 kHz.

The raw data is contained in the track_particle_output files.
However, the routine will take as input the Linked series (for simulation data) and the Dedrifted series (for real data).
If you want to play around with data manipulation (such as linking and dedrifting), try out the GUI_Data_Analysis_Standard routine.

The Gmoduli for the data series have been calculated with the following parameters:

T=298
Beta=5000
Dbead=1.5um (for real data)
     =2 um  (for simulation data)
Log-sampling=1.2 (for real data)
            =1.25 (for simulation data)



