Simulations with ForCEEPS:
- I generate a new climate file (e.g. "forceps.climate.bern+1").

- In the CAPSIS directory (data/forceps/leo), I have a setup file for every site. 
e.g. 
	"forceps.setup_Bern_Temp". In this file I have to change the climate file name:
	from: defaultClimateFileName = forceps.climate.Bern.txt
	to: defaultClimateFileName = forceps.climate.Bern+1.txt
DO NOT FORGET to do it for every simulation.

- I also have a command file for every site:
e.g.
	"cmd2_Bern_Temp.txt"
This cmd file is necessarily linked to one setup file. I thus have to launch a simulation for every temperature condition, with the corresponding climate file;
I cannot automate it.
I can comment the other lines.