Under Windows:

Je vais sur mon dossier Capsis - shift + clic droit - ouvrir PowerShell ici

~\Capsis_dev$ #(= go into Capsis environment)
sh capsis.sh -p script # runs the capsis.sh script found in /
forceps.myscripts.SimulationPG # calls the simulation file found in /src/forceps/myscripts/SimulationPG.java .
#It reads, interprets and runs one simulation per line on the command file and defines the style of the export. 

------
example: run command: cmd2_Bever_random_x.txt
capsis -p script forceps.myscripts.SimulationPG data\forceps\leo\cmd2_Bever_decreasing.txt
ou plutôt \.capsis -p script forceps.myscripts.SimulationPG data\forceps\leo\cmd2_Bever_decreasing.txt

# capsis1 -p script forceps.myscripts.SimulationPG data\forceps\leo\cmd2_Bever_decreasing.txt