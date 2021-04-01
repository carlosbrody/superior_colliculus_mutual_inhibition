If you want to make the figures, you need to move the scripts to the relevant directory so that imports work properly

cp Nature_Neuro_Figures/*.jl ../
cp Nature_Neuro_Figures/*.py ../
cp Nature_Neuro_Figures/EDF_9bcd.jl ../A_REWRITE

To make EDF 9b,c,d
open julia 1.5, in /A_REWRITE
include("EDF_9bcd.jl")

To make EDF 9e
open python3 in superior_colliculus_mutual_inhibition/
import EDF_9e

To make all other figures
open julia 0.6.4 in superior_colliculus_mutual_inhibition/
include("<name of script.jl>")

EDF_8a.jl generates a data file of example model trajectories.
EDF_8c.jl later uses that example file, or you can uncomment line 14 to make examples without running EDF_8a.jl


