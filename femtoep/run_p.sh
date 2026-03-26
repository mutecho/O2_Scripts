aodfile="/Users/allenzhou/ALICE/alidata/23zzh_pass5/AO2D_3.root"

# check AO2D file
if [ ! -f "$aodfile" ]; then
  echo "Error: AOD file '$aodfile' not found!"
  exit 1
fi

o2-analysis-event-selection-service --configuration json://config.json --aod-file ${aodfile} |\
o2-analysis-multcenttable --configuration json://config.json |\
o2-analysis-propagationservice --configuration json://config.json |\
o2-analysis-trackselection --configuration json://config.json |\
o2-analysis-pid-tpc-service --configuration json://config.json |\
o2-analysis-pid-tof-merge --configuration json://config.json |\
o2-analysis-qvector-table --configuration json://config.json |\
o2-analysis-ft0-corrected-table --configuration json://config.json |\
o2-analysis-cf-femtodream-producer --configuration json://config.json


# o2-analysis-timestamp --configuration json://config.json --aod-file ${aodfile} |\
# o2-analysis-event-selection --configuration json://config.json |\
# o2-analysis-multiplicity-table --configuration json://config.json |\
# o2-analysis-centrality-table --configuration json://config.json |\
# o2-analysis-track-propagation --configuration json://config.json |\
# o2-analysis-trackselection --configuration json://config.json |\
# o2-analysis-pid-tpc-service --configuration json://config.json |\
# o2-analysis-pid-tof-merge --configuration json://config.json |\
# o2-analysis-lf-cascadebuilder --configuration json://config.json |\
# o2-analysis-lf-lambdakzerobuilder --configuration json://config.json |\
# o2-analysis-qvector-table --configuration json://config.json |\
# o2-analysis-ft0-corrected-table --configuration json://config.json |\
# o2-analysis-cf-femtodream-producer --configuration json://config.json
# o2-analysis-cf-femtodream-pair-track-track --configuration json://config.json

# o2-analysis-collision-converter --configuration json://config.json |\
# o2-analysis-lf-lambdakzerobuilder --configuration json://config.json |\
# o2-analysis-tracks-extra-v002-converter --configuration json://config.json |\
# o2-analysis-timestamp --configuration json://config.json |\
# o2-analysis-propagationservice --configuration json://config.json |\
# o2-analysis-trackextension --configuration json://config.json |\
# o2-analysis-event-selection-service --configuration json://config.json --aod-file ${aodfile} |\
# o2-analysis-pid-tof-full --configuration json://config.json |\
# o2-analysis-pid-tof-base --configuration json://config.json |\
# o2-analysis-pid-tpc --configuration json://config.json |\
# o2-analysis-pid-tpc-base --configuration json://config.json |\
