#!/bin/bash
#改这个的话config也要改
aodfile="/Users/allenzhou/ALICE/alidata/23zzh_pass5/AO2D_3.root"
aodlist="input_data_local.txt"

# config_json="/Users/allenzhou/ALICE/alidata/hyperloop_setup/po/run_local/configuration.json"
config_json="configuration.json"

# check AO2D file
if [ ! -f "$aodfile" ]; then
  echo "Error: AOD file '$aodfile' not found!"
  exit 1
fi
# check json file
if [ ! -f "$config_json" ]; then
  echo "Error: AOD file '$config_json' not found!"
  exit 1
fi

# export FAIRMQ_SHM_DEFAULT_SIZE=4000000000
# export O2_DPL_TRANSPORT=zeromq
# export FAIRMQ_TRANSPORT=zeromq

# o2-analysis-cf-femtodream-pair-track-cascade --configuration json://${config_json} --aod-file ${aodfile} --shm-segment-size 1000000000 --severity debug

# o2-analysis-event-selection-service --configuration json://${config_json} --aod-file ${aodfile}|\
# o2-analysis-multcenttable --configuration json://${config_json} |\
# o2-analysis-propagationservice --configuration json://${config_json} |\
# o2-analysis-trackselection --configuration json://${config_json} |\
# o2-analysis-pid-tpc-service --configuration json://${config_json} |\
# o2-analysis-pid-tof-merge --configuration json://${config_json} |\
# o2-analysis-ft0-corrected-table --configuration json://${config_json} |\
# o2-analysis-cf-femtodream-producer --configuration json://${config_json} |\
# o2-analysis-cf-femtodream-pair-track-cascade --configuration json://${config_json} --shm-segment-size 1000000000 --severity debug

o2-analysis-pid-tof-merge --configuration json://${config_json} | o2-analysis-ft0-corrected-table --configuration json://${config_json} | o2-analysis-multcenttable --configuration json://${config_json} | o2-analysis-event-selection-service --configuration json://${config_json} | o2-analysis-propagationservice --configuration json://${config_json} | o2-analysis-pid-tpc-service --configuration json://${config_json} | o2-analysis-cf-pid-flow-pt-corr --configuration json://${config_json} | o2-analysis-cf-femtodream-producer --configuration json://${config_json} | o2-analysis-cf-femtodream-pair-track-cascade --configuration json://${config_json} | o2-analysis-trackselection --configuration json://${config_json} --aod-file @input_data_local.txt --readers 1 --aod-memory-rate-limit 471859200 --shm-segment-size 3758096384 --resources-monitoring 5 --resources-monitoring-dump-interval 600 --aod-max-io-rate 1000


# o2-analysis-mc-converter -b --configuration json://${config_json} | o2-analysis-pid-tof-merge -b --configuration json://${config_json} | o2-analysis-ft0-corrected-table -b --configuration json://${config_json} | o2-analysis-multcenttable -b --configuration json://${config_json} | o2-analysis-event-selection-service -b --configuration json://${config_json} | o2-analysis-propagationservice -b --configuration json://${config_json} | o2-analysis-pid-tpc-service -b --configuration json://${config_json} | o2-analysis-cf-pid-flow-pt-corr -b --configuration json://${config_json} | o2-analysis-cf-femtodream-producer -b --configuration json://${config_json} | o2-analysis-cf-femtodream-pair-track-cascade -b --configuration json://${config_json} | o2-analysis-trackselection -b --configuration json://${config_json} --aod-file @input_data.txt --readers 1 --aod-memory-rate-limit 471859200 --shm-segment-size 3758096384 --resources-monitoring 5 --resources-monitoring-dump-interval 600 --aod-max-io-rate 1000


# --timeframes-rate-limit 2 \ --pipeline 1
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
