./clean.sh
rm -rf output
python3 build_sdc.py
mkdir output
cd output
./../job*/run.sh
