#Downloads the data into a zip file
Curl -L0 http://tuvalu.santafe.edu/~aaronc/facultyhiring/replicationData_all.zip > data/replicationData_all.zip

#touch and unzip the data
touch data/replicationData_all.zip
unzip data/replicationData_all.zip -d data/
