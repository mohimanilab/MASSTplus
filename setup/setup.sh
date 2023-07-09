cp binaries/* clustered-index
cp binaries/* unclustered-index

cd clustered-index

# Download clustered index
wget https://dl.dropboxusercontent.com/s/0jhcpmkoxd9dlqi/index-clustered.tar.gz?dl=0
tar -xvf index-clustered.tar.gz
mv config.txt library
rm index-clustered.tar.gz

cd ../unclustered-index

# Download unclustered index
wget -O aa.tar.gz https://dl.dropboxusercontent.com/s/i3j05du4p3dd9ud/aa.tar.gz\?dl\=0
wget -O ai.tar.gz https://dl.dropboxusercontent.com/s/yuzhqefowbui4lx/ai.tar.gz?dl=0
wget -O library.tar.gz https://dl.dropboxusercontent.com/s/6fhizr4byqyp261/library.tar.gz?dl=0
wget -O ah.tar.gz https://dl.dropboxusercontent.com/s/8dj64y0hrzivnop/ah.tar.gz?dl=0
wget -O ag.tar.gz https://dl.dropboxusercontent.com/s/48wfaxua2pcsf0m/ag.tar.gz?dl=0
wget -O af.tar.gz https://dl.dropboxusercontent.com/s/sbs7wxy3tpr87rk/af.tar.gz?dl=0
wget -O ae.tar.gz https://dl.dropboxusercontent.com/s/fytzx82c4pqz2mo/ae.tar.gz?dl=0
wget -O ad.tar.gz https://dl.dropboxusercontent.com/s/r8k9enr0wj8f39r/ad.tar.gz?dl=0
wget -O ac.tar.gz https://dl.dropboxusercontent.com/s/5a0y83jter7hpu2/ac.tar.gz?dl=0
wget -O ab.tar.gz https://dl.dropboxusercontent.com/s/riewsuvl5gb0cti/ab.tar.gz?dl=0

tar -xvf library.tar.gz
rm library.tar.gz
mv a*.tar.gz library/partitions
cd library/partitions
tar -xvf aa.tar.gz
tar -xvf ab.tar.gz
tar -xvf ac.tar.gz
tar -xvf ad.tar.gz
tar -xvf ae.tar.gz
tar -xvf af.tar.gz
tar -xvf ag.tar.gz
tar -xvf ah.tar.gz
tar -xvf ai.tar.gz
rm *.tar.gz
