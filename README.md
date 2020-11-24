# jgi_meta_wdl

![alt text](images/1.png "Titll")
![alt text](images/2.png "Titll")

```
# install wget.
sudo yum install wget -y

# Download and install anaconda.
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/miniconda3
source miniconda3/etc/profile.d/conda.sh && conda activate

# Install git and cromwell.
conda install -c conda-forge git cromwell -y

# Change directory to /tmp and download code.
cd /tmp; git clone https://code.jgi.doe.gov/BFoster/jgi_meta_wdl

# Fetch test data.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR787/004/SRR7877884/SRR7877884_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR787/004/SRR7877884/SRR7877884_1.fastq.gz

# Interleave fastq data.
docker run --volume $PWD:/data -w /data bryce911/bbtools:38.86 reformat.sh in=SRR7877884_1.fastq.gz in2=SRR7877884_2.fastq.gz out=SRR7877884.fastq.gz

# Download rqcfilter dataset (~2hours).
mkdir data; cd data; wget -O - http://portal.nersc.gov/dna/metagenome/assembly/rqcfilter/RQCFilterData.tar | tar -xf - ; cd ..

# Make an inputs.json file containing the the path to the test data.
echo '{"metagenome_filtering_assembly_and_alignment.input_files": ["/tmp/SRR7877884.fastq.gz"]}' > inputs.json

# Run pipeline.
cromwell -Dconfig.file=jgi_meta_wdl/local.conf run -i inputs.json jgi_meta_wdl/metagenome_filtering_assembly_and_alignment.wdl
```

![alt text](images/4.png "Title")