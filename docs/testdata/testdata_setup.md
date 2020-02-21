# Test data download and setup on and off the cluster

#### Requirements:

this guide assumes that you have
1. access to wget

### Download the test data

All test data is stored on Azure in public storage containers. To download test
data sets for different parts of the pipeline, simply use wget + the appropriate URL
for the test data of your choice. 

### pull test data from azure storage

Test data for the different parts of the wgs pipeline is available in the azure storage container `wgstestdata`.

Below are pipeline-specific commands to pull testdata

*alignment:*

```wget https://wgstestsets.blob.core.windows.net/testsets/alignment.tar```

*snv_calling:*

```wget https://wgstestsets.blob.core.windows.net/testsets/snv_calling.tar```

*sv_calling:*

```wget https://wgstestsets.blob.core.windows.net/testsets/sv_calling.tar```

*copy number calling:*

```wget https://wgstestsets.blob.core.windows.net/testsets/cna.tar```

Each tarball contains test data and a premade input.yaml that can be fed into the pipeline.  
