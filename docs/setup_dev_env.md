


1. create working dir and clone code

```
mkdir /devsrc
cd /devsrc
git clone https://github.com/shahcompbio/wgs.git
```

2. build docker container

save the following in `dockerfile` file in `/devsrc`
```
# build on top of out base image
FROM wgspipeline/python_wgs:v0.0.1

# Install any needed packages specified in requirements.txt
RUN rm -rf /opt/conda/lib/python2.7/site-packages/pypeliner* /opt/conda/lib/python2.7/site-packages/wgs* /opt/conda/lib/python2.7/site-packages/biowrappers*
RUN pip install git+https://github.com/shahcompbio/pypeliner.git@master
RUN pip install git+https://bitbucket.org/aroth85/biowrappers.git@singlecell
RUN pip install dill
# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME wgs

ENV PYTHONPATH /devsrc/wgs

# Run app.py when the container launches
CMD ["wgs"]
```

build:
`docker build -t wgs .`

3. test the install

```
(base) root@scdna-dev-node:/devsrc# docker run -v /devsrc:/devsrc -w /devsrc -it wgs bash
(base) root@a8b040c372e3:/devsrc# which python
/opt/conda/bin/python
(base) root@a8b040c372e3:/devsrc# python -c "import wgs; print wgs.__file__"
/devsrc/wgs/wgs/__init__.pyc
```


