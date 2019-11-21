# BioToolBox Docker

Installing complex software projects and keepting them up to date can be an onerous 
and frustrating experience. BioToolBox is not immune to this.

Docker images are becoming an increasingly important part of bioinformatics as a 
solution to this problem. By using a simple Docker image, one can keep 
these software tools contained within their own environment with simple installation.

Here is a Dockerfile for users to build their own BioToolBox image.

## Building the Docker image

On your machine (local or server), you will need [Docker](https://www.docker.com) (or 
a compatible equivalent) installed and running.

This is based on an Ubuntu. It follows the 
[advanced installation](../docs/AdvancedInstallation.md) guide. 

You will need to download the `Dockerfile` from here or reference the raw URL to it. 

	$ docker build -t biotoolbox:1.67 path/to/biotoolbox/Docker/

You can run the image interactively to check things.
	
	$ docker run -it biotoolbox:1.67
	root@b1296410a011:/# which get_datasets.pl
	/usr/bin/get_datasets.pl
	root@b1296410a011:/# exit
	
The image has an empty data directory at `/data`. You can bind a path on your 
machine to this data directory in the image.

	$ docker run -it -v /path/to/my/data:/data biotoolbox:1.67

You can find more information about using docker images elsewhere online.


