# Rstudio_playground

### Info 
https://rocker-project.org/images/versioned/rstudio.html
https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/rstudio_devel.Dockerfile
Extending images    https://rocker-project.org/use/extending.html

# Build image

 Build Image:

```bash
docker build --progress=plain -f rocker/Dockerfile_tidyverse_verse_dockering --tag rocker_verse_rna-seq .
```

Explanation:

cd rocker/: This command changes the current directory to the rocker/ directory, where the Dockerfile and other relevant files are located.

docker build: This command is used to build a Docker image from a Dockerfile.

--progress=plain: This option specifies the format for the progress output during the build process. Setting it to plain ensures that only plain text progress is displayed, which can be helpful for scripting or automation.

-f rocker/Dockerfile_tidyverse_verse_dockering: This option specifies the path to the Dockerfile to use for building the image. In this case, it points to a Dockerfile named Dockerfile_tidyverse_verse_dockering located in the rocker/ directory.

--tag rocker_verse_rna-seq: This option assigns a tag to the built image, allowing it to be easily referenced and identified. Here, the tag rocker_verse_rna-seq is assigned to the image.

# Run Image:

bash
Copy code
docker run --rm -ti -e PASSWORD=password -p 127.0.0.1:8787:8787 -v $(pwd):/home/rstudio -it rocker_verse_rna-seq
Explanation:

docker run: This command is used to run a Docker container based on a specified image.

--rm: This option specifies that the container should be removed automatically after it exits. This helps to clean up resources and avoid cluttering the system with unused containers.

-ti: These options allocate a pseudo-TTY and keep STDIN open even if not attached, allowing interactive terminal access to the container.

-e PASSWORD=password: This option sets the environment variable PASSWORD inside the container to the value password. This is used to set the password for accessing RStudio within the container.

-p 127.0.0.1:8787:8787: This option maps port 8787 on the host machine to port 8787 inside the container. This enables access to RStudio Server running within the container via a web browser on http://localhost:8787.

-v $(pwd):/home/rstudio: This option mounts the current working directory ($(pwd)) on the host machine to the /home/rstudio directory inside the container. This allows files and data from the host to be accessed and manipulated within the container.

-it: These options allocate an interactive terminal session with the container.

rocker_verse_rna-seq: This is the name of the Docker image to run as a container. It refers to the image built in the previous step.





docker run -d -p 127.0.0.2:8787:8787 -v $(pwd):/home/rstudio -e DISABLE_AUTH=true -it iruizdea/rocker-verse-tinytex-rmarkdown
docker run --rm -v $(pwd):/home/rstudio -e USER=rstudio -e DISABLE_AUTH=true -it iruizdea/rocker-verse-tinytex-rmarkdown R

docker run --rm -v $(pwd):/home/rstudio -e USER=rstudio -e DISABLE_AUTH=true -it iruizdea/rocker-verse-tinytex-rmarkdown bash
```

# Clean out any temp or archive not needed to the container itself. From the strudio terminal or with exec

```bash
docker run --rm -v $(pwd):/home/rstudio -e USER=rstudio -e DISABLE_AUTH=true -it iruizdea/rocker-verse-tinytex-rmarkdown bash

rm -r tmp/RtmpCsOuJg/
strip /usr/local/lib/R/site-library/*/libs/*.so
rm -rf /var/lib/apt/lists/*

# From another tab to save a running image

docker ps
CONTAINER ID   IMAGE                                     COMMAND   CREATED         STATUS         PORTS      NAMES
1369933fdf66   iruizdea/rocker-verse-tinytex-rmarkdown   "bash"    6 minutes ago   Up 6 minutes   8787/tcp   reverent_brown

# Now to save this version of the image, in the new terminal window type:

docker commit -m "remove tmp and so* to slim" 1369933fdf66 iruizdea/rocker-verse-tinytex-rmarkdown
```

## Tag & push the Docker image to docker hub

```bash
docker images
docker tag iruizdea/rocker-verse-tinytex-rmarkdown iruizdea/rocker-verse-tinytex-rmarkdown:latest

# Replace your docker hub account Id
docker login -u iruizdea


docker push iruizdea/rocker-verse-tinytex-rmarkdown:latest
```

## Save the docker image into a portable .tar file.
```bash
docker save -o rocker-verse-tinytex-rmarkdown.tar iruizdea/rocker-verse-tinytex-rmarkdown
```

##  Remove all Docker images
```bash
docker system prune -a --volumes
```

# Download and run docker

```bash
docker pull iruizdea/rocker-verse-tinytex-rmarkdown
docker run --rm -ti -e PASSWORD=password -p 127.0.0.1:8788:8787 -v $(pwd):/home/rstudio iruizdea/rocker-verse-tinytex-rmarkdown
```

## Include on the RedHat Quai.io repo

## Export docker to singularity  

```bash
# Point to docker container location 
singularity build fastqc_0.12.1_orad_2.6.1.sif docker-archive://fastqc_0.12.1_orad_2.6.1.tar
```
## Test interactive run

```bash
singularity shell  fastqc_0.12.1_orad_2.6.1.sif 
cat <(orad -c --raw test.fastq.ora ) | fastqc  stdin:Yuuujuuu

```


 Warning message:
#9 21.78 packages ‘airway’, ‘SummarizedExperiment’, ‘DESeq2’, ‘DelayedArray’, ‘vsn’, ‘apeglm’, ‘AnnotationDbi’, ‘org.Hs.eg.db’, ‘clusterProfiler’, ‘enrichplot’ are not available for this version of R



docker build --progress=plain -f Dockerfile_verse --tag iruizdea/rocker-verse-tinytex-rmarkdown .
docker build --progress=plain -f Dockerfile_rbase --tag iruizdea/r-base.rstudio .
docker build Dockerfile_verse --progress=plain --tag iruizdea/rocker-verse-tinytex-rmarkdown


# Works
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 -v $(pwd):/home/rstudio -it iruizdea/rocker-verse-tinytex-rmarkdown
docker run --rm -ti -e PASSWORD=password -p 127.0.0.1:8787:8787 -v $(pwd):/home/rstudio -it iruizdea/rocker-verse-tinytex-rmarkdown
docker run --rm -ti -p 127.0.0.1:8787:8787 iruizdea/rocker-verse-tinytex-rmarkdown


# Create an install inside some packages that fails on the Docker file
https://jsta.github.io/r-docker-tutorial/01-what-and-why.html



docker build -f Dockerfile_tidyverse_verse --progress=plain --tag iruizdea/rocker.tidyverse.verse .
docker run --rm -ti -e PASSWORD=password -p 127.0.0.1:8787:8787 -v $(pwd):/home/rstudio -it iruizdea/rocker.tidyverse.verse

docker run --rm -ti -e PASSWORD=password -p 127.0.0.1:8799:8787 -v $(pwd):/home/rstudio -it iruizdea/rocker.tidyverse.class.bioclite

# Commit Docker container to include installed packages
https://jsta.github.io/r-docker-tutorial/03-install-packages.html

To save this specific version of the image we need to find this containers specific hash. We can see this by typing the following command in the new terminal window, and it will list all running Docker containers:

```
docker ps
``` 

CONTAINER ID   IMAGE                                      COMMAND   CREATED             STATUS             PORTS                      NAMES
8e5bd8b58fb4   iruizdea/rocker.tidyverse.class.bioclite   "/init"   56 minutes ago      Up 56 minutes      127.0.0.1:8799->8787/tcp   sweet_mccarthy
a6119b9170c0   iruizdea/rocker.tidyverse.class            "/init"   About an hour ago   Up About an hour   127.0.0.1:8795->8787/tcp   nice_goldstine
385803535e5f   iruizdea/rocker.tidyverse.class            "/init"   About an hour ago   Up About an hour   127.0.0.2:8791->8787/tcp   infallible_torvalds

# Clean out any temp or archive not needed to the container itself
rm -r tmp/RtmpCsOuJg/
strip /usr/local/lib/R/site-library/*/libs/*.so
rm -rf /var/lib/apt/lists/*

# Now to save this version of the image, in the new terminal window type:
 docker commit -m "install with bioclite and apt-get" 8e5bd8b58fb4 iruizdea/rocker.tidyverse.class.bioclite

