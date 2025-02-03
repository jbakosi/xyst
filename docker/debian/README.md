# This is README of Xyst Docker Image

It was meant to simplify the obtaining process of Xyst, 

## Building process

This tutorial assumes that you have downloaded Docker

1. Clone the repo (if you did not already):
``` bash
git clone git@codeberg.org:xyst/xyst.git 
cd xyst
```

2. Build the docker image

``` bash
docker build docker/debian -t xyst
```

3. Run the installed docker image, with mounted files from your working directory `${PWD}`

``` bash
docker run --rm -it -v ${PWD}:/home/xyst_user/data xyst
```

