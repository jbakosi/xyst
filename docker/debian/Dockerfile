################################################################################
# vim: filetype=dockerfile:
#
# \file      docker/debian
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.,
#            2022-2025 J. Bakosi
#            All rights reserved. See the LICENSE file for details.
# \brief     Install xyst external libraries ready for CI
# \details
#   This image is intended to prepare a docker image that contains the external
#   libraries pre-installed. The image is parked at
#
#   https://hub.docker.com/repository/docker/jbakosi/
#
#   which can then be pulled by continuous integration (CI) to build and test.
#   This ensures the external libs do not have to be rebuilt any time the test
#   suite needs to run which reduces load on CI.
#
#   Build image, squash, and push, after commenting Clone, Build, and Test:
#   docker build -f debian -t xyst-gnu .
#   docker-squash --cleanup xyst-gnu -t jbakosi/xyst-gnu
#   docker push jbakosi/xyst-gnu
#
#   Test in non-SMP mode, after uncommenting Clone, Build, and Test:
#   docker build -f debian -t xyst-gnu .
#
# \note No SMP build because that requires own openmpi to avoid openmpi runtime
#   errors in smp mode, see debian-clang for SMP mode.
#
################################################################################

FROM debian:stable

# Install system-wide prerequisites
RUN apt-get update -y && apt-get install -y git cmake g++ gfortran ninja-build bzip2 zlib1g-dev libnetcdf-mpi-dev

# Install Charm++
WORKDIR /etc/
#RUN git sparse-checkout init
#RUN git sparse-checkout set src cmake CMakeLists.txt build buildcmake &&
#RUN git clone https://github.com/charmplusplus/charm && cd charm && git checkout v7.0.0-rc2
RUN git clone --filter=blob:none -n https://github.com/charmplusplus/charm; cd charm; git sparse-checkout set src cmake contrib tests; git checkout v7.0.0-rc2
WORKDIR /etc/charm/
RUN mkdir nonsmp && cd nonsmp && cmake -DSMP=off -DTARGET=charm -DNETWORK=mpi -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=/etc/charm-non-smp .. && make -sj$(grep -c processor /proc/cpuinfo) install
RUN rm -rf charm

# Clone
ARG COMMIT
RUN git clone --filter=blob:none -n https://codeberg.org/xyst/xyst.git && cd xyst && git sparse-checkout set docker src tests cmake && git checkout $COMMIT && git log -1 HEAD
# # Build
RUN cd xyst && mkdir -p build && cd build && cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -GNinja -DCHARM_INSTALL_DIR=/etc/charm-non-smp -DRUNNER=mpirun -DRUNNER_NCPUS_ARG=-n -DRUNNER_ARGS="--bind-to none --allow-run-as-root --oversubscribe" -DPOSTFIX_RUNNER_ARGS=+setcpuaffinity ../src && ninja

# Test
RUN cd xyst/build && mpirun --oversubscribe --allow-run-as-root -n $(grep -c processor /proc/cpuinfo) Main/unittest -q && ctest -j$(grep -c processor /proc/cpuinfo) --output-on-failure 



RUN install xyst/build/Main/inciter /bin
RUN install xyst/build/Main/meshconv /bin
# RUN rm -rf /xyst

RUN addgroup --gid 1001 xyst_user
RUN adduser --uid 1001 --gid 1001 --disabled-password --gecos "" xyst_user
RUN echo 'xyst_user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers
ADD motd /etc/motd
ADD bash_config /home/xyst_user/.bashrc

VOLUME /home/xyst_user/data
USER xyst_user:xyst_user
WORKDIR /home/xyst_user
ENTRYPOINT bash
