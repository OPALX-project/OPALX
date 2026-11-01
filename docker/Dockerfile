# syntax=docker/dockerfile:1
#
# OPALX container images (CPU / OpenMP).
#
# Build from the repository root:
#   docker build -f docker/Dockerfile -t opalx:cpu --build-arg GIT_SHA=$(git rev-parse HEAD) .
#   docker build -f docker/Dockerfile --target dev -t opalx-dev .
#
# Stages:
#   dev      toolchain only (compilers, MPI, CMake); source is mounted at run time
#   build    compiles OPALX from the build context
#   runtime  the opalx executable and its runtime libraries (default target)

ARG UBUNTU_VERSION=24.04

# ------------------------------------------------------------------------------
# dev: same packages as the CPU OpenMP CI job
# ------------------------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS dev

ENV DEBIAN_FRONTEND=noninteractive \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake git ca-certificates \
      libopenmpi-dev openmpi-bin python3 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# ------------------------------------------------------------------------------
# build: dependencies (IPPL, Kokkos, Heffte, HDF5, H5hut, LAPACK) are fetched
# during configure, so this stage needs network access
# ------------------------------------------------------------------------------
FROM dev AS build

ARG BUILD_TYPE=Release
ARG GIT_SHA=unknown
# Parallel compile jobs; empty means one per CPU
ARG JOBS=

COPY . /opt/opalx/src

# .git is not part of the build context, so write the revision header that
# tools/gen_OPALXrevision would otherwise generate from git
RUN echo "#define GIT_VERSION \"${GIT_SHA}\"" > /opt/opalx/src/src/OPALXrevision.h \
 && cmake -S /opt/opalx/src -B /opt/opalx/build \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DPLATFORMS=OPENMP \
 && cmake --build /opt/opalx/build -j "${JOBS:-$(nproc)}"

# License texts of OPALX and of the libraries compiled into it
RUN python3 /opt/opalx/src/docker/collect_licenses.py \
      /opt/opalx/src /opt/opalx/build /opt/opalx/licenses "${GIT_SHA}"

# ------------------------------------------------------------------------------
# runtime
# ------------------------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS runtime

ARG GIT_SHA=unknown

LABEL org.opencontainers.image.title="OPALX" \
      org.opencontainers.image.description="OPALX, CPU (OpenMP) build" \
      org.opencontainers.image.source="https://github.com/OPALX-project/OPALX" \
      org.opencontainers.image.revision="${GIT_SHA}" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      openmpi-bin libgfortran5 libgomp1 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/opalx/build/src/opalx /usr/local/bin/opalx
COPY --from=build /opt/opalx/licenses /usr/share/doc/opalx
COPY --chmod=0755 docker/entrypoint.sh /usr/local/bin/opalx-entrypoint

# Fail the image build if the executable needs a library that is not installed
RUN ! ldd /usr/local/bin/opalx | grep 'not found'

# Containers often run as root; without the OMPI_ALLOW_* variables mpirun refuses to start
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

WORKDIR /data
ENTRYPOINT ["opalx-entrypoint"]
CMD ["--help"]
