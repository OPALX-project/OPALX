#!/bin/bash

# Copyright (c) 2012, Chris Rogers
# All rights reserved.
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met: 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
# 3. Neither the name of STFC nor the names of its contributors may be used to 
#    endorse or promote products derived from this software without specific 
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.

here=`pwd`/
vers=1.7.0

check_dependencies() {
    dep_list="wget unzip md5sum"
    echo "Checking that dependencies are installed"
    for dep in ${dep_list}; do
        echo "Checking dependency ${dep}..."
        ${dep} --help > /dev/null
        if [ $? -ne 0 ]; then
            echo "FAIL: require ${dep} to run script"
            exit 1
        fi
        echo "                        ok"
    done
}

download_gtest() {
    source_dir=${here}/src
    gtest_dir=gtest-${vers}
    file=${gtest_dir}.zip
    url=http://googletest.googlecode.com/files/${file}
    md5="2d6ec8ccdf5c46b05ba54a9fd1d130d7  ${file}"

    echo "Get a new copy of the gtest library if required"
    mkdir -p ${source_dir}
    cd ${source_dir}
    if [ -e ${source_dir}/${gtest_dir}/${file} ]; then
        mv ${source_dir}/${gtest_dir}/${file} ${source_dir}
    fi
    echo "${md5}" | md5sum -c
    if [ $? -ne 0 ]; then
        wget ${url}
    fi
    echo "${md5}" | md5sum -c
    if [ $? -ne 0 ]; then
        echo "FAIL: failed to download gtest library"
        exit 1
    fi
    echo "Clean existing gtest build"
    rm -rf ${gtest_dir}
    echo "Unzipping gtest source"
    unzip ${file}
    # we clean up by moving the gtest zip file into the gtest build directory
    mv ${file} ${gtest_dir}
}

build_gtest() {
    echo "Build the gtest library"
    cd ${gtest_dir}
    cmake .
    #./configure --prefix=${install_dir}
    make
    if [ $? -ne 0 ]; then
        echo "FAIL: failed to build gtest"
        exit 1
    fi
}

install_gtest() {
    echo "Installing the gtest library"
    install_dir=${here}
    inc_dir=${install_dir}/include/
    lib_dir=${install_dir}/lib/

    # I toyed with the idea of doing a cleanup here, decided against it in the
    # end; so we overwrite but don't clean existing
    mkdir -p ${inc_dir}
    cp -rf ${source_dir}/${gtest_dir}/include/* ${install_dir}/include/
    mkdir -p ${install_dir}/lib/
    cp -r ${source_dir}/${gtest_dir}/*.a ${install_dir}/lib/
}

main() {
    check_dependencies
    download_gtest
    build_gtest
    install_gtest
}

main
