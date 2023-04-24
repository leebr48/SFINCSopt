#!/bin/bash -l

/u/lebra/src/VMEC2000/_skbuild/linux-x86_64-3.10/cmake-build/build/bin/xvmec input.vmec &
wait $!
