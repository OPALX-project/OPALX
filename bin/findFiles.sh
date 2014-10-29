#!/bin/bash

find src -regex '.*\.cpp\|.*\.h' -printf '%h/%f ' |perl bin/mkver2
