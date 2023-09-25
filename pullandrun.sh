#!/bin/bash

cd dETdeta
git pull
cp build_hists.C ..
cd ..
bash build_script.sh
