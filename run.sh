#!/bin/bash

conda activate fairProject

cd apache-ignite-2.15.0-bin/bin

./ignite.sh ../config/ignite-config-1.xml &

./ignite.sh ../config/ignite-config-2.xml &

./ignite.sh ../config/ignite-config-3.xml &

./ignite.sh ../config/ignite-config-4.xml &