#!/bin/bash
make
mv eco-cluster Testes
cd Testes
chmod +x eco_testes.sh
./eco_testes.sh
