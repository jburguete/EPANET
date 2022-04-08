#!/usr/bin/env bash
time for i in ../*/{,*/{,*/}}*.inp; do ./epanet-old $i $i.out.old; done
time for i in ../*/{,*/{,*/}}*.inp; do ./epanet $i $i.out; done
for i in ../*/{,*/{,*/}}*.out; do diff $i $i.old; done
