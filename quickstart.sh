#!/bin/bash

# A script for executing all steps of the arora-ge-ntru program with default 
# parameters given below, which can optionally be changed on the command line.
#
# Usage: in the top directory do
#
#    ./quickstart [-n n] [-q q] [-k k] [-c c] [-r r] [-s s]
#

n=16      # degree
q=31      # prime modulus
keys=10   # number of samples
coeffs=2  # 2 for binary coefficients, 3 for ternary (applies to both numerator and denominator)
ring=1    # 1 for x^n - 1, 2 for x^n + 1.
seed=1    # fixed seed for rng.

while getopts :n:q:k:c:r:s: flag
do
  case "${flag}" in
    n) n=${OPTARG};;
    q) q=${OPTARG};;
    k) keys=${OPTARG};;
    c) coeffs=${OPTARG};;
    r) ring=${OPTARG};;
    s) seed=${OPTARG};;

    #https://stackoverflow.com/questions/22909253
    \?) valid=0
        echo "An invalid option has been entered: $OPTARG"
        exit 1
        ;;

    :)  valid=0
        echo "The additional argument for option $OPTARG was omitted."
        exit 1
        ;;
  esac
done

rm -rf _data
mkdir _data
cd _data

pkout="pk"
skout="sk"
sysout="sys"
resout="res"

#../apps/bin/arora-ge-ntru $n $q -s $seed -c $coeffs -r $ring --verbose all -k $keys

printf "[quickstart] Generating keys...\n"
time ../build/apps/bin/arora-ge-ntru $n $q -s $seed -c $coeffs -r $ring --verbose keygen -k $keys --pk_output=$pkout --sk_output=$skout
printf "[quickstart] Done.\n\n"

if test -f $pkout;
then
  printf "[quickstart] Building linear system...\n"
  time ../build/apps/bin/arora-ge-ntru $n $q -c $coeffs -r $ring --verbose system -i $pkout -o $sysout
  printf "[quickstart] Done.\n\n"
fi

if test -f $sysout;
then
  printf "[quickstart] Recovering secret...\n"
  time ../build/apps/bin/arora-ge-ntru $n $q -c $coeffs -r $ring --verbose recover -i $sysout -o $resout
  printf "[quickstart] Done.\n\n"
fi

if test -f $resout;
then
  printf "[quickstart] Verifying result...\n"
  time ../build/apps/bin/arora-ge-ntru $n $q -c $coeffs -r $ring --verbose verify --sk_input1 $resout --sk_input2 $skout
  printf "[quickstart] Done.\n\n"
fi
