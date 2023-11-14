#!/bin/bash

echo "100 millions in one go"
date +%T:%N
./trivial 0 100000000
date +%T:%N

echo "100 Millions in two go"
date +%T:%N

./trivial 0 50000001 &
./trivial 50000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in three go"
date +%T:%N

./trivial 0 33000001 &
./trivial 33000001 66000001 &
./trivial 66000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in four go"
date +%T:%N

./trivial 0 25000001 &
./trivial 25000001 50000001 &
./trivial 50000001 75000001 &
./trivial 75000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in five go"
date +%T:%N

./trivial 0 20000001 &
./trivial 20000001 40000001 &
./trivial 40000001 60000001 &
./trivial 60000001 80000001 &
./trivial 80000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in six go"
date +%T:%N

./trivial 0 16666667 &
./trivial 16666667 33333334 &
./trivial 33333334 50000001 &
./trivial 50000001 66666668 &
./trivial 66666668 83333335 &
./trivial 83333335 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in seven go"
date +%T:%N

./trivial 0 14290000 &
./trivial 14290000 28580000 &
./trivial 28580000 42870000 &
./trivial 42870000 57160000 &
./trivial 57160000 71450000 &
./trivial 71450000 85740000 &
./trivial 85740000 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in eight go"
date +%T:%N

./trivial 0 12500001 &
./trivial 12500001 25000001 &
./trivial 25000001 37500001 &
./trivial 37500001 50000001 &
./trivial 50000001 62500001 &
./trivial 62500001 75000001 &
./trivial 75000001 87500001 &
./trivial 87500001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in nine go"
date +%T:%N

./trivial 0 11000001 &
./trivial 11000001 22000001 &
./trivial 22000001 33000001 &
./trivial 33000001 44000001 &
./trivial 44000001 55000001 &
./trivial 55000001 65000001 &
./trivial 65000001 76000001 &
./trivial 76000001 86000001 &
./trivial 86000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N

echo "100 Millions in ten go"
date +%T:%N

./trivial 0 11000001 &
./trivial 10000001 20000001 &
./trivial 20000001 30000001 &
./trivial 30000001 40000001 &
./trivial 40000001 50000001 &
./trivial 50000001 60000001 &
./trivial 60000001 70000001 &
./trivial 70000001 80000001 &
./trivial 80000001 90000001 &
./trivial 90000001 100000000 &

wait  # don't measure the time before all background tasks finished
date +%T:%N
