#!/bin/bash

mkdir Output-local
mkdir Output-avg

n=1000  # number of simulations

for i in $(seq 1 $n); do
    sleep 2
    echo "Running simulation $i"

    seed="313$i"
    src_file="activity-Fig4.c"
    run_file="run_$i.c"
    exe_file="sim_$i"

    cp "$src_file" "$run_file"
    sed -i "s/seed999/$seed/g" "$run_file"
    sed -i "s/file9/$i/g" "$run_file"

    gcc "$run_file" -lm -o "$exe_file"

    rem=$(($i % 100))
    if [ $rem -eq 0 ]; then
        echo "Running $exe_file in foreground"
        ./"$exe_file"
    else
        echo "Running $exe_file in background"
        nohup ./"$exe_file" > "output_$i.log" 2>&1 &
    fi
done

wait
echo "All simulations launched."

# ---- CLEANUP SECTION ----
echo "Cleaning up temporary files..."

#rm -f group_layout*.txt run_*.c sim_* output_*.log
rm -f run_*.c sim_* output_*.log group*
echo "Cleanup done."

