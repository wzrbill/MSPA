# programs=("split_dfs" "7" "8" "9")
programs=("1e2_3_cand" "1e2_3_fs")
for pg in "${programs[@]}"; do
    echo "Running $pg"
    bash "script/NY/run_NY_$pg.sh" &
done
