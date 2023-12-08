programs=("split_dfs" "split_bfs" "cand")
for pg in "${programs[@]}"; do
    echo "Running $pg"
    bash "script/NY/run_NY_$pg.sh" &
done
