programs=("split_dfs" "split_bfs" "cand")
for pg in "${programs[@]}"; do
    echo "Running $pg"
    bash "script/Email/run_Email_$pg.sh" &
done
