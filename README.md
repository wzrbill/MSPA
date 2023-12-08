# Intro to MSPA
## I. Run a demo
We publish our demo [here](https://codeocean.com/capsule/8326543/tree).
## II. Run on dataset mentioned in our paper
### 1) Data preparation
1. Graph file format
   
We store our graph in *.gr
```c++
n m       // #nodes #edges
u1 v1 w1  // node id starts from 1 
u2 v2 w2  // in social graph w is 1
...
ux vy wm
``` 

2. Construct index
   
We offer the modified index construction source code in folder "/tools".
Execute the following bash scipt to compile the executable binary. 
```bash
bash run_compile_construct_index.sh
```
And it will generates "h2h_construct_index" and "pll_construct_index"
* For h2h index on road networks
```bash
./h2h_construct_index example-d.gr example-d.index
```
* For pll index on social networks
Note: we set the bit parallel as 64
```bash
./pll_construct_index example.gr example.index
```

3. Generate queryfile
   
The random generator is also offered in folder "/tools".
Execute the following bash scipt to compile the executable binary. 
```bash
bash run_compile_rand_gen.sh
```
And it will generates "rdgen".
Run the following commmand to generate a set of 50 queries in "example_9.Q"
```bash
./rdgen example.gr 9 50
```

4. Batch scripts and config files
   
We offer some batch scipt template easy to be modified to construct indices and query files.
We offer a "CMakeLists.txt" and "settings.json" if to debug in VSCode. Please comment "add_definitions(-DRECOVER_MSPATREE)" to turn off the retrieval function.
Please find the templates and config files in folder "/tools"


## 2) Folders
"-" means the rank of the folder.
```bash
-rootFolder
  --script
    ---example
  --src
  --TKDE_MSPA_DATA
    ---candidates
    ---graphs
    ---log
    ---queryfile
    ---result
  --tools
  --CMakeLists.txt
```
* script
  
We put the scipts to run all our algorithms for graph example.gr in the folder "script/example/".

* src
  
We put the source code of our project in src. And compile the project with "CMakeLists.txt".

* TKDE_MSPA_DATA
  
We put "example.gr" and "example.index" in the folder "TKDE_MSPA_DATA/graphs/example-d/".
We put the candidates selected in folder "TKDE_MSPA_DATA/candidates/".
We put the logs generated in folder "TKDE_MSPA_DATA/log/".
We put the queryfile generated in folder "TKDE_MSPA_DATA/log/".
We put the results generated in folder "TKDE_MSPA_DATA/result/".

* tools
  
We put some scipts and config files in the folder "tools/".

## 3) How to run

For each graph, run scripts for each algorithm.
Take NY for example.
Run "run_NY_advDP.sh" for advDP, |S| from 5 to 8, k=5.
Run "run_NY_apx_cand.sh" for Sim, SimApp, bet, mst, pcp, |S| from 5 to 8, k=5.
Run "run_NY_apx_cand_1e2_3.sh" for Sim, SimApp, bet, mst, pcp, |S| = 100, k=100 and |S| = 1000, k=1000.

## 4) Contact

Feel free to contact me at [wzrbill@e.gzhu.edu.cn](mailto:wzrbill@e.gzhu.edu.cn)
