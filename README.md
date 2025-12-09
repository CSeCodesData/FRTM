# FRTM
## Files
- Codes: the directory of codes
- Examples: examples.zip
- Datasets: BJDATA.zip EURDATA.zip EURDATA.z01 SYNDATA.zip LANDATA.zip LANDATA.z01 LANDATA.z02 LANDATA.z03

## The Environment
- IDE: VS2022
- Operation System: win10 x64
***

## Parameters of The Execution Program
- v : The number of nodes
- e : The number of edges
- t : The number of snapshots
- i : File name of input temporal graph(default: example-temporal-graph.txt, required: information of all edges are sorted by timestamp, and the interval of temporal graph is [0,t-1])
- f : File name of output(default: example-output.txt)
- k : File name of frequency threshold(default:example-k.txt,if not exists k=10)

```
file content format:
     10
     20
     30
```
		
- r : Algorithm Id (default:1)

```	
6:TESTING (testing parameters)

12:FRTM
19:DFRTM
20:FRTMFORDYN (FRTM revision for DFRTM)
	
22:FRTMOPT1 (FRTM OPT1)
23:FRTMOPT1DYN (FRTM OPT1 revision for DFRTM OPT1)
24:DFRTMOPT1 (DFRTM OPT1)

31:FRTMPLUS (FRTM+)
33:FRTMPLUSDYN (FRTM+ revision for DFRTM+)
34:DFRTMPLUS (DFRTM+)
```

- o : Output level of motif (default:0)

```
three levels:

	0:only output motif number, the running time and memory use

	1:except for those outputs mentioned in 0, output the edges number of every motif

	2:except for those outputs mentioned in 1, output the detailed information of motif edges with their labels (output labels)

	3:except for those outputs mentioned in 1, output the ids of motif edges (in one line, for testing)
```

- l : Limit the start time and end time of input data (default:doesn't limit) 

```
format: -l:50, means you only use the interval [0,50] of the input temporal graph in the  the algorithm
```

- n : Limit the end time of input data when snapshot increasing (use -l at the same time) (default:doesn't limit)

```		
format: -l:500 -n:1600 means that the interval of temporal graph is [0,500] before snapshots increase,  and that the interval of temporal graph is [0,1600] after snapshots increase (used in incremental algorithm)
```

- a : The value of the global relaxation bound

```
format: -a:5,value(double)
    e.g.: -a:5,0.04 means the global relaxation bound is 4%
```

- b : The value of the local relaxation bound
	
```
format: -b:1,value(integer), type= 1:  control the maximum length of each gap for noise of all edges

e.g.: -b:1,3 means the gap of each noise duration in each edge of motifs is less than or equal to 3
```

**Examples for static algorithm FRTM**
```
FRTM.exe -i:[graph file] -k:[k file] -v:[vertex number] -e:[edge number] -t:[snapshot number] -a:0.04 -b:3 -r:[algorithm id] -f:output.txt -o:2 
(use DEL-Table, output the motif number, the running time, memory use and detailed information of motif edges with their labels, global relaxation bound = 4%, local relaxation bound = 3)
```

**Examples for incremental algorithm DFRTM**
```
FRTM.exe -i:[graph file] -k:[k file] -v:[vertex number] -e:[edge number] -t:[snapshot number] -l:500 -n:1000 -a:0.04 -b:3 -r:[algorithm id] -f:output.txt -o:0
(use DEL-Table, only output motif number, the running time and memory use, the interval of temporal graph is from [0,500] to [0,1000], global relaxation bound = 4%, local relaxation bound = 3) 
```		


## Datasets
- BJD<font size = 8>ATA</font>: A real-life dataset records road traffic conditions in Beijing. There are three traffic conditions (i.e. 2: congested, 1: slow, -1: fast), and are updated every 5 minutes. 

- EURD<font size = 8>ATA</font>: A real-life dataset for a renewable European electric power system. Each edge represents a transmission line with only static properties, and each node represents a merge point of transmission lines with a dynamic property of the hourly energy demand. We transformed the network by switching edges to nodes and nodes to edges, and used $10$ labels to discretize energy demand signals per $200$ MWh ($1$-$10$: from lowest to highest, such as `$1$' means `$0 \sim 200$ MWh').

- LAND<font size = 8>ATA</font>: A real-life dataset with computer events from the [Los Alamos National Laboratory](https://csr.lanl.gov/data/2017/) enterprise network. Each node represents a computer, and each edge represents the communication of two computers with properties of the start time, the duration, and the total number of bytes sent. We first discretized the temporal network into 5-minute snapshots. For each edge, we uniformly allocated its total bytes sent across all snapshots overlapping with its duration. Then, we used $6$ labels to discretize the number $x$ of bytes sent ($0$: the virtual label means no bytes sent, $1$-$5$: from lowest to highest, such as `$1$' means `$0< x < 4$ KB'). We selected the edges with $<50\%$ virtual labels, and transformed five-day data (Day $3$-$7$) to the network with $1440$ snapshots.

- SYND<font size = 8>ATA</font>: Synthetic datasets produced by the synthetic data generator.

If you want to get the processed data in our paper, please feel free to contact us!   