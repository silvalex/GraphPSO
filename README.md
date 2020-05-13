# GraphPSO

This project provides the implementation for the indirect graph-based particle swarm optimisation approach presented in the paper entitled "Particle Swarm Optimisation with Sequence-Like Indirect Representation for Web Service Composition", which is available [here](https://link.springer.com/chapter/10.1007/978-3-319-30698-8_14). The core idea of this composition approach is to optimise the order of a queue of services, which is then decoded into the corresponding Web service composition graph for fitness calculation. The creation of a separate decoding mechanism allows candidates to be optimised with fewer constraints, since the functional correctness of solutions is ensured during the construction of their corresponding graphs.

## Running the Code

* Clone the repository:
```
git clone git@github.com:silvalex/GraphPSO.git
```

* The code can then be compiled and run, either from an IDE or directly from the command line. The following is an example of how the code can be compiled then run via command line. First we create a bin directory and then compile the java files (the class files are outputted to bin):
```
mkdir bin
javac -d bin src/pso/*
```

* We can then run the code using the example dataset in the repository. This is done as follows:
```
java -classpath bin pso.GraphPSO log.txt histogram.txt testTaskSet.xml testDataset.xml testTaxonomySet.xml 1
```
In the command above, we add the class files we compile to the classpath (`-classpath bin`), then we run from the application's entry point (`pso.GraphPSO`) providing the following arguments: 1) Name of log file, 2) Name of histogram file, 3) Path to dataset file containing composition task, 4) Path to dataset file containing service descriptions, 5) Path to dataset file containing taxonomy information, 6) Seed for random number generator.
> **_Note:_** When running the code from Eclipse, the arguments can be passed by creating a Run Configuration for the code (using `pso.GraphPSO` as the main class), then listing the six arguments in the Arguments tab.

## Running the Code using WSC Datasets

In the paper, experiments are conducted using the WSC-2008 and WSC-2009 datasets, which contains a variety of composition tasks with varying levels of complexity. These datasets can be found in ZIP format [in this repository](https://github.com/silvalex/WSC2008_2009). Once unzipped, WSC-2008 contains 8 directories with different composition tasks and WSC-2009 contains five. Each of these directories contains a composition task (problem.xml), a service description file (services-output.xml), and a taxonomy information file (taxonomy.xml). In order to run the code using one of these tasks, simply  change the paths to composition task, service description, and taxonomy information files. For instance, assuming that WSC-2008 has been unzipped to the directory `˜/git/wsc2008`, the first composition task can be run as follows:
```
java -classpath bin pso.GraphPSO log.txt histogram.txt ˜/git/wsc2008/Set01MetaData/problem.xml ˜/git/wsc2008/Set01MetaData/services-output.xml ˜/git/wsc2008/Set01MetaData/taxonomy.xml 1
```

## Interpreting the Output

Once the code is successfully run, it will generate two output files: A log file (in our example above, `log.txt`) and a histogram file (in our example, `histogram.txt`).

### Log File

In the log file, each line contains information on a single iteration of the PSO algorithm, and the last line of the log file contains a DOT representation of the graph corresponding to the best global solution identified throughout the run. Each line has 18 columns, with each column value corresponding to the following element:

1. Iteration number for PSO (zero-based).
2. Initialisation time in milliseconds before this iteration (in practice, this value is only greater than 0 in the first line of the file, when the particle swarm is being created).
3. Time in milliseconds elapsed in this iteration.
4. Mean fitness across all particles in this swarm in this iteration.
5. Highest fitness in the swarm in this iteration.
6. Highest fitness identified by the swarm so far (including previous iterations).
7. Mean availability value (not normalised) across all corresponding graph solutions in this iteration.
8. Mean reliability value (not normalised) across all corresponding graph solutions in this iteration.
9. Mean time value (not normalised) across all corresponding graph solutions in this iteration.
10. Mean cost value (not normalised) across all corresponding graph solutions in this iteration.
11. Highest availability value (not normalised) in the graph solutions in this iteration.
12. Highest reliability value (not normalised) in the graph solutions in this iteration.
13. Highest time value (not normalised) in the graph solutions in this iteration.
14. Highest cost value (not normalised) in the graph solutions in this iteration.
15. Highest availability value (not normalised) identified in the solutions so far (including previous iterations).
16. Highest reliability value (not normalised) identified in the solutions so far (including previous iterations).
17. Highest time value (not normalised) identified in the solutions so far (including previous iterations).
18. Highest cost value (not normalised) identified in the solutions so far (including previous iterations).

The log file is structured in this way so that it can be easily statistically analysed and/or visualised using a language such as R.

> **_NOTE:_** The DOT representation in the last line of the log file can be converted to a PDF/PNG file, which contains a rendering of the graph. This is done by saving that line as a `.dot` file and then using the utilities provided by GraphViz (please see [their guide](http://www.graphviz.org/pdf/dotguide.pdf) for further information).

### Histogram File

The histogram contains four lines, defined as follows:

1. The names of all services that appear at least once in graph solutions generated throughout the PSO run (note that each name is a single token).
2. The number of times each service in the line above has appeared in the graph solutions generated throughout the PSO run.
3. The names of all edges that appear at least once in graph solutions generated throughout the PSO run (each edge name is single-token combination of the origin node name and the destination node name, connected by "->").
4. The number of times each edge in the line above has appeared in the graph solutions generated throughout the PSO run.

## Changing PSO Parameter Settings

Parameters for the algorithm are set as constants in the GraphPSO class, and these can be modified when new settings are desired:
* `MAX_NUM_ITERATIONS`: Controls the number of iterations in the search (default: 100).
* `NUM_PARTICLES`: Controls the size of the swarm (default: 30).
* `C1`: Controls the first acceleration coefficient in PSO (default: 1.49618).
* `C2`: Controls the second acceleration coefficient in PSO (default: 1.49618).
* `W`: Controls the inertia coefficient in PSO (default: 0.7298).
* `dynamicNormalisation`: Boolean that specifies whether normalisation boundaries for QoS attributes should be recalculated at each iteration. If true, the highest/lowest QoS values across all graphs in one iteration are used as the upper/lower limits for normalisation. If false, the upper/lower limits are estimated once, at the beginning of the run (default: true).
* `W1`: Used in fitness function to control the weight of the availability component (default: 0.25).
* `W2`: Used in fitness function to control the weight of the reliability component (default: 0.25).
* `W3`: Used in fitness function to control the weight of the time component (default: 0.25).
* `W4`: Used in fitness function to control the weight of the cost component (default: 0.25).

> **_NOTE:_** `W1`, `W2`, `W3`, and `W4` should always add up to 1.

## Running the code on the ECS Computing Grid

This code was developed at the School of Engineering and Computer Science (ECS) at Victoria University of Wellington, New Zealand. ECS has a computing grid that allows for parallel execution of many program instances, each of which is run by submitting a bash script. Details on how to use the ECS grid can be seen [here](https://ecs.wgtn.ac.nz/Support/TechNoteEcsGrid). Three files in this repository were written for use with the grid:

* `graph_pso.sh`: This is the bash script that defines the run of a single instance of the algorithm. It copies the code and dataset to the computer selected by the grid, for running locally. It then copies the log and histogram outputs to a directory for collection. Note that the paths throughout this script must be modified to match your own directory tree.
* `run_pso_grid_2008.sh`: This script repeatedly submits `graph_pso.sh` to the grid, testing all the different tasks in WSC-2008.
* `run_pso_grid_2009.sh`: This script repeatedly submits `graph_pso.sh` to the grid, testing all the different tasks in WSC-2009.
