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
In the command above, we add the class files we compile to the classpath (-classpath bin), then we run from the application's entry point (pso.GraphPSO) providing the following arguments: 1) Name of log file, 2) Name of histogram file, 3) Path to dataset file containing composition task, 4) Path to dataset file containing service descriptions, 5) Path to dataset file containing taxonomy information, 6) Seed for random number generator. Note: When running the code from Eclipse, the arguments can be passed by creating a Run Configuration for the code (using pso.GraphPSO as the main class), then listing the six arguments in the Arguments tab.

## Running the Code using WSC Datasets

In the paper, experiments are conducted using the WSC-2008 and WSC-2009 datasets, which contains a variety of composition tasks with varying levels of complexity. These datasets can be found in ZIP format [in this repository](https://github.com/silvalex/WSC2008_2009). Once unzipped, WSC-2008 contains 8 directories with different composition tasks and WSC-2009 contains five. Each of these directories contains a composition task (problem.xml), a service description file (services-output.xml), and a taxonomy information file (taxonomy.xml). In order to run the code using one of these tasks, simply  change the paths to composition task, service description, and taxonomy information files. For instance, assuming that WSC-2008 has been unzipped in the directory ˜/git/wsc2008, the first composition task can be run as follows:
```
java -classpath bin pso.GraphPSO log.txt histogram.txt ˜/git/wsc2008/Set01MetaData/problem.xml ˜/git/wsc2008/Set01MetaData/services-output.xml ˜/git/wsc2008/Set01MetaData/taxonomy.xml 1
```

## Interpreting the Output
TBC

