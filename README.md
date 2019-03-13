# Heat Transfer Simulation

<p>
This project simulates the heat transfer equation for a specified amount of steps, as a part of a university assignment.
It is written in C, and parallelised using OpenMP to support multiple processors.
</p>

## Getting Started

### Prerequisites

<ul>
<li>gnuplot</li>
</ul>

<p>
For use with OpenMP on Mac OS X, install GCC 8 with:
</p>

```
$ brew install gcc@8
```

### Building

<p>
At the project root, build the executable with:
</p>

#### Linux

```
$ gcc -o hte -fopenmp hte.c
```

#### Mac OS X

```
$ gcc-8 -o hte -fopenmp hte.c
```

### Running

<p>
The executable can be run as shown:
</p>

```
# Run the simulation with default settings (single threaded, run for 10000 steps)
$ ./hte

# Run the simulation specifying number of threads and steps
$ ./hte --num-threads 4 --steps 2000

# Run the simulation multiple times, output a graph of computational time against number of threads
$ ./hte --investigate
```

<p>
Once completed, the output is written to a gnuplot script and can be viewed with:
</p>

```
# Run script in gnuplot
$ gnuplot heatmap.gnu
```
