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

```bash
$ brew install gcc@8
```

### Building

<p>
At the project root, build the executable with:
</p>

#### Linux

```bash
$ gcc -o hte -fopenmp hte.c
```

#### Mac OS X

```bash
$ gcc-8 -o hte -fopenmp hte.c
```

### Running

<p>
The executable can be run as shown:
</p>

```bash
# Run the simulation with default settings (uses all available processors, run for 10000 steps)
$ ./hte

# Run the simulation specifying number of threads and steps (and write the output to file)
$ ./hte --num-threads 2 --steps 3000 --write-to-file

# Run the simulation multiple times while varying the number of threads (and write output to file)
$ ./hte --investigate-threading
```

<p>
Once completed, the output is written to a gnuplot script and can be viewed with:
</p>

```bash
# Run script in gnuplot to view the heatmap
$ gnuplot heatmap.gnu
```
