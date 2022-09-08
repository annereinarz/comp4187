---
title: "Parallel Hello World"
weight: 1
katex: false
---

# Hello, World!

Every programming course has to start with "hello world", this is no
exception. The goal of this is to familiarise you with compiling and
running code using MPI, the parallel library we'll be using, either on
Hamilton, or your own machine. So take a look at the [setup guide]({{<
ref "/setup/mpi.md" >}}) if you haven't already.

## A Python version {#mpi}

MPI is a specification for a library-based programming model. The
standard specifies Fortran and C/C++ interfaces, and there are
wrappers for many popular programming languages including
[Python](https://mpi4py.readthedocs.io/en/stable/) and
[Julia](https://github.com/JuliaParallel/MPI.jl).

We're using `mpi4py` in this course, but we'll need to install it
(either on our machines, or Hamilton).

On Hamilton, we need to load the right modules. We'll need a modern
version of Python, along with an MPI library. If you're
running on your own system and have set things up right, you don't
need modules.

To load a module (which sets up the environment) we need to run
`module load MODULENAME`. So log in and load

```
python/3.6.8
gcc/8.2.0
intelmpi/gcc/2019.6
```

Our hello world code looks like this

{{< code-include "parallel/hello/hello.py" "py" >}}

This doesn't look that different from a serial hello world, except
there are a bunch of additional calls to library functions in the
`MPI` package.

We need to install `mpi4py` to be able to run this, so let's make a
virtual environment

```
$ python -m venv comp4187
```

and activate it

```
$ . comp4187/bin/activate
(comp4187) $
```
Your prompt should have changed.

Now install `mpi4py`

```
(comp4187) $ pip install mpi4py
```

Let's confirm that it works using only one process by running

```
(comp4187) $ python hello.py
```

To run in parallel, we need to use the `mpirun` launcher. This takes
care of allocating parallel resources and setting up the processes
such that they can communicate with one another.

We can run it on the login node using four parallel processes with

```
(comp4187) $ mpirun -n 4 python hello.py
```

{{< hint warning >}}

You should not use the login node for large-scale computation, but
instead use a batch script and submit to the batch system. You should
also use the batch system when benchmarking.

{{< /hint >}}

If we want to submit to the backend, we need to write a batch script

{{< code-include "parallel/hello/hello.slurm" "sh" >}}


{{< hint "info" >}}

On some systems you need to specify the number of processes you want
to use when executing `mpirun`. However, on Hamilton, the metadata in
the scheduler is used to determine the number of processes. Here we
request 1 node and 24 tasks per node.

Hence you should only explicitly specify if you want to run with an
amount of parallelism different to that specified in your submission
script.
{{< /hint >}}

{{< question >}}

Try running on two compute nodes, by changing the `--nodes=1` line to
`--nodes=2` in the batch script. How many total processes do you now
have? What do you notice about the node names?

{{< /question >}}
