# Contents

What each file contains:

- `orbitmesy.py`: Exc and its inverse, rowmotion for 321 and 123-avoiding permutations, conversions from Dyck Paths to these types of permutations, and the other way around.
- `knuth_moves.py`: Implementations of Knuth moves K1, K2, and KB, detection of whether a tableau is good or bad, a function to generate the superstandard tableau of a given shape, and a function that plots how permutations are related via Knuth moves.
- `tableau_heights.py`: Possibly not needed anymore. Calculates how many bad tableaux of a certain size have height equal to their soliton decomposition's height.
- `superstandard_words.py`: Superstandard words for tableaux, testing some hypotheses about number of abc(d) patterns and heights of SD(w).
- `graph.png`: Sample output from knuth_moves.py.
- `Dockerfile`: Don't change this, but I think it needs to be here.
- `run_sage`: A script to make running the Docker image easier.

# Docker

We encountered some difficulties installing `graphviz` and `dot2tex` on OSes other than Linux, and running Sage on Windows at all is apparently quite painful. Docker provides "containers," which are little virtual OSes that contain nothing but the exact software needed for a specific application. It makes running software like this cross-platform much easier, hopefully.

Instructions on how to install Docker are on their website. Once you install it, make sure you have a command called `docker` accessible from the command line.

To run Sage via Docker, download and run the `run_sage` script. You may have to run `chmod +x run_sage` to ensure it's actually executable first. Running just `./run_sage` on its own will open a Sage interactive prompt. Following it with a filename of a Python script (`./run_sage orbitmesy.py`) will run that script. The Docker image has access to your local filesystem, so if you run a script that creates a new file, it will put it in the current directory.

If you're using Windows, use WSL. The `run_sage` script is written in Bash and thus won't work on the Windows command line. I have no computers with Windows installed and thus can't test anything related to that OS on my own.