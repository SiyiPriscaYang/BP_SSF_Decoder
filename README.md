# BP_SSF_Decoder
This software implements the BP+Small Set Flip (SSF) decoder proposed in the following paper, with modifications tailored for the quaternary BP decoder:
```
Grospellier, A., Grouès, L., Krishna, A., & Leverrier, A. (2021). "Combining hard and soft decoders for hypergraph product codes." Quantum Journal, volume 5, page 432. [DOI:10.22331/q-2021-04-15-432](https://doi.org/10.22331/q-2021-04-15-432)  
```
For QLDPC codes with large codeword lengths, directly running BP decoders with post-processing schemes like OSD can be computationally expensive. Alternatively, using a BP+SSF decoder first to collect a small set of error vectors that result in decoding failures, followed by a more advanced BP+OSD decoder, can significantly reduce simulation time.

# Compilation of the Decoder
The file `nonbinary_decoder_simulation_degenerate2.cpp` is the main file. Run the following code at the terminal to compile it and obtain an executable `qldpc_nb_deg2_th`:
```
g++ nonbinary_decoder_simulation_degenerate2.cpp -o qldpc_nb_deg2_th
```

# Example of a SLURM job submission script 
Below is an example of a SLURM job submission script. 
```
#!/bin/bash

#SBATCH -p common
#SBATCH --mem=4G
#SBATCH -o errors_%a.out
#SBATCH -a 1-200

./qldpc_nb_deg2_th Code2_3_7.txt 0.08 10 500 0.0 $SLURM_ARRAY_TASK_ID
```

## Instructions for Customizing the Script:

- Parity Check Matrix:

  Replace `Code2_3_7.txt` with the path to your .txt file containing the parity check matrix for your QLDPC code. The matrix should be formatted as follows:
  
  - Use 0, 1, 2, 3 to represent I, X, Y, and Z operations, respectively.
  - Each row in the file should correspond to a stabilizer, with entries separated by spaces.

- Depolarization Error (p):

  Replace 0.08 with the desired depolarization error probability (p). The error distribution for the operations (I, X, Y, Z) is defined as:

  - (I, X, Y, Z) = (1 - p, p/3, p/3, p/3)

- Number of Trials:

  Replace 10 with the number of trials to run for each job. This value controls how many times the simulation will be repeated within each individual job.

- Job Array Indices:

  Modify `#SBATCH -a 1-200` to set the range of job indices to be submitted. For example, if you want to run 300 jobs, change this to `#SBATCH -a 1-300` or `#SBATCH -a 201-500`.

- Maximum BP Iterations:

  Replace 500 with the upper limit for the number of belief propagation (BP) iterations. This defines the maximum number of iterations to attempt during the decoding process.

## Output Files

The output files `job_k.txt`, `xv_k.txt`, and `ev_k.txt` (where `k` corresponds to the job index) store the following data:
- `job_k.txt`: Error statistics
- `xv_k.txt`: Input error vectors
- `ev_k.txt`: Output error vectors (the difference between the estimated and input error vectors)

## Gathering Errors and Obtaining Statistics

The file `error_analysis.cpp` is used to combine the error statistics, input error vectors, and output error vectors collected from different jobs that run the decoder. To compile it, use the following command:

```bash
g++ error_analysis.cpp -o error_analysis
```
Run the command below to generate the error statistics:
```
./error_analysis Code2_3_7_0.08 1 200
```
# Example Terminal Output:
```
Average FER: 0.18 // this is the LER we use to plot the simulation result
total errors: 360
total convergent errors: 12
total non-convergent errors: 348
total trials: 2000
Average BER: 0.000694914
Total number of erroneous bits: 8061
Number of columns: 5800
```
Additionally, two .txt files will be generated:

- `Code2_3_7_0.08_input.txt`: Contains the input error vectors corresponding to decoding failures. You could use this file as input of more advanced decoders like BP+OSD decoder to save simulation time.
- `Code2_3_7_0.08_error.txt`: Contains the output error vectors corresponding to decoding failures. You could use this file to analyze the error profile.



