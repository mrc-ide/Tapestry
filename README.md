 <!-- badges: start -->
  [![checks_develop](https://github.com/mrc-ide/Tapestry/workflows/checks_develop/badge.svg)](https://github.com/mrc-ide/Tapestry/actions)
  [![checks_master](https://github.com/mrc-ide/Tapestry/workflows/checks_master/badge.svg)](https://github.com/mrc-ide/Tapestry/actions)
  <!-- badges: end -->

# Tapestry

Malaria infections often contain multiple genotypes, and when sequenced together
these produce a complex signal that is a mixture of the individual genotypes.
Building on the framework of earlier programs like DEploid and DEploidIBD,
Tapestry attempts to pull these individual genotypes apart by exploiting allele
frequency imbalances within a sample, while simultaneously estimating segments
of identity by descent (IBD) between sequences. Unlike previous programs,
Tapestry uses advanced MCMC methods to ensure that results are robust even for
high complexity of infection (COI).

## Install
### Build from source

**Step 1:** Install [HTSlib](https://github.com/samtools/htslib) dependency.
e.g.
```
git clone https://github.com/samtools/htslib
cd htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install
```

**Step 2:** Clone Tapestry.

```
git clone https://github.com/mrc-ide/Tapestry.git
cd Tapestry
```
**Step 3:** Compile Tapestry using [CMake](https://cmake.org/).

For a (slower) debugging version:

```
cd build
cmake ..
make .
```

For a (faster) release version:

```
mkdir release
cd release
cmake .. -DCMAKE_BUILD_TYPE="Release"
make .
```

**Optional**: By default the tests are not compiled. To compile them, open [CMakeLists.txt](https://github.com/mrc-ide/Tapestry/blob/feature/likelihood/CMakeLists.txt) in a text editor, and change line 16:

```
set(COMPILE_TESTS OFF)
```

...to...

```
set(COMPILE_TESTS ON)
```

Repeat Step 3 and then run tests with the executable `./test_tapestry`. The testing framework is [GoogleTest](https://github.com/google/googletest).

## Quickstart
The executable `tapestry` will be in your `/build` or `/release`.

```
bash_prompt$ ./release/tapestry infer --help
Run inference from an filtered VCF.
Usage: ./release/tapestry infer [OPTIONS]

Options:
  -h,--help                   Print this help message and exit


Input and output:
  -i,--input_vcf TEXT:FILE REQUIRED
                              Path to input VCF file.
  -s,--target_sample TEXT REQUIRED
                              Target sample in VCF.
  -o,--output_dir TEXT        Output directory.


Model Hyperparameters:
  -K,--COI INT:INT in [1 - 6] Complexity of infection.
  -e,--error_ref FLOAT:INT in [0 - 1]
                              Probability of REF->ALT error.
  -E,--error_alt FLOAT:INT in [0 - 1]
                              Probability of ALT->REF error.
  -v,--var_wsaf FLOAT:POSITIVE
                              Controls dispersion in WSAF. Larger is less dispersed.
  -r,--recomb_rate FLOAT:POSITIVE
                              Recombination rate in kbp/cM.
  -b,--n_wsaf_bins INT:INT in [100 - 10000]
                              Number of WSAF bins in Betabin lookup table.


MCMC Parameters:
  -w,--w_proposal FLOAT:POSITIVE
                              Controls variance in proportion proposals.
```

## Helper scripts
[scripts](https://github.com/mrc-ide/Tapestry/tree/feature/likelihood/scripts) contains python scripts for:

1. Quickly running Tapestry over a set of simulated *P. falciparum* infections
2. Visualising outputs

You can get help with `python scripts/<script_name> --help`:

```
bash_prompt$ python scripts/infer_multiple_samples.py --help
Usage: infer_multiple_samples.py [OPTIONS]

  Run `Tapestry` over a selection of samples

Options:
  -i, --input_vcf PATH     Path to input VCF.
  -c, --summary_csv PATH   Path to summary CSV.
  -n, --n_samples INTEGER  Number of samples.
  -K, --coi INTEGER        COI to run.
  -w, --w_proposal FLOAT   Proportion-titre proposal SD ~N(0, w).
  -v, --var_wsaf FLOAT     Variance in WSAF.
  --help                   Show this message and exit.

```

Here is an example workflow:

```
input_vcf="example_data/simulated_infections.DRCongo.K03.vcf"
summary_csv="simulated_infections.DRCongo.K03.summary.csv"

python scripts/infer_multiple_samples.py -i $input_vcf -c $summary_csv
```

Outputs will be deposited in a directory `/results`. Visualise with:

```
python scripts/plot_mcmc_diagnostics -c $summary_csv
```

<p align="center"><img src="images/plot.parameters.SMI002-K03-B03-M00.png" width="1000"></p>

## Model

We imagine a scenario where our sequencing data is generated by a fixed number of genetically distinct strains, which may or may not share regions of identity by descent. Each strain is imagined to comprise a fixed proportion of the infection. The fraction of reads that are derived from each strain is influenced by this proportion.

### Likelihood

The likelihood is formulated as a Hidden Markov Model:

$$ P(\vec{X}, \vec{S} | \Theta) = P(S_1) P(X_1 | S_1) \prod_{i=1}^{L} P(S_i|S_{i-1}) P(X_i | S_i) $$

As such, we can compute it by defining initiation, transition, and emission probabilities.

#### Emission probabilities, $P(X_i|S_i)$

Given a set of proportions and haplotype states, we can compute the proportion of the sample comprised of the alternative allele (or rather the expected WSAF without accounting for error), as:

$$ q_i = \sum_{j=1}^{K} w_jh_{ij} $$

We then adjust for sequencing error by assuming two fixed error rates, $e_0$ and $e_1$:

$$ \pi_i = q_i(1-e_1) + (1 - q_i)e_0 $$

If we sequenced every parasite genome in the host, our resultant observed WSAF, $x_i$, would be exactly $\pi$. However, in practice we generate a finite number of reads sampled from the infection through a random process. The sampling variation is modelled using a [Beta-binomial distribution](https://en.wikipedia.org/wiki/Beta-binomial_distribution):

$$ X_i \sim BetaBin(a_i + r_i, \alpha, \beta) $$

We re-parameterise the distribution to have better control over its mean and variance:

$$\alpha=v\pi$$ 

$$\beta=v(1-\pi)$$

With this parameterisation, the expected error-adjusted WSAF is:

$$ E[X_i|v, \pi_i] = \frac{\alpha}{\alpha+\beta} = \pi_i $$

as we sought. Additionally we have:

$$Var(X_i|v, \pi) \propto \frac{1}{v}$$

which gives us good control over the variance.

If reads were sampled independently, randomly, and with replacement from the underlying strain proportions, we would expect the observed WSAF to be binomially distributed. However, we would like to allow additional dispersion across SNPs. For example, genomic context can influence sequencing performance, and that context will be different for each SNP; SNPs with identical haplotype configurations should have more than just binomial variance in their observed WSAF. The $v$ term in the beta-binomial allows us to capture this additional dispersion.

Rather than explicitly inferring haplotypes at each site, we marginalise over all possible haplotype configurations by making our final emission probability a finite mixture of beta-binomial distributions. Each haplotype configuration, $\vec{h}$, corresponds to a particular subset of the $K$ strains carrying the alternative allele. Note that $|\vec{h}| = 2^K$, and we index these configurations with a superscript $b$. For example, if $K=2$, the $\vec{h}$ would be $\Set{0,0}, \Set{0, 1}, \Set{1, 0}, \Set{1, 1}$. In essence, each $b$ generates a mode in the multi-modal WSAF distribution.

Our emission probability at each site becomes:

$$ P(X_i|S_i, w_j, p_i, v, e_0, e_1) = \sum_{b=1}^{2^{K}} P(a_i, r_i | \vec{h^{b}}, S_i, w_j, v, e_0, e_1) P(\vec{h^{b}}|S_i, p)$$

The IBD configuration $S_i$ limits which haplotype configurations are possible, by restricting strains in IBD to have the same haplotype state. In addition, we assume each group of strains in IBD is sampled independently and randomly. Defining the number of IBD groups carrying the alternative allele as $C_{h=1|S_i}$. Then,

$$ P(\vec{h^{b}}|S_i, p)=p^{C_{h=1|S_i}}(1-p)^{K-C_{h=1|S_i}}$$

#### Initiation probabilities, $P(S_1)$

TODO

#### Transition probabilities, $P(S_i|S_{i-1})$

TODO




### Parameters
#### Data
| Parameter | Description |
| --------- | ----------- |
| $L$ | Number of SNPs. |
| $i$ | SNP index, $i \in \Set{1, ..., L}$. |
| $r_{i}$ | Reference (REF) allele count. |
| $a_{i}$ | Alternative  (ALT) allele count. |
| $x_{i} := \frac{a_{i}}{r_{i}+a_{i}} $ | Observed within-sample alternative allele frequency (WSAF). |
| $p_{i}$ | Estimated population-level alternative allele frequency. |
| $d_{i,i+1}$ | Physical distance between SNPs, in basepairs. |

#### Model
| Parameter | Description |
| --------- | ----------- |
| $K$ | Number of strains in sample, i.e. complexity of infection (COI). |
| $j$ | Strain index, $j \in \Set{1, ..., K}$. |
| $w_j$ | Abundance of strain $j$ as a fraction; proportion of sample comprised of strain $j$. |
| $h_{ij}$ | Haplotype state with $h_{ij} \in \Set{0, 1}$ for REF and ALT, respectively.|
| $q$ | The expected WSAF, without error adjustment. |
| $\pi$ | The error-adjusted expected WSAF. |
| $s_i$ | Index for the IBD state, $S_i \in \Set{1, ..., B_K}$, where $B_K$ is the [Bell number](https://en.wikipedia.org/wiki/Bell_number). |

#### Hyper-parameters
| Parameter | Description | Value | Reference |
| --------- | ----------- | ----- | --------- |
| $\rho$ | Recombination rate. | 13.5 kbp per centiMorgan | [Miles et al. (2016)](https://genome.cshlp.org/content/26/9/1288.full) |
| $e_0$ | REF to ALT read count error rate. | 0.01 | Calibrated from Pf3k |
| $e_1$ | ALT to REF read count error rate. | 0.05 | Calibrated from Pf3k |
| $v$ | Term setting dispersion in WSAF. | 500 | Calibrated from Pf3k |








