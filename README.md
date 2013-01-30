BEASTifier
====================
Generate a heap of BEAST xml files for simulation experiments

Compile
---------------

In a  terminal prompt, type:

	make

Usage
---------------

Type:

	./BEASTifier -h

for help.

Usage:

BEASTifier utilizes a configuration file for all analysis parameters. Call as:

	./BEASTifier -config config_filename

where 'config_filename' contains all analysis settings. Parameters are listed one per line, in any order. The character '#' is used for comments.

###Arguments:

   -alist: filename
      - name of text file listing alignment filenames.
      - one alignment filename per line.
      - Required; all other arguments are optional.
   -mods: list substitution model(s) to analyze data
      - supported models: JC, K80, HKY, TrNef, TrN, K3P, K3Puf, TIMef, TIM, TVMef, TVM, SYM, GTR.
      - if more than one model, separate by spaces.
      - models themselves must contain no spaces and at most one '+'.
         - e.g. 'GTR+IG' = good; 'GTR+I+G' no es bueno.
      - default: -mods JC HKY GTR JC+G HKY+G GTR+G
   -clock: list of flavour(s) of clock model to implement.
      - supported models: 'strict' or 'ucln' of 'uced'
      - default = -clock ucln
   -mcmc: the number of mcmc generations to run analysis.
      - default: -mcmc 20000000
   -tsamp: the interval (in generations) for sampling trees.
      - default: -ts 5000
   -psamp: the interval (in generations) for sampling parameter values.
      - default: -ps 1000
   -ssamp: the interval (in generations) for printing results to standard output.
      - default: -ss 500
   -logphy: turn on logging of phylograms (in addition to chronograms).
      -  default: don't log.
   -tprior: specify tree prior.
      - supported: 'bd', 'yule', 'concoal', 'expcoal', 'logcoal'.
      - default = -tprior bd
   -fixtree: turn off topology manipulation operators (i.e. fix to input topology).
      - default = estimate topology.
   -rprior: specify prior for root age.
      - supported: 'unif'' or 'norm'.
      - if 'unif', expecting '-rprior unif min_value max_value'.
      - if 'norm', expecting -rprior norm mean_value stdev_value'.
   -overwrite: overwrite existing files.
      - default = don't overwrite; warn instead.

Consult 'config.example' as a, well, example.