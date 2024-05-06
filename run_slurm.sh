## ER stats
# sbatch ./runner.slurm netres.py -g er -n 2 -s diam -r rand -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g er -n 2 -s diam -r targ -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g er -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g er -s frag -r targ -p --nperms 5 --seed 1000

## PA stats
# sbatch ./runner.slurm netres.py -g pa -n 2 -s diam -r rand -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g pa -n 2 -s diam -r targ -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g pa -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g pa -s frag -r targ -p --nperms 5 --seed 1000

## ER2 stats
# sbatch ./runner.slurm netres.py -g er2 -n 2 -s diam -r rand -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g er2 -n 2 -s diam -r targ -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g er2 -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g er2 -s frag -r targ -p --nperms 5 --seed 1000

## PA2 stats
# sbatch ./runner.slurm netres.py -g pa2 -n 2 -s diam -r rand -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g pa2 -n 2 -s diam -r targ -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g pa2 -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g pa2 -s frag -r targ -p --nperms 5 --seed 1000

## ERPA stats
# sbatch ./runner.slurm netres.py -g erpa -n 2 -s diam -r rand -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g erpa -n 2 -s diam -r targ -p --nperms 2 --seed 1000
# sbatch ./runner.slurm netres.py -g erpa -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g erpa -s frag -r targ -p --nperms 5 --seed 1000

## internet stats
# sbatch ./runner.slurm netres.py -g int -n 1 -s diam -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g int -n 1 -s diam -r targ -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g int -n 1 -s frag -r rand -p --nperms 5 --seed 1000
# sbatch ./runner.slurm netres.py -g int -n 1 -s frag -r targ -p --nperms 5 --seed 1000

## world wide web stats
# sbatch ./runner.slurm netres.py -g www -n 1 -s diam -r rand -p --nperms 1 --seed 1000
# sbatch ./runner.slurm netres.py -g www -n 1 -s diam -r targ -p --nperms 1 --seed 1000
# sbatch ./runner.slurm netres.py -g www -n 1 -s frag -r rand -p --nperms 1 --seed 1000
# sbatch ./runner.slurm netres.py -g www -n 1 -s frag -r targ -p --nperms 1 --seed 1000
