## Prep
# create env -name tcrdist3
# pip install tcrdist3==0.2.2

## Way to run
# python3 tcrdist3.py
  
## subpackages
#palmotif for CDR3 logo generation
#tcrsampler for generating V-J gene-matched background receptor sets (non productive from human/mouse)
#pwseqdist for efficient and parallelizable pairwise distance computation

from tcrsampler.setup_db import install_all_next_gen
install_all_next_gen(dry_run = False)


  
