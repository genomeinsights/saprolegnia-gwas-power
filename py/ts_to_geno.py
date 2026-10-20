"""Recapitate a clonal_metapop.slim tree sequence, add neutral mutations and
export diploid genotypes of the sampled isolates.

Recapitation is done with msprime directly (equivalent to pyslim.recapitate):
the K SLiM demes merge, at the first SLiM tick, into one sexual ancestral
population of size N*K.

Outputs (prefix = OUT):
  OUT.geno.u8   uint8 matrix, isolates x sites, row-major (0/1/2 derived dosage)
  OUT.sites.tsv chr, pos (1-based within chromosome), global position
  OUT.inds.tsv  row order of isolates: pedigree_id, deme, clone_id
"""
import argparse, csv, json
import numpy as np
import msprime, tskit

ap = argparse.ArgumentParser()
ap.add_argument("--trees", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--mu", type=float, default=3e-8)
ap.add_argument("--r", type=float, default=1e-8)
ap.add_argument("--nchr", type=int, default=4)
ap.add_argument("--chrlen", type=int, default=15_750_000)
ap.add_argument("--n", type=int, default=1000, help="SLiM deme size")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

ts = tskit.load(a.trees)
L = a.nchr * a.chrlen
assert int(ts.sequence_length) == L, (ts.sequence_length, L)
tick = ts.metadata["SLiM"]["tick"]

# --- recapitate --------------------------------------------------------------
demog = msprime.Demography.from_tree_sequence(ts)
slim_pops = []
for pop in demog.populations:
    pop.initial_size = a.n
    if ts.population(pop.id).metadata is not None and ts.population(pop.id).metadata.get("name"):
        slim_pops.append(pop.name)
demog.add_population(name="ancestral", initial_size=a.n * len(slim_pops))
demog.add_population_split(time=tick, derived=slim_pops, ancestral="ancestral")
pos = [0]
rates = []
for c in range(1, a.nchr + 1):
    pos.append(c * a.chrlen - 1); rates.append(a.r)
    if c < a.nchr:
        pos.append(c * a.chrlen); rates.append(0.5)
pos[-1] = L
rmap = msprime.RateMap(position=pos, rate=rates)
rts = msprime.sim_ancestry(initial_state=ts, demography=demog, recombination_rate=rmap,
                           random_seed=a.seed)

# --- keep sampled isolates only, then mutate -------------------------------------
want = {}
with open(a.samples) as fh:
    for row in csv.DictReader(fh):
        want[int(row["pedigree_id"])] = row
rows, nodes = [], []
for ind in rts.individuals():
    pid = ind.metadata["pedigree_id"]
    if pid in want and len(ind.nodes) == 2 and rts.node(ind.nodes[0]).time == 0:
        rows.append(want[pid]); nodes.extend(ind.nodes)
assert len(rows) == len(want), (len(rows), len(want))
sts = rts.simplify(samples=nodes, keep_input_roots=False)
mts = msprime.sim_mutations(sts, rate=a.mu, model=msprime.BinaryMutationModel(),
                            random_seed=a.seed, keep=False)

# --- genotypes ---------------------------------------------------------------------
G = mts.genotype_matrix()                      # sites x haplotypes (0/1)
D = (G[:, 0::2] + G[:, 1::2]).astype(np.uint8)  # sites x isolates
keep = (D.sum(1) > 0) & (D.sum(1) < 2 * D.shape[1])
D = D[keep]
gpos = mts.tables.sites.position[keep].astype(np.int64)
chrom = gpos // a.chrlen + 1
D.T.copy(order="C").tofile(a.out + ".geno.u8")
with open(a.out + ".sites.tsv", "w") as fh:
    fh.write("chr\tpos\tgpos\n")
    for c, p in zip(chrom, gpos):
        fh.write(f"{c}\t{p - (c - 1) * a.chrlen + 1}\t{p}\n")
with open(a.out + ".inds.tsv", "w") as fh:
    fh.write("pedigree_id\tdeme\tclone_id\n")
    for r in rows:
        fh.write(f"{r['pedigree_id']}\t{r['deme']}\t{r['clone_id']}\n")
print(json.dumps({"isolates": D.shape[1], "sites": int(D.shape[0]),
                  "mean_het": float((D == 1).mean())}))
