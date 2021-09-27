import msprime

def twopulse_upd(T=[10,10], M=[0.2, 0.2], sample_sizes=[100, 100, 100], mu=1.25e-8, rho = 1.6e-9, N_haploid = 1000, Tdiv = 4000, lenght_m = 1, seed=1):
    length=int(lenght_m/rho)

    dem = msprime.Demography()
    dem.add_population(name='H', initial_size=N_haploid)
    dem.add_population(name='F', initial_size=N_haploid)
    dem.add_population(name='G', initial_size=N_haploid)
    dem.add_population(name='old', initial_size=N_haploid)

    dem.add_mass_migration(time=T[1], source = "H", dest = "G", proportion=M[1])
    dem.add_admixture(time=T[0]+T[1], ancestral=["F","G"], derived="H", proportions=[1-M[0], M[0]])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()

    ts = msprime.sim_ancestry(samples={"H": sample_sizes[0]/2, "F": sample_sizes[1]/2, "G": sample_sizes[2]/2}, demography=dem,
        sequence_length = length, recombination_rate=rho, model=[msprime.DiscreteTimeWrightFisher(duration=50) ,msprime.StandardCoalescent(duration=3950)], random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    return mts
