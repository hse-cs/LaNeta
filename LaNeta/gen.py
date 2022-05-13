import msprime

def twopulse_upd(T=[10,10], M=[0.2, 0.2], sample_sizes=[100, 100, 100], mu=1.25e-8, rho = 1.6e-9, N_haploid = [1000, 1000, 1000, 1000], Tdiv = 4000, lenght_m = 1, seed=1):
    length=int(lenght_m/rho)

    dem = msprime.Demography()
    dem.add_population(name='H', initial_size=N_haploid[0])
    dem.add_population(name='F', initial_size=N_haploid[1])#, default_sampling_time=T[0]+T[1])
    dem.add_population(name='G', initial_size=N_haploid[2])
    dem.add_population(name='old', initial_size=N_haploid[3])

    dem.add_mass_migration(time=T[1] + 1, source = "H", dest = "G", proportion=M[1])
    dem.add_admixture(time=T[0]+T[1] + 1, ancestral=["F","G"], derived="H", proportions=[1-M[0], M[0]])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()

    ts = msprime.sim_ancestry(
        samples={"H": sample_sizes[0], "F": sample_sizes[1], "G": sample_sizes[2]},
        demography=dem, sequence_length = length, recombination_rate=rho, ploidy=2,
        model=[msprime.DiscreteTimeWrightFisher(duration=50), msprime.StandardCoalescent(duration=3950)],
        #model='hudson',
        random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    return mts

def continuous_adm(g_s=10, m=0.5, sample_sizes=[100, 100, 100], mu=1.25e-8, rho = 1.6e-9, N_haploid = [1000, 1000, 1000, 1000], Tdiv = 4000, lenght_m = 1, seed=1):
    length=int(lenght_m/rho)

    ms = 1-(1-m)**(1/(g_s-5))

    dem = msprime.Demography()
    dem.add_population(name='H', initial_size=N_haploid[0])
    dem.add_population(name='F', initial_size=N_haploid[1])#, default_sampling_time=T[0]+T[1])
    dem.add_population(name='G', initial_size=N_haploid[2])
    dem.add_population(name='old', initial_size=N_haploid[3])

    dem.add_migration_rate_change(time=6, rate=ms, source='H', dest='G')
    dem.add_migration_rate_change(time=g_s, rate=0, source='H', dest='G')
    dem.add_admixture(time=g_s, ancestral=["F","G"], derived="H", proportions=[1-ms, ms])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()
    dem


    ts = msprime.sim_ancestry(
        samples={"H": sample_sizes[0], "F": sample_sizes[1], "G": sample_sizes[2]},
        demography=dem, sequence_length = length, recombination_rate=rho, ploidy=2,
        model=#[msprime.DiscreteTimeWrightFisher(duration=50), msprime.StandardCoalescent(duration=3950)],
        msprime.DiscreteTimeWrightFisher(),
        #model='hudson',
        random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    return mts



def twopulse_sub_upd(T=[10,10], M=[0.2, 0.2], sample_sizes=[100, 100, 100], mu=1.25e-8, rho = 1.6e-9, N_haploid = [1000, 1000, 1000, 1000], Tdiv = 4000, lenght_m = 1, seed=1):
    length=int(lenght_m/rho)

    dem = msprime.Demography.island_model([N_haploid[0]/10]*10, 0)
    dem.add_population(name='H', initial_size=N_haploid[0])
    dem.add_population(name='F', initial_size=N_haploid[1])#, default_sampling_time=T[0]+T[1])
    dem.add_population(name='G', initial_size=N_haploid[2])
    dem.add_population(name='old', initial_size=N_haploid[3])


    dem.add_population_split(time=T[1], derived=['pop_0',
                                                 'pop_1',
                                                 'pop_2',
                                                 'pop_3',
                                                 'pop_4',
                                                 'pop_5',
                                                 'pop_6',
                                                 'pop_7',
                                                 'pop_8',
                                                 'pop_9'
                                                 ],
                                                 ancestral="H")

    dem.add_mass_migration(time=T[1], source = "H", dest = "G", proportion=M[1])
    dem.add_admixture(time=T[0]+T[1], ancestral=["F","G"], derived="H", proportions=[1-M[0], M[0]])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()

    print(dem.debug())


    ts = msprime.sim_ancestry(
        samples={'pop_0':sample_sizes[0]/10,
                 'pop_1':sample_sizes[0]/10,
                 'pop_2':sample_sizes[0]/10,
                 'pop_3':sample_sizes[0]/10,
                 'pop_4':sample_sizes[0]/10,
                 'pop_5':sample_sizes[0]/10,
                 'pop_6':sample_sizes[0]/10,
                 'pop_7':sample_sizes[0]/10,
                 'pop_8':sample_sizes[0]/10,
                 'pop_9':sample_sizes[0]/10,
                 "F": sample_sizes[1],
                 "G": sample_sizes[2]},
        demography=dem, sequence_length = length, recombination_rate=rho, ploidy=2,
        model=[msprime.DiscreteTimeWrightFisher(duration=50), msprime.StandardCoalescent(duration=3950)],
        #model='hudson',
        random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    return mts
