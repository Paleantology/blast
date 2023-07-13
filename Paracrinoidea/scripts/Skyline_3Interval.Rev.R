
taxa <- readTaxonData("data/ages.tsv")
taxa
morpho <- readDiscreteCharacterData("data/matrix_only.nex")
morpho.addMissingTaxa(taxa)
moves = VectorMoves()
monitors = VectorMonitors()

     n_taxa <- taxa.size()

     num_branches <- 2 * n_taxa - 2



      timeline <- v(15, 32);
      # Diversification Rates based on Echinodermata
      speciation_rate ~ dnExponential(1.471);
      # NOTE: If it gets stuck in this script, then set origination & extinction to 1.0
      moves.append(mvScale(speciation_rate, lambda=0.01, weight=5));
      moves.append(mvScale(speciation_rate, lambda=0.10, weight=3));
      moves.append(mvScale(speciation_rate, lambda=1.00, weight=1));

      turnover ~ dnUnif(0.9, 1.05);
      moves.append(mvSlide(turnover, delta=0.01, weight=5));
      moves.append(mvSlide(turnover, delta=0.10, weight=3));
      moves.append(mvSlide(turnover, delta=1.00, weight=1));
      extinction_rate := turnover*speciation_rate;
      diversification := speciation_rate - extinction_rate;

      # old extinction stuff. We should not use this, as extinction should not be independent of origination!
      #extinction_rate ~ dnExponential(1.471);
      #moves.append(mvScale(extinction_rate, lambda=0.01, weight=5));
      #moves.append(mvScale(extinction_rate, lambda=0.10, weight=3));
      #moves.append(mvScale(extinction_rate, lambda=1.00, weight=1));
      #turnover := extinction_rate/speciation_rate;

      # Fossil Sampling Rates based on collection occupied by Echinodermata
      psi ~ dnExponential(3.892);
      completeness := psi/(extinction_rate+psi);
      moves.append(mvScale(psi, lambda=0.01, weight=5));
      moves.append(mvScale(psi, lambda=0.10, weight=3));
      moves.append(mvScale(psi, lambda=1.00, weight=1));

      # Proportional Taxon Sampling of Youngest Time Slice
      rho <- 0;	# 'extant' sampling.

     # Establish Basal Divergence Time
     origin_time ~ dnUnif(35, 50);
     moves.append(mvSlide(origin_time, delta=0.01, weight=5));
     moves.append(mvSlide(origin_time, delta=0.10, weight=3));
     moves.append(mvSlide(origin_time, delta=1.00, weight=1));


     fbd_dist = dnFBDP(originAge=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, timeline=timeline, taxa=taxa, condition="sampling")



  #   outgroup = clade("Cheirocystis_fultonensis");

 #constraints = v(outgroup)


     fbd_tree ~ dnConstrainedTopology(fbd_dist)


     moves.append(mvFNPR(fbd_tree, weight=15.0))

     moves.append(mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6.0))

     moves.append(mvNodeTimeSlideUniform(fbd_tree, weight=40.0))

     moves.append(mvRootTimeSlideUniform(fbd_tree, origin_time, weight=5.0))



 # Setup the fossil tip sampling #

 # Use a for loop to create a uniform distribution on the occurence time for each fossil #

 # The boundaries of the uniform distribution are specified in the tsv file #

 fossils = fbd_tree.getFossils()
 for(i in 1:fossils.size())
{
    fossils[i]
    t[i] := tmrca(fbd_tree, clade(fossils[i]))
    print("t")
    t[i]
    a_i = fossils[i].getMinAge()
    print("a")
    a_i
    b_i = fossils[i].getMaxAge()
    print("b")
    b_i
    F[i] ~ dnUniform(t[i] - b_i, t[i] - a_i)
        print("F")
    F[i]
 #   F[i].clamp( 0 )
}

 # Add a move to sample the fossil times #
 moves.append( mvFossilTimeSlideUniform(fbd_tree, origin_time, weight=5.0) )


     num_samp_anc := fbd_tree.numSampledAncestors()


     clock_morpho ~ dnExponential(1.0)

     moves.append( mvScale(clock_morpho, lambda=0.01, weight=4.0) )
     moves.append( mvScale(clock_morpho, lambda=0.1,  weight=4.0) )
     moves.append( mvScale(clock_morpho, lambda=1,    weight=4.0) )



     alpha_morpho ~ dnUniform( 0, 1E6 )

     rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )

     #Moves on the parameters of the Gamma distribution.

     moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))





     n_max_states <- 3
     idx = 1
     for (i in 1:n_max_states) {
            morpho_bystate[i] <- morpho
        morpho_bystate[i].setNumStatesPartition(i)
         nc = morpho_bystate[i].nchar()
         # for non-empty character blocks
         if (nc > 0) {
             # make i-by-i rate matrix
             q[idx] <- fnJC(i)
     # create model of evolution for the character block
             m_morph[idx] ~ dnPhyloCTMC( tree=fbd_tree,
                                         Q=q[idx],
                                         nSites=nc,
                                         siteRates=rates_morpho,
                                         branchRates=clock_morpho,
                                         type="Standard")

             # attach the data
     	    m_morph[idx].clamp(morpho_bystate[i])

             # increment counter
             idx = idx + 1
     idx
     }
     }

 #    n_max_states <- 4
  #   for (i in 1:n_max_states) {
   #                 morpho_nf_bystate[i] <- morpho_nf
    #morpho_nf_bystate[i].setNumStatesPartition(i)
     #    nc = morpho_nf_bystate[i].nchar()
         # for non-empty character blocks
      #   if (nc > 0) {
             # make i-by-i rate matrix
       #      q[idx] <- fnJC(i)
     # create model of evolution for the character block
        #     m_morph[idx] ~ dnPhyloCTMC( tree=fbd_tree,
         #                                Q=q[idx],
          #                               nSites=nc,
           #                              siteRates=rates_morpho,
            #                             branchRates=clock_morpho,
             #                            type="Standard")

             # attach the data
     	    #m_morph[idx].clamp(morpho_nf_bystate[i])

             # increment counter
             #idx = idx + 1
     #idx
     #}
     #}



     mymodel = model(fbd_tree)



     monitors.append(mnModel(filename="output/para.log", printgen=10))



     monitors.append(mnFile(filename="output/para.trees", printgen=10, fbd_tree))



     monitors.append(mnScreen(printgen=10, num_samp_anc, origin_time))



     mymcmc = mcmc(mymodel, monitors, moves)


     mymcmc.run(generations=1000000)



     q()
