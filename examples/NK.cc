/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019
 *
 *  @file  NK.cc
 *  @brief Implementation of NK-evolution using MABE.
 */

#include <iostream>

#include "../source/core/MABE.h"
#include "../source/core/DirectEncoding.h"

#include "config/ArgManager.h"
#include "tools/BitVector.h"
#include "tools/Random.h"

#include "../source/evaluate/EvalNK.h"
#include "../source/select/SelectElite.h"

// EMP_BUILD_CONFIG( NKConfig,
//   GROUP(DEFAULT, "Default settings for NK model"),
//   VALUE(K, uint32_t, 10, "Level of epistasis in the NK model"),
//   VALUE(N, uint32_t, 200, "Number of bits in each organisms (must be > K)"), ALIAS(GENOME_SIZE),
//   VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
//   VALUE(POP_SIZE, uint32_t, 1000, "Number of organisms in the popoulation."),
//   VALUE(MAX_GENS, uint32_t, 2000, "How many generations should we process?"),
//   VALUE(MUT_COUNT, uint32_t, 3, "How many bit positions should be randomized?"), ALIAS(NUM_MUTS),
// )

using BitOrg = emp::BitVector;

int main(int argc, char* argv[])
{
  // NKConfig config;
  // config.Read("NK.cfg");

  // auto args = emp::cl::ArgManager(argc, argv);
  // if (args.ProcessConfigOptions(config, std::cout, "NK.cfg", "NK-macros.h") == false) exit(0);
  // if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  const uint32_t N = 30;  //config.N();
  const uint32_t K = 4;   //config.K();
  // [[maybe_unused]] const uint32_t POP_SIZE = config.POP_SIZE();
  // [[maybe_unused]] const uint32_t MAX_GENS = config.MAX_GENS();
  // [[maybe_unused]] const uint32_t MUT_COUNT = config.MUT_COUNT();

  mabe::MABE control;
  control.AddOrganismType< mabe::DirectEncoding<emp::BitVector> >("Bit Orgs");
  control.AddModule<mabe::EvalNK>(N, K, "bits", "fitness");
  control.AddModule<mabe::SelectElite>("fitness");
  control.Setup();
  control.Update(100);

/*
  // emp::EAWorld<BitOrg, emp::FitCacheOff> pop(random, "NKWorld");
  emp::World<BitOrg> pop(random, "NKWorld");
  pop.SetupFitnessFile().SetTimingRepeat(10);
  pop.SetupSystematicsFile().SetTimingRepeat(10);
  pop.SetupPopulationFile().SetTimingRepeat(10);
  pop.SetPopStruct_Mixed(true);
  pop.SetCache();

  // Build a random initial population
  for (uint32_t i = 0; i < POP_SIZE; i++) {
    BitOrg next_org(N);
    for (uint32_t j = 0; j < N; j++) next_org[j] = random.P(0.5);
    pop.Inject(next_org);
  }

  // Setup the mutation function.
  std::function<size_t(BitOrg &, emp::Random &)> mut_fun =
    [MUT_COUNT, N](BitOrg & org, emp::Random & random) {
      size_t num_muts = 0;
      for (uint32_t m = 0; m < MUT_COUNT; m++) {
        const uint32_t pos = random.GetUInt(N);
        if (random.P(0.5)) {
          org[pos] ^= 1;
          num_muts++;
        }
      }
      return num_muts;
    };
  pop.SetMutFun( mut_fun );
  pop.SetAutoMutate();

  std::cout << 0 << " : " << pop[0] << " : " << landscape.GetFitness(pop[0]) << std::endl;

  std::function<double(BitOrg&)> fit_fun =
    [&landscape](BitOrg & org){ return landscape.GetFitness(org); };
  pop.SetFitFun( fit_fun );

  // Loop through updates
  for (uint32_t ud = 0; ud < MAX_GENS; ud++) {
    // Print current state.
    // for (uint32_t i = 0; i < pop.GetSize(); i++) std::cout << pop[i] << std::endl;
    // std::cout << std::endl;

    // Keep the best individual.
    emp::EliteSelect(pop, 1, 1);

    // Run a tournament for the rest...
    TournamentSelect(pop, 5, POP_SIZE-1);
    pop.Update();
    std::cout << (ud+1) << " : " << pop[0] << " : " << landscape.GetFitness(pop[0]) << std::endl;
  }

  // pop.PrintLineage(0);

//  std::cout << MAX_GENS << " : " << pop[0] << " : " << landscape.GetFitness(pop[0]) << std::endl;

  // pop.GetSignalControl().PrintNames();
  */
}
