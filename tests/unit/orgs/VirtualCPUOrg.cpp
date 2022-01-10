/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2022.
 *
 *  @file  VirtualCPUOrg.cpp
 *  @brief Test all functionality of the VirtualCPU organism. 
 */

// CATCH
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
// MABE
#include "orgs/VirtualCPUOrg.hpp"
#include "orgs/instructions/VirtualCPU_Inst_Nop.hpp"
#include "orgs/instructions/VirtualCPU_Inst_IO.hpp"

//
// TODO
//  [ ] GenerateOutput
//    [ ] ???
//  [ ] SetupConfig
//    [ ] ???
//  [ ] ProcessStep calls VCPU Process
//
//  Called in various other test cases
//  [X] SetupModule
//    [X] Manager mutation vars are set
//    [X] Traits are setup
//    [X] SetupInstLib was called
//  [X] SetupInstLib
//    [X] All expected instructions appear in inst. lib.

TEST_CASE("VirtualCPUOrg_Placeholder", "[core]"){
  // Initialize the instruction library, which only needs done once
  mabe::MABE control(0, nullptr);
  control.GetRandom().ResetSeed(100);
  control.AddPopulation("test_pop", 0);
  mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
  // Add NopA, NopB, and NopC
  mabe::VirtualCPU_Inst_Nop nop_inst_module(control);
  // Add IO (for required traits) 
  mabe::VirtualCPU_Inst_IO io_inst_module(control);
  mabe::VirtualCPUOrg tmp_org(manager);
  control.GetTraitManager().Unlock();
  nop_inst_module.SetupModule();
  io_inst_module.SetupModule();
  tmp_org.SetupModule();
  control.GetTraitManager().Lock();
  control.GetTraitManager().Verify(true);
  std::cout << "InstLib size: " << tmp_org.GetInstLib().GetSize() << std::endl;
  emp::DataMap data_map = control.GetOrganismDataMap();
  control.GetTraitManager().RegisterAll(data_map);
  data_map.LockLayout();          
  { // ManagerData
    // [X] Defaults correctly
    // [X] Can be changed
    // [X] Changes are seen by other orgs 
    mabe::MABE control(0, nullptr);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org_a(manager);
    mabe::VirtualCPUOrg org_b(manager);
    // Check that all shared data default properly (other than mutation datastructs)
    CHECK(org_a.SharedData().mut_prob == 0.01);
    CHECK(org_a.SharedData().init_length == 100);
    CHECK(org_a.SharedData().init_random == true);
    CHECK(org_a.SharedData().eval_time == 500);
    CHECK(org_a.SharedData().input_name == "input");
    CHECK(org_a.SharedData().output_name == "output");
    CHECK(org_a.SharedData().merit_name == "merit");
    CHECK(org_a.SharedData().genome_name == "genome");
    CHECK(org_a.SharedData().child_merit_name == "child_merit");
    CHECK(org_a.SharedData().initial_merit == 0); 
    CHECK(org_a.SharedData().verbose == false); 
    CHECK(org_a.SharedData().initial_genome_filename == "ancestor.org");
    CHECK(org_a.SharedData().expanded_nop_args == false);
    // Values can be changed
    org_a.SharedData().mut_prob = 0.05;
    CHECK(org_a.SharedData().mut_prob == 0.05);
    org_a.SharedData().init_length = 1000;
    CHECK(org_a.SharedData().init_length == 1000);
    org_a.SharedData().init_random = false;
    CHECK(org_a.SharedData().init_random == false);
    org_a.SharedData().eval_time = 100;
    CHECK(org_a.SharedData().eval_time == 100);
    org_a.SharedData().input_name = "input_2";
    CHECK(org_a.SharedData().input_name == "input_2");
    org_a.SharedData().output_name = "output_2";
    CHECK(org_a.SharedData().output_name == "output_2");
    org_a.SharedData().merit_name = "merit_2";
    CHECK(org_a.SharedData().merit_name == "merit_2");
    org_a.SharedData().genome_name = "genome_2";
    CHECK(org_a.SharedData().genome_name == "genome_2");
    org_a.SharedData().child_merit_name = "child_merit_2";
    CHECK(org_a.SharedData().child_merit_name == "child_merit_2");
    org_a.SharedData().initial_merit = 1; 
    CHECK(org_a.SharedData().initial_merit == 1); 
    org_a.SharedData().verbose = true; 
    CHECK(org_a.SharedData().verbose == true); 
    org_a.SharedData().initial_genome_filename = "ancestor_expanded.org";
    CHECK(org_a.SharedData().initial_genome_filename == "ancestor_expanded.org");
    org_a.SharedData().expanded_nop_args = true;
    CHECK(org_a.SharedData().expanded_nop_args == true);
    // Value changes are reflected in other organisms
    CHECK(org_b.SharedData().mut_prob == 0.05);
    CHECK(org_b.SharedData().init_length == 1000);
    CHECK(org_b.SharedData().init_random == false);
    CHECK(org_b.SharedData().eval_time == 100);
    CHECK(org_b.SharedData().input_name == "input_2");
    CHECK(org_b.SharedData().output_name == "output_2");
    CHECK(org_b.SharedData().merit_name == "merit_2");
    CHECK(org_b.SharedData().genome_name == "genome_2");
    CHECK(org_b.SharedData().child_merit_name == "child_merit_2");
    CHECK(org_b.SharedData().initial_merit == 1); 
    CHECK(org_b.SharedData().verbose == true); 
    CHECK(org_b.SharedData().initial_genome_filename == "ancestor_expanded.org");
    CHECK(org_b.SharedData().expanded_nop_args == true);
  }
  { // Mutate
    //  [X] Mutate
    //    [X] Returned number of mutations matches changed positions
    //    [X] Mutation probability is used 
    //    [X] Genome trait is updated
    control.GetRandom().ResetSeed(100);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org(manager);
    org.SetupMutationDistribution();
    emp::DataMap data_map = control.GetOrganismDataMap();
    control.GetTraitManager().RegisterAll(data_map);
    data_map.LockLayout();          
    org.SetDataMap(data_map);
    org.Initialize(control.GetRandom());

    while(true){
      std::string original_genome = org.GetGenomeString();
      // Keep mutating organism until we get two mutations at once
      size_t expected_num_muts = org.Mutate(control.GetRandom()); 
      if(expected_num_muts == 2){ 
        std::cout << "Orig. genome:" << std::endl;
        std::cout << original_genome << std::endl;
        std::string new_genome = org.GetGenomeString();
        size_t realized_num_muts = 0;
        for(size_t idx = 5; idx < new_genome.size(); ++idx){
          if(original_genome[idx] != new_genome[idx]){
            realized_num_muts++;
          }
        }
        std::cout << "Mutated genome:" << std::endl;
        std::cout << new_genome << std::endl;
        CHECK(realized_num_muts == 2); // We should have two different instructions!
        break;
      }
    }

    // Ensure mutation probability is used
    std::string original_genome = org.GetGenomeString();
    org.SharedData().mut_prob = 1.0;
    org.SetupMutationDistribution();
    CHECK(org.Mutate(control.GetRandom()) == 100);
    std::string new_genome = org.GetGenomeString();
    size_t realized_num_muts = 0;
    for(size_t idx = 5; idx < new_genome.size(); ++idx){
      if(original_genome[idx] != new_genome[idx]){
        realized_num_muts++;
      }
    }
    CHECK(realized_num_muts == 100);
  }
  {// Randomize
    control.GetRandom().ResetSeed(102);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org(manager);
    org.SetupMutationDistribution();
    emp::DataMap data_map = control.GetOrganismDataMap();
    control.GetTraitManager().RegisterAll(data_map);
    data_map.LockLayout();          
    org.SetDataMap(data_map);
    org.Initialize(control.GetRandom());
    std::string original_genome = org.GetGenomeString();
    org.Randomize(control.GetRandom());
    std::string randomized_genome = org.GetGenomeString();
    size_t realized_num_muts = 0;
    for(size_t idx = 5; idx < randomized_genome.size(); ++idx){
      if(original_genome[idx] != randomized_genome[idx]){
        realized_num_muts++;
      }
    }
    CHECK(realized_num_muts >= 50); // Basic check to make sure its randomized, but 
                                        // some instructions may not actually change
  }
  { // Initialize
    //  [X] If init random is true, init randomly
    //  [X] Else, init to given genome
    //  [X] Genome trait is set
    //  [X] Merit trait is initialized
    //  [X] Child merit is initialized
    //  [X] Underlying VCPU is initialized
    //  [X] Underlying VCPU curates nops
    
    // Random initialization
    control.GetRandom().ResetSeed(102);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org(manager);
    org.SetupMutationDistribution();
    emp::DataMap data_map = control.GetOrganismDataMap();
    control.GetTraitManager().RegisterAll(data_map);
    data_map.LockLayout();          
    org.SetDataMap(data_map);
    org.SharedData().init_length = 50;
    org.SharedData().initial_merit = 20;;
    org.Initialize(control.GetRandom());
    CHECK(org.GetGenomeSize() == 50);
    CHECK(org.GetTrait<std::string>("genome") == org.GetGenomeString());
    CHECK(org.GetTrait<double>("merit") == 20);
    CHECK(org.GetTrait<double>("child_merit") == 20);
    CHECK(org.nops_need_curated == false);
    CHECK(org.are_nops_counted == true);
    // File initialization
    mabe::VirtualCPUOrg org_2(manager);
    org_2.SetupMutationDistribution();
    org_2.SetDataMap(data_map);
    org_2.SharedData().init_random = false;
    org_2.SharedData().initial_genome_filename = "org_nops.org";
    org_2.Initialize(control.GetRandom());
    CHECK(org_2.GetGenomeSize() == 50);
    CHECK(org_2.GetTrait<std::string>("genome") == "[50]abcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcab");
    CHECK(org_2.GetTrait<double>("merit") == 20);
    CHECK(org_2.GetTrait<double>("child_merit") == 20);
    CHECK(org_2.nops_need_curated == false);
    CHECK(org_2.are_nops_counted == true);
  }
  {// CloneOrganism 
    //  [X] All offspring traits match original traits
    //  [X] Offspring genome is _identical_
    control.GetRandom().ResetSeed(103);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org(manager);
    org.SharedData().init_random = true;
    org.SharedData().init_length = 50;
    org.SetupMutationDistribution();
    emp::DataMap data_map = control.GetOrganismDataMap();
    control.GetTraitManager().RegisterAll(data_map);
    data_map.LockLayout();          
    org.SetDataMap(data_map);
    org.Initialize(control.GetRandom());
    org.SetTrait<double>("merit", 2);
    org.SetTrait<double>("child_merit", 3);
    org.AdvanceIP(1);
    std::string original_genome = org.GetGenomeString();
    emp::Ptr<mabe::VirtualCPUOrg> child_org = org.CloneOrganism().DynamicCast<mabe::VirtualCPUOrg>();
    std::string child_genome =child_org->GetGenomeString();
    std::cout << "Parent: " << original_genome << std::endl;
    std::cout << "Child:  " << child_genome << std::endl;
    CHECK(original_genome == child_genome);
    CHECK(child_org->GetTrait<double>("merit") == org.GetTrait<double>("merit")); 
    CHECK(child_org->GetTrait<double>("child_merit") == org.SharedData().initial_merit); 
    CHECK(child_org->inst_ptr == 0); 
    child_org.Delete();
  }
  { // MakeOffspringOrganism
    // [X] Offspring can be mutated
    // [X] Low mutation rate -> offspring might not be mutated
    // [X] Offspring traits are set
    //  [X] Parent's child merit -> offspring's merit
    control.GetRandom().ResetSeed(104);
    mabe::OrganismManager<mabe::VirtualCPUOrg> manager(control, "name", "desc");
    mabe::VirtualCPUOrg org(manager);
    org.SharedData().init_random = true;
    org.SharedData().init_length = 50;
    // First offspring -> no mutations
    org.SharedData().mut_prob = 0;
    org.SetupMutationDistribution();
    emp::DataMap data_map = control.GetOrganismDataMap();
    control.GetTraitManager().RegisterAll(data_map);
    data_map.LockLayout();          
    org.SetDataMap(data_map);
    org.Initialize(control.GetRandom());
    org.SetTrait<double>("merit", 2);
    org.SetTrait<double>("child_merit", 3);
    std::string original_genome = org.GetGenomeString();
    org.AdvanceIP(1);
    emp::Ptr<mabe::VirtualCPUOrg> child_org_1 = 
        org.MakeOffspringOrganism(control.GetRandom()).DynamicCast<mabe::VirtualCPUOrg>();
    std::string child_genome_1 =
        child_org_1->GetGenomeString();
    std::cout << "Parent: " << original_genome << std::endl;
    std::cout << "Child:  " << child_genome_1 << std::endl;
    CHECK(original_genome == child_genome_1);
    CHECK(child_org_1->GetTrait<double>("merit") == org.GetTrait<double>("child_merit")); 
    CHECK(child_org_1->GetTrait<double>("child_merit") == org.SharedData().initial_merit); 
    CHECK(child_org_1->inst_ptr == 0); 
    child_org_1.Delete();
    // Second offspring -> guaranteed mutations
    org.SharedData().mut_prob = 1;
    org.SetupMutationDistribution();
    emp::Ptr<mabe::VirtualCPUOrg> child_org_2 = 
        org.MakeOffspringOrganism(control.GetRandom()).DynamicCast<mabe::VirtualCPUOrg>();
    std::string child_genome_2 =
        child_org_2->GetGenomeString();
    std::cout << "Parent: " << original_genome << std::endl;
    std::cout << "Child:  " << child_genome_2 << std::endl;
    CHECK(original_genome != child_genome_2);
    CHECK(child_org_2->GetTrait<double>("merit") == org.GetTrait<double>("child_merit")); 
    CHECK(child_org_2->GetTrait<double>("child_merit") == org.SharedData().initial_merit); 
    CHECK(child_org_2->inst_ptr == 0); 
    child_org_2.Delete();
  }
}
