/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021-2022.
 *
 *  @file  EvalAntagonistic.hpp
 *  @brief MABE Evaluation module for counting the number of ones (or zeros) in an output.
 * 
 *  Developer notes:
 *  - Can allow vals_trait to also be a vector.
 */

#ifndef MABE_EVAL_ANTAGONISTIC_H
#define MABE_EVAL_ANTAGONISTIC_H

#include <emp/math/constants.hpp>

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"

namespace mabe {

  class EvalAntagonistic : public Module {
  private:
    size_t num_vals = 100;                     // Cardinality of the problem space.
    RequiredMultiTrait<double> vals_trait{this, "vals", "Set of values to evaluate.", AsConfig(num_vals)};
    OwnedMultiTrait<double> scores_trait{this, "scores", "Set of scores for each value.", AsConfig(num_vals)};
    OwnedTrait<double> total_trait{this, "total", "A single value totalling all scores."};
    OwnedTrait<size_t> first_trait{this, "first", "Location of first active positions."};
    OwnedTrait<size_t> active_count_trait{this, "active_count", "Number of activation positions."};


  public:
    EvalAntagonistic(mabe::MABE & control,
                   const std::string & name="EvalAntagonistic",
                   const std::string & desc="Evaluate value sets using the antagonistic fitness function.")
      : Module(control, name, desc) { SetEvaluateMod(true); }
    ~EvalAntagonistic() { }

    // Setup member functions associated with this class.
    static void InitType(emplode::TypeInfo & info) {
      info.AddMemberFunction(
        "EVAL",
        [](EvalAntagonistic & mod, Collection orgs) { return mod.Evaluate(orgs); },
        "Evaluate organisms using the specified diagnostic."
      );
    }

    void SetupConfig() override {
      LinkVar(num_vals, "N", "Cardinality of the problem (number of values to analyze)");
    }

    void SetupModule() override {
      // Nothing needed here yet...
    }

    double Evaluate(Collection orgs) {
      // Track the organism with the highest total score.
      double max_total = 0.0;
      emp::Ptr<Organism> max_org = nullptr;

      // Loop through the living organisms in the target collection to evaluate each.
      mabe::Collection alive_collect( orgs.GetAlive() );
      for (Organism & org : alive_collect) {        
        // Make sure this organism has its values ready for us to access.
        org.GenerateOutput();

        // Get access to the data_map elements that we need.
        std::span<double> vals = vals_trait(org);
        std::span<double> scores = scores_trait(org);
        double & total_score = total_trait(org);
        size_t & first_active = first_trait(org);
        size_t & active_count = active_count_trait(org);

        // Initialize output values.
        total_score = 0.0;
        size_t pos = 0;

        // Only count highest value
        pos = emp::FindMaxIndex(vals);  // Find the sole active position.
        total_score = scores[pos] = vals[pos];
        first_active = pos;
        active_count = 1;

        // All others are subtracted from max and divided by two, creating a
        // pressure to minimize.
        for (size_t i = 0; i < vals.size(); i++) {
          // if (i == pos) {
          //   scores[i] = vals[i];
          // } else {
          scores[i] = vals[i] - (emp::Sum(vals) / 2.0) + vals[i]/2;
          // }
          total_score += scores[i];
        }

        if (total_score > max_total || !max_org) {
          max_total = total_score;
          max_org = &org;
        }
      }
      return max_total;
    }
  };

  MABE_REGISTER_MODULE(EvalAntagonistic, "Evaluate set of values with the antagonistic problem.");
}

#endif
