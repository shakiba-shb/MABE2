/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021.
 *
 *  @file  FactoryModule.hpp
 *  @brief Base module to manage a selection of objects that share a common configiguration.
 */

#ifndef MABE_FACTORY_MODULE_H
#define MABE_FACTORY_MODULE_H

#include "emp/meta/TypeID.hpp"

#include "../config/Config.hpp"

#include "MABE.hpp"
#include "Module.hpp"

namespace mabe {

  // Pre-declarations...
  class MABE;
  template <typename OBJ_T, typename BASE_T> class FactoryModule;

  /// Base class for factory products that uses "curiously recursive templates" to fill
  /// out default functionality for when you know the derived type.
  template <typename OBJ_T, typename BASE_T>
  class ProductTemplate : public BASE_T {
  public:
    ProductTemplate(ModuleBase & _man) : BASE_T(_man) { ; }

    using obj_t = OBJ_T;
    using manager_t = FactoryModule<OBJ_T, BASE_T>;

    /// Get the manager for this type of organism.
    manager_t & GetManager() {
      return (manager_t &) BASE_T::GetManager();
    }
    const manager_t & GetManager() const {
      return (const manager_t &) BASE_T::GetManager();
    }

    auto & SharedData() { return GetManager().data; }
    const auto & SharedData() const { return GetManager().data; }
  };


  /// @param OBJ_T the object type being managed by the factory.
  /// @param BASE_T the base object category being mnagaed by the factory.
  template <typename OBJ_T, typename BASE_T>
  class FactoryModule : public Module {
    /// Allow factory products to access private shared data in their own manager only.
    friend class ProductTemplate<OBJ_T, BASE_T>;

  private:
    /// Locate the specification for the data that we need for management in the factory module.
    using data_t = typename OBJ_T::ManagerData;
    
    /// Shared data across all objects that use this factory.
    data_t data;

    /// Maintain a prototype for the objects being created.
    emp::Ptr<BASE_T> obj_prototype;

  public:
    FactoryModule(MABE & in_control, const std::string & in_name, const std::string & in_desc="")
      : Module(in_control, in_name, in_desc)
    {
      SetManageMod(); // @CAO should specify what type of object is managed.
      obj_prototype = emp::NewPtr<obj_t>(*this);
    }
    virtual ~FactoryModule() { obj_prototype.Delete(); }

    /// Save the object type that uses this manager.
    using obj_t = OBJ_T;

    /// Also get the TypeID for this object for more run-time type management.
    emp::TypeID GetObjType() const override { return emp::GetTypeID<obj_t>(); }

    /// Create a clone of the provided object; default to using copy constructor.
    emp::Ptr<BASE_T> CloneObject(const BASE_T & obj) override {
      return emp::NewPtr<obj_t>( (const obj_t &) obj );
    }

    /// Create a random object from scratch.  Default to using the obj_prototype object.
    emp::Ptr<BASE_T> Make() override {
      auto obj_ptr = obj_prototype->Clone();
      return obj_ptr;
    }

    /// Create a random object from scratch.  Default to using the obj_prototype object
    /// and then randomize if a random number generator is provided.
    emp::Ptr<BASE_T> Make(emp::Random & random) override {
      auto obj_ptr = obj_prototype->Clone();
      obj_ptr->Initialize(random);
      return obj_ptr;
    }


    void SetupModule() override {
      obj_prototype->SetupModule();
    }

    void SetupDataMap(emp::DataMap & in_dm) override {
      obj_prototype->SetDataMap(in_dm);
    }

    void SetupConfig() override {
      obj_prototype->SetupConfig();
    }

  };

  /// Build a class that will automatically register modules when created (globally)
  template <typename FACTORY_T>
  struct FactoryModuleRegistrar {
    FactoryModuleRegistrar(const std::string & type_name, const std::string & desc) {
      ModuleInfo new_info;
      new_info.name = type_name;
      new_info.desc = desc;
      new_info.init_fun = [desc](MABE & control, const std::string & name) -> ConfigType & {
        return control.AddModule<FACTORY_T>(name, desc);
      };
      GetModuleInfo().insert(new_info);
    }
  };


  /// MACRO for quickly adding new factory modules.
  #define MABE_REGISTER_FACTORY_MODULE(TYPE, BASE_TYPE, DESC) \
        mabe::FactoryModuleRegistrar<FactoryModule<TYPE, BASE_TYPE>> MABE_ ## TYPE ## _Registrar(#TYPE, DESC)

}

#endif
