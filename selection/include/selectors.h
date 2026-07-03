/**
 * @file selectors.h
 * @brief Header file for the selectors used in the SPINE analysis framework.
 * @details This file contains the definitions of selectors which can be used
 * to select a single particle (by index) within an interaction. This is a
 * useful feature for reducing down a collection of particles in a final state
 * to just a single one, which allows a user to broadcast a particle-level
 * variable "upwards" to the interaction level (i.e., it can be placed in a
 * branch of a tree that is otherwise filled with interaction-level variables).
 * @author mueller@fnal.gov
 */
#ifndef SELECTORS_H
#define SELECTORS_H
#include <vector>

#include "framework.h"
#include "include/particle_cuts.h"
#include "include/particle_variables.h"

/**
 * @namespace selectors
 * @brief Namespace for organizing selectors which act on interactions.
 * @details This namespace is intended to be used for organizing selectors
 * which act on interactions. Each selector is implemented as a function
 * which takes an interaction object as an argument and returns the index
 * of the selected particle. The function should be templated on the type
 * of interaction object if the selector is intended to be used on both
 * true and reconstructed interactions.
 */
namespace selectors
{
    /**
     * @brief Finds the index corresponding to the leading particle of the
     * specified particle type.
     * @details The leading particle is defined as the particle with the
     * highest kinetic energy. The method of calculating kinetic energy is
     * inherited by the @ref pvars::ke function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
     */
    template <class T>
    size_t leading_particle_index(const T & obj, uint16_t pid)
    {
        double leading_ke(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::ke(p));
            if(pvars::pid(p) == pid && energy > leading_ke)
            {
                leading_ke = energy;
                index = i;
            }
        }
        return index;
    }

    /**
     * @brief Finds the index corresponding to the leading particle of the
     * specified type, ranked by deposited energy rather than kinetic energy.
     * @details Identical in structure to @ref leading_particle_index, but
     * uses @ref pvars::total_depositions (the sum of non-ghost spacepoint
     * depositions) as the energy metric rather than kinetic energy. This
     * variant is preferred in reco-level selections where deposited energy
     * is the primary observable and kinetic energy may be unreliable.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid the SPINE PID index of the particle type to select.
     * @return the index of the particle of the specified type with the
     * highest deposited energy, or kNoMatch if none is found.
     */    
    template <class T>
    size_t leading_depot_index(const T & obj, uint16_t pid)
    {
        double leading_depot(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::total_depositions(p));
            if(pvars::pid(p) == pid && energy > leading_depot)
            {
                leading_depot = energy;
                index = i;
            }
        }
        return index;
    }

    /**
     * @brief Finds the index corresponding to the sub-leading particle of the
     * specified type, ranked by kinetic energy.
     * @details The sub-leading particle is defined as the particle with the
     * second-highest kinetic energy among all particles of the specified PID
     * in the interaction. The kinetic energy is calculated via @ref pvars::ke.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid the SPINE PID index of the particle type to select.
     * @return the index of the sub-leading particle (second-highest KE), or
     * kNoMatch if fewer than two particles of the specified type are found.
     */
    template <class T>
    size_t sub_leading_particle_index(const T & obj, uint16_t pid)
    {
        double sub_leading_ke(0);
        double leading_ke(0);
        size_t index_l(kNoMatch);
        size_t index_s(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::ke(p));
            if(pvars::pid(p) == pid && energy > sub_leading_ke)
            {
                if(energy > leading_ke)
                {
                    sub_leading_ke = leading_ke;
                    index_s = index_l;

                    leading_ke = energy;
                    index_l = i;
                }
                else{
                    sub_leading_ke = energy;
                    index_s = i;
                }
                
            }
        }
        return index_s;
    }

    /**
     * @brief Finds the index corresponding to the sub-leading particle of the
     * specified type, ranked by deposited energy.
     * @details Identical in structure to @ref sub_leading_particle_index, but
     * uses @ref pvars::total_depositions (the sum of non-ghost spacepoint
     * depositions) as the energy metric rather than kinetic energy. This
     * variant is preferred in reco-level selections where deposited energy is
     * the primary observable.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid the SPINE PID index of the particle type to select.
     * @return the index of the sub-leading particle (second-highest deposited
     * energy), or kNoMatch if fewer than two particles of the specified type
     * are found.
     */
    template <class T>
    size_t sub_leading_depot_index(const T & obj, uint16_t pid)
    {
        double sub_leading_depot(0);
        double leading_depot(0);
        size_t index_l(kNoMatch);
        size_t index_s(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::total_depositions(p));
            if(pvars::pid(p) == pid && energy > sub_leading_depot)
            {
                if(energy > leading_depot)
                {
                    sub_leading_depot = leading_depot;
                    index_s = index_l;

                    leading_depot = energy;
                    index_l = i;
                }
                else{
                    sub_leading_depot = energy;
                    index_s = i;
                }
                
            }
        }
        return index_s;
    }

    /**
     * @brief Finds the index corresponding to the longest track.
     * @details The longest track is defined as the track with the longest
     * length, which is calculated upstream in SPINE. The particle instance is
     * required to have a semantic type of 1 (track) and have a start point
     * within 6 cm of the interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the longest track.
     */
    template<class T>
    size_t longest_track(const T & obj)
    {
        double longest_length(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];

            // Distance between interaction vertex and particle start.
            double vertex_distance = std::sqrt(
                std::pow(pvars::start_x(p) - obj.vertex[0], 2) +
                std::pow(pvars::start_y(p) - obj.vertex[1], 2) +
                std::pow(pvars::start_z(p) - obj.vertex[2], 2)
            );

            // Skip particles that are not tracks or are too far from the
            // interaction vertex.
            if(pvars::semantic_type(p) != 1 || vertex_distance >= 6)
                continue;
            
            // Update the longest length and index if the current particle
            // is longer than the longest found so far.
            if(pvars::length(p) > longest_length)
            {
                longest_length = pvars::length(p);
                index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(longest_track, longest_track);

    /**
     * @brief Finds the index corresponding to the second longest track.
     * @details The second longest track is defined as the track with the
     * second longest length, which is calculated upstream in SPINE. The
     * particle instance is required to have a semantic type of 1 (track) and
     * have a start point within 6 cm of the interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the second longest track.
    */
    template<class T>
    size_t second_longest_track(const T & obj)
    {
        double longest_length(0);
        double second_longest_length(0);
        size_t index(kNoMatch), second_index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];

            // Distance between interaction vertex and particle start.
            double vertex_distance = std::sqrt(
                std::pow(pvars::start_x(p) - obj.vertex[0], 2) +
                std::pow(pvars::start_y(p) - obj.vertex[1], 2) +
                std::pow(pvars::start_z(p) - obj.vertex[2], 2)
            );

            // Skip particles that are not tracks or are too far from the
            // interaction vertex.
            if(pvars::semantic_type(p) != 1 || vertex_distance >= 6)
                continue;

            // Check if the current particle is longer than the longest found
            // so far. If so, update the longest and second longest lengths.
            if(pvars::length(p) > longest_length)
            {
                second_longest_length = longest_length;
                longest_length = pvars::length(p);
                second_index = index;
                index = i;
            }

            // If the current particle is not longer than the longest but
            // is longer than the second longest, update the second longest.
            else if(pvars::length(p) > second_longest_length)
            {
                second_longest_length = pvars::length(p);
                second_index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(second_longest_track, second_longest_track);
    
    /**
     * @brief Finds the index corresponding to the sub-leading photon, ranked
     * by deposited energy.
     * @details Returns the index of the photon with the second-highest
     * deposited energy (depositions_sum) in the interaction, via
     * @ref sub_leading_depot_index. This is the deposited-energy analogue of
     * @ref sub_leading_photon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the sub-leading photon by deposited energy, or
     * kNoMatch if fewer than two photons are found.
     */
    template<class T>
    size_t sub_leading_depot_photon(const T & obj)
    {
        return sub_leading_depot_index(obj, pvars::kPhoton);
    }
    REGISTER_SELECTOR(sub_leading_depot_photon, sub_leading_depot_photon);

    /**
     * @brief Finds the index corresponding to the leading and sub-leading shower, 
     * ranked by deposited energy.
     * @details Returns the index of the shower with the highest and second-highest
     * deposited energy (depositions_sum) in the interaction
     * @note created for single shower selection
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading and sub-leading shower by deposited energy, or
     * kNoMatch if showers are not found.
     */
    template<class T>
    std::tuple<size_t, size_t> lead_and_sub_shower_index_depot(const T & obj)
    {
        double sub_photon_depot(0);
        double leading_photon_depot(0);
        double sub_electron_depot(0);
        double leading_electron_depot(0);

        size_t index_lead_photon(kNoMatch);
        size_t index_sub_photon(kNoMatch);
        size_t index_lead_electron(kNoMatch);
        size_t index_sub_electron(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::total_depositions(p));
            if(pvars::pid(p) == pvars::kPhoton && energy > sub_photon_depot)
            {
                if(energy > leading_photon_depot)
                {
                    sub_photon_depot = leading_photon_depot;
                    index_sub_photon = index_lead_photon;

                    leading_photon_depot = energy;
                    index_lead_photon = i;
                }
                else{
                    sub_photon_depot = energy;
                    index_sub_photon = i;
                }
                
            }

            else if(pvars::pid(p) == pvars::kElectron && energy > sub_electron_depot){
                if(energy > leading_electron_depot)
                {
                    sub_electron_depot = leading_electron_depot;
                    index_sub_electron = index_lead_electron;

                    leading_electron_depot = energy;
                    index_lead_electron = i;
                }
                else{
                    sub_electron_depot = energy;
                    index_sub_electron = i;
                }
                
            }
        }

        size_t index_lead(kNoMatch);
        size_t index_sub(kNoMatch);
        if (leading_photon_depot>leading_electron_depot){
            index_lead = index_lead_photon;
            if (sub_photon_depot>leading_electron_depot){
                index_sub = index_sub_photon;
            }
            else{
                index_sub = index_lead_electron;
            }
        }
        else{
            index_lead = index_lead_electron;
            if (sub_electron_depot>leading_photon_depot){
                index_sub = index_sub_electron;
            }
            else{
                index_sub = index_lead_photon;
            }
        }
        
        return {index_lead, index_sub};
    }

    /**
     * @brief Finds the index corresponding to the two highest deposition showers 
     * between the max and min threshold ranked by deposited energy.
     * @details Returns the index of the shower with the highest and second-highest
     * deposited energy (depositions_sum) bewtween the two energy thresholds in the interaction
     * @note created for single shower selection - used to remove electrons between 10-20MeV
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the two highest energy showers between a threshold
     */
    template<class T>
    std::tuple<size_t, size_t> mid_and_submid_shower_index_depot(const T & obj)
    {
        double sub_photon_depot(0);
        double leading_photon_depot(0);
        double sub_electron_depot(0);
        double leading_electron_depot(0);

        size_t index_lead_photon(kNoMatch);
        size_t index_sub_photon(kNoMatch);
        size_t index_lead_electron(kNoMatch);
        size_t index_sub_electron(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::total_depositions(p));
            if(pvars::pid(p) == pvars::kPhoton && energy > sub_photon_depot && energy < 20 && energy > 10)
            {
                if(energy > leading_photon_depot)
                {
                    sub_photon_depot = leading_photon_depot;
                    index_sub_photon = index_lead_photon;

                    leading_photon_depot = energy;
                    index_lead_photon = i;
                }
                else{
                    sub_photon_depot = energy;
                    index_sub_photon = i;
                }
                
            }

            else if(pvars::pid(p) == pvars::kElectron && energy > sub_electron_depot && energy < 20 && energy > 10){
                if(energy > leading_electron_depot)
                {
                    sub_electron_depot = leading_electron_depot;
                    index_sub_electron = index_lead_electron;

                    leading_electron_depot = energy;
                    index_lead_electron = i;
                }
                else{
                    sub_electron_depot = energy;
                    index_sub_electron = i;
                }
                
            }
        }

        size_t index_lead(kNoMatch);
        size_t index_sub(kNoMatch);
        if (index_lead_electron == kNoMatch){
            index_lead = index_lead_photon;
            index_sub = index_sub_photon;
        }
        else if(index_sub_electron==kNoMatch){
            index_lead = index_lead_electron;
            index_sub = index_sub_photon;
        }
        else{
            index_lead = index_lead_electron;
            index_sub = index_sub_electron;
        }
        
        return {index_lead, index_sub};
    }


    /**
     * @brief Finds the index corresponding to the leading and sub-leading shower, 
     * ranked by kinetic energy.
     * @details Returns the index of the shower with the highest and second-highest
     * kinetic energy (ke) in the interaction
     * @note created for single shower selection
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading and sub-leading shower by kinetic energy, or
     * kNoMatch if showers are not found.
     */
    template<class T>
    std::tuple<size_t, size_t> lead_and_sub_shower_index(const T & obj)
    {
        double sub_photon(0);
        double leading_photon(0);
        double sub_electron(0);
        double leading_electron(0);

        size_t index_lead_photon(kNoMatch);
        size_t index_sub_photon(kNoMatch);
        size_t index_lead_electron(kNoMatch);
        size_t index_sub_electron(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::ke(p));
            if(pvars::pid(p) == pvars::kPhoton && energy > sub_photon)
            {
                if(energy > leading_photon)
                {
                    sub_photon = leading_photon;
                    index_sub_photon = index_lead_photon;

                    leading_photon = energy;
                    index_lead_photon = i;
                }
                else{
                    sub_photon = energy;
                    index_sub_photon = i;
                }
                
            }

            else if(pvars::pid(p) == pvars::kElectron && energy > sub_electron){
                if(energy > leading_electron)
                {
                    sub_electron = leading_electron;
                    index_sub_electron = index_lead_electron;

                    leading_electron = energy;
                    index_lead_electron = i;
                }
                else{
                    sub_electron = energy;
                    index_sub_electron = i;
                }
                
            }
        }

        size_t index_lead(kNoMatch);
        size_t index_sub(kNoMatch);
        if (leading_photon>leading_electron){
            index_lead = index_lead_photon;
            if (sub_photon>leading_electron){
                index_sub = index_sub_photon;
            }
            else{
                index_sub = index_lead_electron;
            }
        }
        else{
            index_lead = index_lead_electron;
            if (sub_electron>leading_photon){
                index_sub = index_sub_electron;
            }
            else{
                index_sub = index_lead_photon;
            }
        }
        
        return {index_lead, index_sub};
    }


    /**
     * @brief Finds the index corresponding to the leading shower, ranked by
     * kinetic energy.
     * @details Returns the index of the shower with the highest kinetic
     * energy (ke) in the interaction, via
     * @ref lead_and_sublead_shower_index. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading shower by kinetic energy, or
     * kNoMatch if no showers are found.
     */
    template<class T>
    size_t leading_shower(const T & obj)
    {
        auto [index_lead, index_sub] = lead_and_sub_shower_index(obj);
        return index_lead;
    }
    REGISTER_SELECTOR(leading_shower, leading_shower);

    /**
     * @brief Finds the index corresponding to the sub-leading shower, ranked by
     * kinetic energy.
     * @details Returns the index of the shower with the second highest kinetic
     * energy (ke) in the interaction, via
     * @ref lead_and_sublead_shower_index. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the sub-leading shower by kinetic energy, or
     * kNoMatch if no showers are found.
     */
    template<class T>
    size_t sub_leading_shower(const T & obj)
    {
        auto [index_lead, index_sub] = lead_and_sub_shower_index(obj);
        return index_sub;
    }
    REGISTER_SELECTOR(sub_leading_shower, sub_leading_shower);

    /**
     * @brief Finds the index corresponding to the leading shower, ranked by
     * deposited energy.
     * @details Returns the index of the shower with the highest deposited
     * energy (depositions_sum) in the interaction, via
     * @ref lead_and_sublead_shower_index_depot. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading shower by deposited energy, or
     * kNoMatch if no showers are found.
     */
    template<class T>
    size_t leading_depot_shower(const T & obj)
    {
        auto [index_lead, index_sub] = lead_and_sub_shower_index_depot(obj);
        return index_lead;
    }
    REGISTER_SELECTOR(leading_depot_shower, leading_depot_shower);


    template<class T>
    size_t mid_depot_shower(const T & obj)
    {
        auto [index_lead, index_sub] = mid_and_submid_shower_index_depot(obj);
        return index_lead;
    }
    REGISTER_SELECTOR(mid_depot_shower, mid_depot_shower);
    

    template<class T>
    size_t submid_depot_shower(const T & obj)
    {
        auto [index_lead, index_sub] = mid_and_submid_shower_index_depot(obj);
        return index_lead;
    }
    REGISTER_SELECTOR(submid_depot_shower, submid_depot_shower);


    /**
     * @brief Finds the index corresponding to the sub-leading shower, ranked by
     * deposited energy.
     * @details Returns the index of the shower with the second highest deposited
     * energy (depositions_sum) in the interaction, via
     * @ref lead_and_sublead_shower_index_depot. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the sub-leading shower by deposited energy, or
     * kNoMatch if no showers are found.
     */
    template<class T>
    size_t sub_leading_depot_shower(const T & obj)
    {
        auto [index_lead, index_sub] = lead_and_sub_shower_index_depot(obj);
        return index_sub;
    }
    REGISTER_SELECTOR(sub_leading_depot_shower, sub_leading_depot_shower);

    /**
     * @brief Finds the index corresponding to the leading photon, ranked by
     * deposited energy.
     * @details Returns the index of the photon with the highest deposited
     * energy (depositions_sum) in the interaction, via
     * @ref leading_depot_index. This is the deposited-energy analogue of
     * @ref leading_photon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading photon by deposited energy, or
     * kNoMatch if no photons are found.
     */
    template<class T>
    size_t leading_depot_photon(const T & obj)
    {
        return leading_depot_index(obj, pvars::kPhoton);
    }
    REGISTER_SELECTOR(leading_depot_photon, leading_depot_photon);

    /**
     * @brief Finds the index corresponding to the sub-leading photon, ranked
     * by kinetic energy.
     * @details Returns the index of the photon with the second-highest kinetic
     * energy in the interaction, via @ref sub_leading_particle_index.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the sub-leading photon by kinetic energy, or
     * kNoMatch if fewer than two photons are found.
     */
    template<class T>
    size_t sub_leading_photon(const T & obj)
    {
        return sub_leading_particle_index(obj, pvars::kPhoton);
    }
    REGISTER_SELECTOR(sub_leading_photon, sub_leading_photon);

    /**
     * @brief Finds the index corresponding to the leading photon.
     * @details The leading photon is defined as the photon with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading photon (highest KE).
     */
    template<class T>
    size_t leading_photon(const T & obj)
    {
        return leading_particle_index(obj, pvars::kPhoton);
    }
    REGISTER_SELECTOR(leading_photon, leading_photon);

    /**
     * @brief Finds the index corresponding to the leading electron.
     * @details The leading electron is defined as the electron with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading electron (highest KE).
     */
    template<class T>
    size_t leading_electron(const T & obj)
    {
        return leading_particle_index(obj, pvars::kElectron);
    }
    REGISTER_SELECTOR(leading_electron, leading_electron);

    /**
     * @brief Finds the index corresponding to the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading muon (highest KE).
     */
    template<class T>
    size_t leading_muon(const T & obj)
    {
        return leading_particle_index(obj, pvars::kMuon);
    }
    REGISTER_SELECTOR(leading_muon, leading_muon);

    /**
     * @brief Finds the index corresponding to the leading pion.
     * @details The leading pion is defined as the pion with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading pion (highest KE).
     */
    template<class T>
    size_t leading_pion(const T & obj)
    {
        return leading_particle_index(obj, pvars::kPion);
    }
    REGISTER_SELECTOR(leading_pion, leading_pion);
    
    /**
     * @brief Finds the index corresponding to the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading proton (highest KE).
     */
    template<class T>
    size_t leading_proton(const T & obj)
    {
        return leading_particle_index(obj, pvars::kProton);
    }
    REGISTER_SELECTOR(leading_proton, leading_proton);
    
    /**
     * @brief Finds the index corresponding to the target Michel.
     * @details The target Michel is defined as the Michel with the most
     * depositions in the interaction. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the target Michel (largest).
     */
    template<class T>
    size_t target_michel(const T & obj)
    {
        double largest_size(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double size(p.size);
            if(pvars::semantic_type(p) == 2 && size > largest_size)
            {
                largest_size = size;
                index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(target_michel, target_michel);
}
#endif // SELECTORS_H
