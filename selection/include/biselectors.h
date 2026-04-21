/**
 * @file biselectors.h
 * @brief Header file for biselectors used in the SPINE analysis framework.
 * @details This file contains definitions of biselectors which can be used
 * to select a pair of particles (by index) within an interaction. This is a
 * useful feature for reducing down a collection of particles in a final state
 * to just a pair, which allows a user to broadcast a two-particle variable
 * "upwards" to the interaction level.
 * @author mueller@fnal.gov
 */
#ifndef BISELECTORS_H
#define BISELECTORS_H
#include <vector>
#include <utility>

#include "framework.h"
#include "include/selectors.h"

/**
 * @namespace biselectors
 * @brief Namespace for organizing biselectors which act on interactions.
 * @details This namespace is intended to be used for organizing biselectors
 * which act on interactions. Each biselector is implemented as a function
 * which takes an interaction object as an argument and returns a pair of
 * indices corresponding to the two selected particles. The function should
 * be templated on the type of interaction object if the biselector is
 * intended to be used on both true and reconstructed interactions.
 */
namespace biselectors
{
    /**
     * @brief Selects the leading muon and leading proton.
     * @details The leading muon and proton are defined as the particles with
     * the highest kinetic energy of their respective types.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {leading_muon, leading_proton}.
     */
    template<class T>
    std::pair<size_t, size_t> muon_proton(const T & obj)
    {
        return { selectors::leading_muon(obj), selectors::leading_proton(obj) };
    }
    REGISTER_BISELECTOR(muon_proton, muon_proton);

    /**
     * @brief Selects the two longest tracks in the interaction.
     * @details The longest and second longest tracks are defined by their
     * track length as calculated upstream in SPINE.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {longest_track, second_longest_track}.
     */
    template<class T>
    std::pair<size_t, size_t> two_longest_tracks(const T & obj)
    {
        return { selectors::longest_track(obj), selectors::second_longest_track(obj) };
    }
    REGISTER_BISELECTOR(two_longest_tracks, two_longest_tracks);

    /**
     * @brief Selects the leading muon and leading pion.
     * @details The leading muon and pion are defined as the particles with
     * the highest kinetic energy of their respective types.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {leading_muon, leading_pion}.
     */
    template<class T>
    std::pair<size_t, size_t> muon_pion(const T & obj)
    {
        return { selectors::leading_muon(obj), selectors::leading_pion(obj) };
    }
    REGISTER_BISELECTOR(muon_pion, muon_pion);

    
    /**
     * @brief Selects the leading and sub-leading shower.
     * @details The leading and sub-leading shower are defined as the particles with
     * the highest and second highest deposited_sum of the combined photon and electron types.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {leading_shower, sub_leading_shower}.
     */
    template<class T>
    std::pair<size_t, size_t> lead_and_sub_shower(const T & obj)
    {
        auto [lead_shower, sub_shower] = selectors::lead_and_sub_shower_index_depot(obj);
        return {lead_shower, sub_shower};
    }
    REGISTER_BISELECTOR(lead_and_sub_shower, lead_and_sub_shower);

}
#endif // BISELECTORS_H
