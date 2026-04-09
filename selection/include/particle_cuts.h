/**
 * @file particle_cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select particles. Each cut is implemented as a function which takes an
 * particle object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections/interaction cuts.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_CUTS_H
#define PARTICLE_CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "particle_variables.h"

/**
 * @namespace pcuts
 * @brief Namespace for organizing generic cuts which act on particles.
 * @details This namespace is intended to be used for organizing cuts which act
 * on particles. Each cut is implemented as a function which takes a
 * particle object as an argument and returns a boolean. The function should
 * be templated on the type of particle object if the cut is intended to be
 * used on both true and reconstructed particles.
 */
namespace pcuts
{
    /**
     * @brief The per-particle thresholds for final state signal particles.
     */
    static std::vector<double> final_state_signal_thresholds = {25.0, 25.0, 143.425, 25.0, 50.0};

    /**
     * @brief Check if the particle is a primary particle.
     * @details This function checks if the particle is a primary particle.
     * Primary designation is handled upstream in SPINE and is based on the
     * softmax primary/secondary scores that are assigned to each particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a primary particle.
     */
    template<class T>
    bool is_primary(const T & p)
    {
        return pvars::primary_classification(p) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, is_primary, is_primary);
    
    /**
     * @brief Check if the particle is a secondary particle.
     * @details A particle is classified as secondary if its SPINE softmax
     * primary/secondary score does not meet the upstream threshold for primary
     * classification. Secondary particles are produced in re-interactions or
     * decays downstream of the neutrino interaction vertex, rather than
     * originating at the vertex itself. This is the logical complement of
     * @ref is_primary.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is not classified as primary.
     * @note The primary classification variable is a softmax score; any value
     * other than 1 is treated as secondary.
     */
    template<class T>
    bool is_secondary(const T & p)
    {
        return pvars::primary_classification(p) != 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, is_secondary, is_secondary);

    /**
     * @brief Check if the particle's start point (first interaction point) is
     * contained within the TPC volume by a given margin.
     * @details Returns true if the particle's start point lies strictly inside
     * the SBND TPC active volume, inset by a configurable margin on all six
     * faces.
     *
     * This cut is used to assess the containment of photons and other shower
     * particles based solely on their first interaction point, which is
     * meaningful for particles that may deposit energy across a large volume
     * (e.g. photons converting far from the vertex).
     *
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @param params a single-element vector whose entry is the inset margin
     * (in cm) applied to all TPC faces. Defaults to 3.0 cm.
     * @return true if the particle's start point is within the TPC volume by
     * at least the specified margin on all sides.
     */
    template<class T>
    bool fip_contained(const T & p, std::vector<double> params={3.0}) {
        if(p.start_point[0] > SBND_XMIN + params[0] && p.start_point[0] < SBND_XMAX-params[0] && p.start_point[1] > SBND_YMIN+params[0] && p.start_point[1] < SBND_YMAX-params[0]
                                && p.start_point[2] > SBND_ZMIN + params[0] && p.start_point[2] < SBND_ZMAX - params[0])
            return true;
        else 
            return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, fip_contained, fip_contained); 

    /**
     * @brief Check if the particle is contained.
     * @details This function checks if the particle is contained. The
     * containment of a particle is defined upstream in SPINE using the
     * full set of points in the particle's track. A particle is considered
     * contained if all of its points are within the detector volume by some
     * margin (e.g., 5 cm) and all of its points are reconstructed within the
     * TPC that physically reads out the hits comprising the spacepoints of
     * the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is contained.
     */
    template<class T>
    bool containment_cut(const T & p) { return p.is_contained == 1; }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, containment_cut, containment_cut);

    /**
     * @brief Place a size cut on the particle.
     * @details This function places a size cut on the particle. The size
     * of the particle is defined as the number of spacepoints in the
     * particle. This cut is intended to be used to remove very small
     * particles which are likely to be noise or reconstruction
     * artifacts.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @param params the parameters for the cut. In this case, this sets the
     * minimum number of spacepoints required for the particle to pass
     * the cut. Defaults to 20 spacepoints.
     * @return true if the particle has a size above the threshold.
     */
    template<class T>
    bool size_cut(const T & p, std::vector<double> params={20.0,})
    {
        if(params.size() != 1)
            throw std::invalid_argument("size_cut requires exactly one parameter: the minimum number of spacepoints.");
        return p.size > params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, size_cut, size_cut);

    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details The particle must be marked as a primary particle and have a
     * kinetic energy above a certain threshold. The per-particle thresholds
     * are stored in a global vector and filled with user-provided values
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
    bool final_state_signal(const T & p)
    {
        return pvars::primary_classification(p) && pvars::ke(p) > final_state_signal_thresholds[pvars::pid(p)];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, final_state_signal, final_state_signal);

    /**
     * @brief Check if the particle is throughgoing.
     * @details This function checks if the particle is throughgoing. A
     * throughgoing particle is defined as a particle which has both ends
     * of the track near the boundary of the detector. This is only applicable
     * to tracks as it is somewhat nonsensical for showers.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is throughgoing.
     */
    template<class T>
    bool throughgoing(const T & p)
    {
        utilities::three_vector start_point = {p.start_point[0], p.start_point[1], p.start_point[2]};
        utilities::three_vector end_point = {p.end_point[0], p.end_point[1], p.end_point[2]};
        return pvars::pid(p) > 1 && utilities::near_boundary(start_point) && utilities::near_boundary(end_point);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, throughgoing, throughgoing);

    /**
     * @brief Check if the particle is of the given type.
     * @details This function checks if the particle is of the given type.
     * The type is determined by the particle ID, which is abstracted away
     * using the PIDFUNC function.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @param params the parameters for the cut. In this case, this sets the
     * particle ID to check against. Defaults to 0, which corresponds to
     * a photon.
     */
    template<class T>
    bool is_pid(const T & p, std::vector<double> params={0.0,})
    {
        return pvars::pid(p) == static_cast<size_t>(params[0]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, is_pid, is_pid);

    /**
     * @brief Check if the particle is of the given semantic type.
     * @details This function checks if the particle is of the given semantic
     * type. The semantic type is determined by the shape variable, which is
     * assigned upstream in SPINE based on the pixel-level semantic
     * segmentation.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @param params the parameters for the cut. In this case, this sets the
     * semantic type to check against. Defaults to 0, which corresponds to
     * a shower.
     * @return true if the particle is of the given semantic type.
     */
    template<class T>
    bool is_semantic_type(const T & p, std::vector<double> params={0.0,})
    {
        if(params.size() != 1)
            throw std::invalid_argument("is_semantic_type requires exactly one parameter: the semantic type to check against.");
        return pvars::semantic_type(p) == static_cast<int>(params[0]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::BothParticle, is_semantic_type, is_semantic_type);
}
#endif // PARTICLE_CUTS_H