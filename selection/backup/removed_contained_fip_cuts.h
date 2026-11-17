/**
 * @file removed_contained_fip_cuts.h
 * @brief BACKUP ONLY - not included or compiled anywhere.
 * @details Photon-only first-interaction-point (FIP) containment cuts that
 * were removed from include/cuts.h on 2026-11-16 in favour of the gOre
 * (photon or electron) versions, contained_gore_fip_cut and
 * contained_gore_fip_cut_depot.
 *
 * Difference from the gOre versions: these check the leading and sub-leading
 * PHOTON only (selectors::leading_photon / sub_leading_photon, or the
 * deposited-energy equivalents); electrons are ignored entirely.
 *
 * To restore: paste the two functions back into namespace cuts in
 * include/cuts.h (they lived just after fiducial_cut). The code below is
 * verbatim from cuts.h at the time of removal.
 */

namespace cuts
{
    /**
     * @brief Apply a fiducial volume cut on photons based on their first interaction
     * point (FIP) and kinetic energy threshold.
     * @details Classifies an interaction as contained or not based on whether
     * the leading (and sub-leading, if present) photon's first interaction point
     * lies within the TPC fiducial margin. Only photons with a kinetic energy
     * at or above 20 MeV are subject to the containment check; photons below
     * this threshold are treated as contained regardless of their FIP. If no
     * photons are found in the interaction, the interaction is also treated as
     * contained.
     *
     * The containment check is delegated to @ref pcuts::fip_contained with a
     * margin parameter of 1.0 cm.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the photons fip is within fiducial volume
     */
    template<class T>
    bool contained_fip_cut(const T & obj, std::vector<double> params){
        if(params.size() != 1)
            throw std::invalid_argument("contained_fip_cut requires exactly one parameter: the inset margin (cm).");
        size_t lead_phot_index = selectors::leading_photon(obj);
        size_t sub_lead_phot_index = selectors::sub_leading_photon(obj);
        // want it to be contained whenever there are no photons or the photons are too low of an energy 
        // not contained whenever the photons with suitable energy have fip within the TPC margin

        if (lead_phot_index == kNoMatch){
            return true;
        }

        else if (sub_lead_phot_index == kNoMatch) {
            auto & lead_phot(obj.particles[lead_phot_index]);
            if (lead_phot.ke>=20){
                if (pcuts::fip_contained(lead_phot, params)==1){
                return true;
                }
                else{
                    return false;
                }
            }
            else{
                return true;
            }
        }
        else {
            auto & lead_phot(obj.particles[lead_phot_index]);
            auto & sub_lead_phot(obj.particles[sub_lead_phot_index]);
            if (lead_phot.ke>=20 && sub_lead_phot.ke>=20){
                if (pcuts::fip_contained(lead_phot,params)==1 && pcuts::fip_contained(sub_lead_phot,params)==1){
                    return true;
                }
                else{
                    return false;
                }
            }
            else if (lead_phot.ke>=20){
                if (pcuts::fip_contained(lead_phot,params)==1){
                return true;
                }
                else{
                    return false;
                }
            }
            else{
                return true;
            }
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, contained_fip_cut, contained_fip_cut);


    /**
     * @brief Apply a fiducial volume cut on photons based on their first interaction
     * point (FIP) and deposit_sum threshold.
     * @details Classifies an interaction as contained or not based on whether
     * the leading (and sub-leading, if present) photon's first interaction point
     * lies within the TPC fiducial margin. Only photons with a deposit_sum
     * at or above 20 MeV are subject to the containment check; photons below
     * this threshold are treated as contained regardless of their FIP. If no
     * photons are found in the interaction, the interaction is also treated as
     * contained.
     *
     * The containment check is delegated to @ref pcuts::fip_contained with a
     * margin parameter of 1.0 cm.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the photons fip is within fiducial volume
     */
    template<class T>
    bool contained_fip_cut_depot(const T & obj, std::vector<double> params){
        if(params.size() != 1)
            throw std::invalid_argument("contained_fip_cut_depot requires exactly one parameter: the inset margin (cm).");
        size_t lead_phot_index = selectors::leading_depot_photon(obj);
        size_t sub_lead_phot_index = selectors::sub_leading_depot_photon(obj);
        // want it to be contained whenever there are no photons or the photons are too low of an energy 
        // not contained whenever the photons with suitable energy have fip within the TPC margin

        if (lead_phot_index == kNoMatch){
            return true;
        }

        else if (sub_lead_phot_index == kNoMatch) {
            auto & lead_phot(obj.particles[lead_phot_index]);
            if (lead_phot.depositions_sum>=20){
                if (pcuts::fip_contained(lead_phot,params)==1){
                return true;
                }
                else{
                    return false;
                }
            }
            else{
                return true;
            }
        }
        else {
            auto & lead_phot(obj.particles[lead_phot_index]);
            auto & sub_lead_phot(obj.particles[sub_lead_phot_index]);
            if (lead_phot.depositions_sum>=20 && sub_lead_phot.depositions_sum>=20){
                if (pcuts::fip_contained(lead_phot,params)==1 && pcuts::fip_contained(sub_lead_phot,params)==1){
                    return true;
                }
                else{
                    return false;
                }
            }
            else if (lead_phot.depositions_sum>=20){
                if (pcuts::fip_contained(lead_phot,params)==1){
                return true;
                }
                else{
                    return false;
                }
            }
            else{
                return true;
            }
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, contained_fip_cut_depot, contained_fip_cut_depot);

} // namespace cuts
