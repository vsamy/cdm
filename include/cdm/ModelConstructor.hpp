#pragma once

#include "cdm/Model.hpp"

namespace cdm {

/*! \brief Class to construct a model from link and joint connections
 * Providing a set of links and connecting them together will allow the class to generate the model.
 * Only availablefor open kinematic chain. Closed loop will not work.s
 */
class CDM_DLLAPI ModelConstructor {
public:
    /*! \brief Default constructor. */
    ModelConstructor() = default;

    /*! \brief Add a link connected to a joint to the model.
     * \param joint New joint.
     * \param link Successor link of the joint.
     */
    void addPair(const Joint& joint, const Link& link);
    /*! \brief Connect a lin
     * \param linkName Name of the link to connect.
     * \param nextJointName Joint connected to the link (link becomes a predecessor of the joint).
     * \param T_l_j Transformation from the joint to the link \f${^l}T_j\f$.
     */
    void connect(const std::string& linkName, const std::string& nextJointName, const Transform& T_l_j);
    /*! \brief Build the model.
     * \param rootName Name of the root joint.
     * \param T_w_r Transformation from root to world \f${^w}T_r\f$. Default is the identity.
     * \return The built model.
     */
    Model build(const std::string& rootName, const Transform& T_w_r = Transform::Identity()) const;

private:
    /*! \brief Link basic information */
    struct LinkBase {
        Link link; /*!< Link */
        std::vector<std::pair<Index, Transform>> next; /*!< Link successor indices */
    };
    struct JointBase {
        Joint joint; /*!< Joint */
        Index previous; /*!< Joint predecessor index */
    };

private:
    std::vector<LinkBase> m_linkBases; /*!< Link basic information */
    std::vector<JointBase> m_jointBases; /*!< Joint basic information */
    std::unordered_map<std::string, Index> m_jointIndexFromName; /*!< Joint index from name */
    std::unordered_map<std::string, Index> m_linkIndexFromName; /*!< Link index from name */
};

} // namespace cdm
