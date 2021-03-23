#pragma once

#include "cdm/Model.hpp"

namespace cdm {

/*! \brief Class to construct a model from body and joint connections
 * Providing a set of bodies and connecting them together will allow the class to generate the model.
 * Only availablefor open kinematic chain. Closed loop will not work.s
 */
class CDM_DLLAPI ModelConstructor {
public:
    /*! \brief Default constructor. */
    ModelConstructor() = default;

    /*! \brief Add a body connected to a joint to the model.
     * \param joint New joint.
     * \param body Successor body of the joint.
     */
    void addLink(const Joint& joint, const Body& body);
    /*! \brief Connect a lin
     * \param linkName Name of the body to connect.
     * \param nextJointName Joint connected to the body (body becomes a predecessor of the joint).
     * \param T_l_j Transformation from the joint to the body \f${^l}T_j\f$.
     */
    void connectLink(const std::string& bodyLinkName, const std::string& nextJointName, const Transform& T_l_j);
    /*! \brief Build the model.
     * \param rootName Name of the root joint.
     * \param T_w_r Transformation from root to world \f${^w}T_r\f$. Default is the identity.
     * \return The built model.
     */
    Model build(const std::string& rootName, const Transform& T_w_r = Transform::Identity()) const;

private:
    /*! \brief Body basic information */
    struct LinkBase {
        Body body; /*!< Body */
        std::vector<std::pair<Index, Transform>> next; /*!< Body successor indices */
    };
    struct JointBase {
        Joint joint; /*!< Joint */
        Index previous; /*!< Joint predecessor index */
    };

private:
    std::vector<LinkBase> m_bodyBases; /*!< Body basic information */
    std::vector<JointBase> m_jointBases; /*!< Joint basic information */
    std::unordered_map<std::string, Index> m_jointIndexFromName; /*!< Joint index from name */
    std::unordered_map<std::string, Index> m_bodyIndexFromName; /*!< Body index from name */
};

} // namespace cdm
