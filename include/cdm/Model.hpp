#pragma once

#include "cdm/Joint.hpp"
#include "cdm/Link.hpp"

namespace cdm {

/*! \brief Representation of the model.
 * The model is not modifiable but provide access to all the model information.
 */
class CDM_DLLAPI Model {
public:
    /*! \brief Default constructor. */
    Model() = default;
    /*! \brief Construct the model.
     * \param joints Set of model joints in order.
     * \param links Set of model links in order.
     * \param jointParents Mapping of joint parents.
     * \param jointChilds Mapping of joint childs.
     * \param T0 Relative transformation from previous links to joints.
     */
    Model(std::vector<Joint> joints, std::vector<Link> links, std::vector<Index> jointParents,
        std::vector<Index> jointChilds, std::vector<Transform> T0);

    /*! \brief Get number of links. */
    Index nLinks() const noexcept;
    /*! \brief Get number of joint parameters. */
    Index nParam() const noexcept;
    /*! \brief Get number of DoF. */
    Index nDof() const noexcept;
    /*! \brief Get parent link index from joint index. */
    Index jointParent(Index jointIndex) const noexcept;
    /*! \brief Get parent link index from joint index. Can throw. */
    Index jointParentAt(Index jointIndex) const;
    /*! \brief Get child link index from joint index. */
    Index jointChild(Index jointIndex) const noexcept;
    /*! \brief Get child link index from joint index. Can throw. */
    Index jointChildAt(Index jointIndex) const;
    /*! \brief Position of the joint index in vector of parameters size. */
    Index jointPosInParam(Index jointIndex) const;
    /*! \brief Position of the joint index in vector of DoF size. */
    Index jointPosInDof(Index jointIndex) const;
    /*! \brief Get joint from its index. */
    const Joint& joint(Index jointIndex) const noexcept;
    /*! \brief Get joint from its index. Can throw. */
    const Joint& jointAt(Index jointIndex) const;
    /*! \brief Get joint from its index. */
    const Link& link(Index bodyIndex) const noexcept;
    /*! \brief Get joint from its index. Can throw. */
    const Link& linkAt(Index bodyIndex) const;
    /*! \brief Get transformation from previous link to current joint index. */
    const Transform& T0(Index jointIndex) const noexcept;
    /*! \brief Get transformation from previous link to current joint index. Can throw. */
    const Transform& T0At(Index jointIndex) const;
    /*! \brief Get set of all joints. */
    const std::vector<Joint>& joints() const noexcept;
    /*! \brief Get set of all links. */
    const std::vector<Link>& links() const noexcept;
    /*! \brief Get set of all joint parents. */
    const std::vector<Index>& jointParents() const noexcept;
    /*! \brief Get set of all joint childs. */
    const std::vector<Index>& jointChilds() const noexcept;
    /*! \brief Get set of all joint position in vector of parameters size. */
    const std::vector<Index>& jointPosInParam() const noexcept;
    /*! \brief Get set of all joint position in vector of DoF size. */
    const std::vector<Index>& jointPosInDof() const noexcept;
    /*! \brief Get set of all relative transformation from previous links to joints. */
    const std::vector<Transform>& T0s() const noexcept;

private:
    Index m_nLinks; /*!< Number of links */
    Index m_nParam; /*!< Number of joint parameters */
    Index m_nDof; /*!< Number of degrees-of-freedom */

    std::vector<Joint> m_joints; /*!< Set of all joints in order */
    std::vector<Link> m_links; /*!< Set of all links in order */
    std::vector<Index> m_jointParents; /*!< Set of all joint parents */
    std::vector<Index> m_jointChilds; /*!< Set of all joint childs */
    std::vector<Transform> m_T0; /*!< Set of all relative transformation */
    std::vector<Index> m_jointPosInParam; /*!< Set of joint position in vector of parameters */
    std::vector<Index> m_jointPosInDof; /*!< Set of joint position in vector of DoF */
    std::unordered_map<std::string, Index> m_jointIndexByName; /*!< Map of joint name to index */
    std::unordered_map<std::string, Index> m_linkIndexByName; /*!< Map of link name to index */
};

} // namespace cdm
