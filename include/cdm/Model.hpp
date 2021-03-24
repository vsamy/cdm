#pragma once

#include "cdm/Body.hpp"
#include "cdm/Joint.hpp"

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
     * \param bodies Set of model bodies in order.
     * \param jointParents Mapping of joint parents.
     * \param jointChilds Mapping of joint childs.
     * \param T0 Relative transformation from previous bodies to joints.
     */
    Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<Index> jointParents,
        std::vector<Index> jointChilds, std::vector<Transform> T0);

    /*! \brief Get number of bodies. */
    inline Index nLinks() const noexcept { return m_nLinks; }
    /*! \brief Get number of joint parameters. */
    inline Index nParam() const noexcept { return m_nParam; }
    /*! \brief Get number of DoF. */
    inline Index nDof() const noexcept { return m_nDof; }
    /*! \brief Get parent body index from joint index. */
    inline Index jointParent(Index jointIndex) const noexcept { return m_jointParents[jointIndex]; }
    /*! \brief Get parent body index from joint index. Can throw. */
    inline Index jointParentAt(Index jointIndex) const { return m_jointParents.at(jointIndex); }
    /*! \brief Get child body index from joint index. */
    inline Index jointChild(Index jointIndex) const noexcept { return m_jointChilds[jointIndex]; }
    /*! \brief Get child body index from joint index. Can throw. */
    inline Index jointChildAt(Index jointIndex) const { return m_jointChilds.at(jointIndex); }
    /*! \brief Position of the joint index in vector of parameters size. */
    inline Index jointPosInParam(Index jointIndex) const { return m_jointPosInParam.at(jointIndex); }
    /*! \brief Position of the joint index in vector of DoF size. */
    inline Index jointPosInDof(Index jointIndex) const { return m_jointPosInDof.at(jointIndex); }
    /*! \brief Get body index from its name. */
    inline Index bodyIndexByName(const std::string& bodyName) const noexcept { return m_bodyIndexByName.find(bodyName)->second; }
    /*! \brief Get body index from its name. */
    inline Index bodyIndexByNameAt(const std::string& bodyName) const { return m_bodyIndexByName.at(bodyName); }
    /*! \brief Get joint index from its name. */
    inline Index jointIndexByName(const std::string& jointName) const noexcept { return m_jointIndexByName.find(jointName)->second; }
    /*! \brief Get joint index from its name. */
    inline Index jointIndexByNameAt(const std::string& jointName) const { return m_jointIndexByName.at(jointName); }
    /*! \brief Get joint from its index. */
    inline const Joint& joint(Index jointIndex) const noexcept { return m_joints[jointIndex]; }
    /*! \brief Get joint from its index. Can throw. */
    inline const Joint& jointAt(Index jointIndex) const { return m_joints.at(jointIndex); }
    /*! \brief Get joint from its index. */
    inline const Body& body(Index bodyIndex) const noexcept { return m_bodies[bodyIndex]; }
    /*! \brief Get joint from its index. Can throw. */
    inline const Body& bodyAt(Index bodyIndex) const { return m_bodies.at(bodyIndex); }
    /*! \brief Get transformation from previous body to current joint index. */
    inline const Transform& T0(Index jointIndex) const noexcept { return m_T0[jointIndex]; }
    /*! \brief Get transformation from previous body to current joint index. Can throw. */
    inline const Transform& T0At(Index jointIndex) const { return m_T0.at(jointIndex); }
    /*! \brief Get set of all joints. */
    inline const std::vector<Joint>& joints() const noexcept { return m_joints; }
    /*! \brief Get set of all bodies. */
    inline const std::vector<Body>& bodies() const noexcept { return m_bodies; }
    /*! \brief Get set of all joint parents. */
    inline const std::vector<Index>& jointParents() const noexcept { return m_jointParents; }
    /*! \brief Get set of all joint childs. */
    inline const std::vector<Index>& jointChilds() const noexcept { return m_jointChilds; }
    /*! \brief Get set of all joint position in vector of parameters size. */
    inline const std::vector<Index>& jointPosInParam() const noexcept { return m_jointPosInParam; }
    /*! \brief Get set of all joint position in vector of DoF size. */
    inline const std::vector<Index>& jointPosInDof() const noexcept { return m_jointPosInDof; }
    /*! \brief Get set of all relative transformation from previous bodies to joints. */
    inline const std::vector<Transform>& T0s() const noexcept { return m_T0; }

private:
    Index m_nLinks; /*!< Number of bodies */
    Index m_nParam; /*!< Number of joint parameters */
    Index m_nDof; /*!< Number of degrees-of-freedom */

    std::vector<Joint> m_joints; /*!< Set of all joints in order */
    std::vector<Body> m_bodies; /*!< Set of all bodies in order */
    std::vector<Index> m_jointParents; /*!< Set of all joint parents */
    std::vector<Index> m_jointChilds; /*!< Set of all joint childs */
    std::vector<Transform> m_T0; /*!< Set of all relative transformation */
    std::vector<Index> m_jointPosInParam; /*!< Set of joint position in vector of parameters */
    std::vector<Index> m_jointPosInDof; /*!< Set of joint position in vector of DoF */
    std::unordered_map<std::string, Index> m_jointIndexByName; /*!< Map of joint name to index */
    std::unordered_map<std::string, Index> m_bodyIndexByName; /*!< Map of body name to index */
};

} // namespace cdm
