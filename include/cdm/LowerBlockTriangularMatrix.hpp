#pragma once

#include "cdm/typedefs.hpp"
#include <vector>

namespace cdm {

template <typename BlockType>
class LowerBlockTriangularMatrix {
    using block_t = BlockType;
    static inline constexpr int block_rows = block_t::RowsAtCompileTime;
    static inline constexpr int block_cols = block_t::ColsAtCompileTime;
    using Scalar = typename block_t::Scalar;
    using matX_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    static_assert(block_rows != Eigen::Dynamic && block_cols != Eigen::Dynamic, "Only available for fixed-size matrix");

public:
    LowerBlockTriangularMatrix() = default;
    explicit LowerBlockTriangularMatrix(Eigen::Index size1d)
        : m_size1d(size1d)
        , m_rows((static_cast<Eigen::Index>(std::sqrt(1 + 8 * size1d)) - 1) / 2) // Solve n^2 + n - 2 * size1d = 0 and takepositive value
        , m_cols(m_rows)
        , m_blocks(m_size1d)
    {}
    LowerBlockTriangularMatrix(Eigen::Index rows, Eigen::Index cols)
    {
        resize(rows, cols);
    }

    void resize(Eigen::Index rows, Eigen::Index cols)
    {
        assert(rows == cols); // For now
        m_size1d = rows * (rows + 1) / 2;
        m_rows = rows;
        m_cols = cols;
        m_blocks.resize(m_size1d);
    }
    LowerBlockTriangularMatrix& setZero(Eigen::Index rows, Eigen::Index cols)
    {
        resize(rows, cols);
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            m_blocks[k].setZero();
        }

        return *this;
    }
    LowerBlockTriangularMatrix& setIdentity(Eigen::Index rows, Eigen::Index cols)
    {
        resize(rows, cols);
        Eigen::Index pastCols = 0;
        while (r > 0) {
            m_blocks[pastCols].setIdentity();
            for (Eigen::Index k = 1; k < r; ++k) {
                m_blocks[pastCols + k].setZero();
            }
            pastCols += r;
            r--;
        }

        return *this;
    }

    const block_t& operator()(Eigen::Index k) const { return m_blocks[k]; }
    block_t& operator()(Eigen::Index k) { return m_blocks[k]; }
    const block_t& operator()(Eigen::Index i, Eigen::Index j) const { return m_blocks[blockIndexFromij(i, j)]; }
    block_t& operator()(Eigen::Index i, Eigen::Index j) { return m_blocks[blockIndexFromij(i, j)]; }
    block_t blockAt(Eigen::Index i, Eigen::Index j)
    {
        assert(i <= m_rows && j <= m_cols);
        if (i <= j) {
            return m_blocks[blockIndexFromij(i, j)];
        } else {
            return block_t::Zero();
        }
    }
    matX_t matrix() const
    {
        matX_t out(block_rows * m_rows, block_cols * m_cols);
        out.setZero();
        for (Eigen::Index i = 0; i < m_rows; ++i) {
            for (Eigen::Index j = 0; j <= i; ++j) {
                out.template block<block_rows, block_cols>(i * block_rows, j * block_cols) = m_blocks[blockIndexFromij(i, j)];
            }
        }

        return out;
    }
    LowerBlockTriangularMatrix inverse() const
    {
        LowerBlockTriangularMatrix out;
        out.setZero(m_rows, m_cols);
        for (Eigen::Index i = 0; i < m_rows; ++i) {
            out(i, i) = m_blocks[blockIndexFromij(i, i)].inverse();
            for (Eigen::Index j = 0; j < i; ++j) {
                for (Eigen::Index k = j; k < i; ++k) {
                    out(i, j) -= m_blocks[blockIndexFromij(i, k)] * out(k, j);
                }
                out(i, j) = out(i, i) * out(i, j);
            }
        }

        return out;
    }

    // Operator
    friend LowerBlockTriangularMatrix operator+(LowerBlockTriangularMatrix lhs, const LowerBlockTriangularMatrix& rhs)
    {
        lhs += rhs;
        return lhs;
    }
    friend LowerBlockTriangularMatrix operator-(LowerBlockTriangularMatrix lhs, const LowerBlockTriangularMatrix& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
    LowerBlockTriangularMatrix& operator+=(const LowerBlockTriangularMatrix& rhs)
    {
        assert(m_size1d == rhs.m_size1d);
        for (size_t i = 0; i < m_blocks.size(); ++i) {
            m_blocks[i] += rhs.m_blocks[i];
        }

        return *this;
    }
    template <int NVec>
    LowerBlockTriangularMatrix& operator+=(const DiInertia<NVec>& rhs)
    {
        assert(m_rows == NVec);
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            m_blocks[k] += rhs.block();
        }

        return *this;
    }
    LowerBlockTriangularMatrix& operator-=(const LowerBlockTriangularMatrix& rhs)
    {
        assert(m_size1d == rhs.m_size1d);
        for (size_t i = 0; i < m_blocks.size(); ++i) {
            m_blocks[i] -= rhs.m_blocks[i];
        }

        return *this;
    }
    friend LowerBlockTriangularMatrix operator*(const LowerBlockTriangularMatrix& lhs, const LowerBlockTriangularMatrix& rhs)
    {
        assert(lhs.m_size1d == rhs.m_size1d);
        LowerBlockTriangularMatrix out(lhs.m_size1d);
        Eigen::Index curColSize = lhs.m_rows;
        Eigen::Index curCol = 0;
        Eigen::Index pastColSize = 0;
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            Eigen::Index e = 0;
            out.m_blocks[k].setZero();
            for (Eigen::Index j = 0; j <= curCol; ++j) {
                Eigen::Index lhsIndex = pastColSize + j * curColSize - e + curCol;
                out.m_blocks[k] += lhs.m_blocks[lhsIndex] * rhs.m_blocks[pastColSize + j];
                e += j + 1;
            }
            curCol++;
            if (curCol % curColSize == 0) {
                pastColSize += curColSize;
                curColSize--;
                curCol = 0;
            }
        }

        return out;
    }
    template <int NVec>
    friend LowerBlockTriangularMatrix operator*(const LowerBlockTriangularMatrix& lhs, const DiMotionSubspace<NVec>& rhs)
    {
        assert(lhs.m_rows == NVec);
        LowerBlockTriangularMatrix out(NVec);
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            out.m_blocks[k] = lhs.m_block[k] * rhs.block();
        }

        return out;
    }
    template <int NVec>
    friend LowerBlockTriangularMatrix operator*(const DiMotionSubspace<NVec>& lsh, const LowerBlockTriangularMatrix& rhs)
    {
        assert(lhs.m_rows == NVec);
        LowerBlockTriangularMatrix out(NVec);
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            out.m_blocks[k] = lhs.block() * rhs.m_block[k];
        }

        return out;
    }
    friend std::ostream& operator<<(std::ostream& os, const LowerBlockTriangularMatrix& mat)
    {
        for (Eigen::Index k = 0; k < m_size1d; ++k) {
            auto ij = mat.ijFromBlockIndex(k);
            os << "\nBlock (" << ij.first << "," << ij.second << ")\n" << mat.m_blocks[k];
        }
        return os;
    }
    

private:
    Eigen::Index blockIndexFromij(Eigen::Index i, Eigen::Index j) const
    {
        assert(i <= m_rows && j <= m_cols && j <= i);
        return j * (m_rows - 1) - j * (j - 1) / 2 + i;
    }
    std::pair<Eigen::Index, Eigen::Index> ijFromBlockIndex(Eigen::Index k) const
    {
        std::pair<Eigen::Index, Eigen::Index> ij{ k, 0 };
        while (ij.first >= m_rows - ij.second) {
            ij.first += ij.second - m_rows;
            ij.second++;
        }
        ij.first += ij.second;

        return ij;
    }

private:
    Eigen::Index m_size1d;
    Eigen::Index m_rows;
    Eigen::Index m_cols;
    std::vector<block_t> m_blocks;
};

} // namespace cdm
