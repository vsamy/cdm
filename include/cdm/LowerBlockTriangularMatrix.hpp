#pragma once

#include <Eigen/Core>
#include <vector>

namespace cdm {

template <typename EigenType>
class LowerBlockTriangularMatrix {
    using block_t = Eigen::MatrixBase<EigenType>;
    static_assert(EigenType::RowsAtCompileTime != Eigen::Dynamic && EigenType::ColsAtCompileTime != Eigen::Dynamic, "Only available for fixed-size matrix");

public:
    LowerBlockTriangularMatrix() = default;
    explicit LowerBlockTriangularMatrix(int size1d)
        : m_size1d(size1d)
        , m_rows((std::sqrt(1 + 8 * size1d) - 1) / 2) // Solve n^2 + n - 2 * size1d = 0 and takepositive value
        , m_cols(m_rows)
        , m_blocks(m_size1d)
    {}
    LowerBlockTriangularMatrix(Eigen::Index rows, Eigen::Index cols)
        : m_size1d(rows * (rows + 1) / 2)
        , m_rows(rows)
        , m_cols(cols)
        , m_blocks(m_size1d)
    {
        assert(rows == cols); // For now
    }

    const block_t& operator()(Eigen::Index i, Eigen::Index j) const
    {
        return m_blocks[blockIndexFromij(i, j)];
    }
    block_t& operator()(Eigen::Index i, Eigen::Index j)
    {
        return m_blocks[blockIndexFromij(i, j)];
    }
    block_t blockAt(Eigen::Index i, Eigen::Index j)
    {
        assert(i <= m_size2d && j <= m_size2d);
        if (i <= j) {
            return operator()(i, j);
        } else {
            return block_t::Zero();
        }
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
        Eigen::Index curColSize = m_rows;
        Eigen::Index curCol = 0;
        Eigen::Index pastColSize = 0;
        for (Eigen::Index k = 0; k < lhs.blocks.size(); ++k) {
            Eigen::Index e = 0;
            out.blocks[k].setZero();
            for (Eigen::Index j = 0; j <= curCol; ++j) {
                Eigen::Index lhsIndex = pastColSize + j * curColSize - e + curCol;
                out.blocks[k] += lhs.m_blocks[lhsIndex] * rhs.m_blocks[pastColSize + j];
                e += j + 1;
            }
            curCol++;
            if (curCol % curColSize == 0) {
                pastColSize += curColSize;
                curColSize--;
                curCol = 0;
            }
        }
    }

private:
    Eigen::Index blockIndexFromij(Eigen::Index i, Eigen::Index j)
    {
        assert(i <= m_size2d && j <= m_size2d && i <= j);
        return j * (m_size2d - 1) - j * (j - 1) / 2 + i;
    }
    std::pair<Eigen::Index, Eigen::Index> ijFromBlockIndex(Eigen::Index blockInd)
    {
        std::pair<Eigen::Index, Eigen::Index> ij{ blockInd, 0 };
        while (ij.first >= m_rows - ij.second) {
            ij.first += ij.second - m_rows;
            ij.second++;
        }
        return ij;
    }

private:
    Eigen::Index m_size1d;
    Eigen::Index m_rows;
    Eigen::Index m_cols;
    std::vector<block_t> m_blocks;
};

} // namespace cdm
