#include "../../include/linalg/sparse_vector.h"

#include <algorithm>
#include <limits>

SparseVector::SparseVector() = default;

SparseVector::SparseVector(const std::vector<Scalar> &x_dense) {
    for (Index i = 0; i < Index(x_dense.size()); i++) {
        if (x_dense[i] != Scalar(0)) {
            PushBack(i, x_dense[i]);
        }
    }
}

SparseVector::SparseVector(const Scalar *x_dense, size_t size) {
    for (Index i = 0; i < Index(size); i++) {
        if (x_dense[i] != Scalar(0)) {
            PushBack(i, x_dense[i]);
        }
    }
}

SparseVector SparseVector::operator-() const {
    SparseVector res = *this;
    for (auto &val : res.values) {
        val = -val;
    }
    return res;
}

SparseVector SparseVector::operator*(const Scalar x) const {
    SparseVector res = *this;
    for (auto &val : res.values) {
        val *= x;
    }
    return res;
}

Scalar SparseVector::operator*(const SparseVector &x) const {
    Scalar res(0);

    Index i_le = 0, i_ri = 0;
    while (i_le < this->size() and i_ri < x.size()) {
        auto pos_le = this->index[i_le];
        auto pos_ri = x.index[i_ri];
        if (pos_le < pos_ri) {
            i_le++;
            continue;
        } else if (pos_le > pos_ri) {
            i_ri++;
            continue;
        } else {
            res += this->values[i_le] * x.values[i_ri];
            i_le++;
            i_ri++;
        }
    }

    return res;
}

Scalar SparseVector::operator*(const Scalar *x) const {
    auto dot = Scalar(0);
    for (Index k = 0; k < size(); ++k) {
        dot += values[k] * x[index[k]];
    }
    return dot;
}

SparseVector SparseVector::operator+(const SparseVector &x) const {
    SparseVector res;
    res.Reserve(this->size() + x.size());

    Index i_le = 0, i_ri = 0;
    while (true) {
        auto pos_le = std::numeric_limits<Index>::max();
        auto pos_ri = std::numeric_limits<Index>::max();

        if (i_le < this->size()) {
            pos_le = this->index[i_le];
        }
        if (i_ri < x.size()) {
            pos_ri = x.index[i_ri];
        }

        if (pos_le < pos_ri) {
            res.PushBack(pos_le, this->values[i_le]);
            i_le++;
        } else if (pos_le > pos_ri) {
            res.PushBack(pos_ri, x.values[i_ri]);
            i_ri++;
        } else {
            if (pos_ri == std::numeric_limits<Index>::max()) {
                break;
            }
            if (this->values[i_le] != -x.values[i_ri]) {
                res.PushBack(pos_le, this->values[i_le] + x.values[i_ri]);
            }
            i_le++;
            i_ri++;
        }
    }
    return res;
}

SparseVector SparseVector::operator-(const SparseVector &x) const {
    SparseVector res;
    res.Reserve(this->size() + x.size());

    Index i_le = 0, i_ri = 0;
    while (true) {
        auto pos_le = std::numeric_limits<Index>::max();
        auto pos_ri = std::numeric_limits<Index>::max();

        if (i_le < this->size()) {
            pos_le = this->index[i_le];
        }
        if (i_ri < x.size()) {
            pos_ri = x.index[i_ri];
        }

        if (pos_le < pos_ri) {
            res.PushBack(pos_le, this->values[i_le]);
            i_le++;
        } else if (pos_le > pos_ri) {
            res.PushBack(pos_ri, -x.values[i_ri]);
            i_ri++;
        } else {
            if (pos_ri == std::numeric_limits<Index>::max()) {
                break;
            }
            if (this->values[i_le] != x.values[i_ri]) {
                res.PushBack(pos_le, this->values[i_le] - x.values[i_ri]);
            }
            i_le++;
            i_ri++;
        }
    }
    return res;
}

void SparseVector::operator*=(Scalar x) {
    for (auto &val : values) {
        val *= x;
    }
}

size_t SparseVector::RemoveElements(const std::vector<Index> &indices) {
    // Remove if el.pos is in indices, correspondingly adjust others' el.pos
    // indices must be sorted
    size_t n_removed = 0;
    Index i_remove = 0;
    for (Index j = 0; j < this->size(); j++) {
        Index ind = -1;
        while (true) {
            if (i_remove < Index(indices.size())) {
                ind = indices[i_remove];
                if (ind < index[j]) {
                    i_remove++;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        if (ind == index[j]) {
            n_removed++;
        } else {
            index[j - n_removed] = index[j] - i_remove;
            values[j - n_removed] = values[j];
        }
    }
    Resize(this->size() - n_removed);
    return n_removed;
}

Scalar SparseVector::operator*(const std::vector<Scalar> &x) const {
    auto dot = Scalar(0);
    for (Index k = 0; k < size(); ++k) {
        dot += values[k] * x[index[k]];
    }
    return dot;
}

void SparseVector::Sort() {
    std::vector<std::pair<Index, Scalar>> tmp(size());

    for (Index i = 0; i < Index(size()); i++) {
        tmp[i] = {index[i], values[i]};
    }

    std::sort(tmp.begin(), tmp.end(),
              [](const std::pair<Index, Scalar> &lhs, const std::pair<Index, Scalar> &rhs) {
                  return lhs.first < rhs.first;
              });

    for (Index i = 0; i < Index(size()); i++) {
        index[i] = tmp[i].first;
        values[i] = tmp[i].second;
    }
}
