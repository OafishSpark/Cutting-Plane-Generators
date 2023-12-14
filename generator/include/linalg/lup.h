#ifndef ESOLVER_LINALG_LU_H
#define ESOLVER_LINALG_LU_H

#include <memory>
#include <vector>

#include "../scalars/scalar.h"

class LUP {
   public:
    LUP() : blocked_(false){};
    explicit LUP(bool blocked) : blocked_(blocked){};

    Index lu(Index, Index, Scalar *, Index, Index *);
    Index inverse(Index, Scalar *, Index, const Index *);

   private:
    Index blocked_lu(Index, Index, Scalar *, Index, Index *);
    Index unblocked_lu(Index, Index, Scalar *, Index, Index *);
    Index blocked_inv_u(Index, Scalar *, Index);
    Index unblocked_inv_u(Index, Scalar *, Index);
    void solve_lu_l(Index, Index, const Scalar *, Index, Scalar *, Index);
    void solve_lu_u(Index, Index, const Scalar *, Index, Scalar *, Index);

    Scalar *work_;
    bool blocked_;
};

#endif  // ESOLVER_LINALG_LU_H
