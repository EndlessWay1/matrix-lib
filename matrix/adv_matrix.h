
// Depedency
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "../dtypes/dtypes.c"
#include "../Arena/arena.h"
#include "mat_base.h"
#include "mat_base.c"


const double EPS = 1E-12;

b32 mat_inverse(matrix* dst, const matrix* src);
b32 mat_adj(matrix* dst, matrix* src);
f128 mat_det(matrix* a);
f128 mat_mul_ij(matrix* a);

u32 _argmax(matrix* mat, u32 r, u32 c, b8 absolute);
void _swap(matrix* mat, u32 r1, u32 r2);

b32 mat_ref(
    const matrix* eq_src, const matrix* sol_src,
    matrix* eq_dst, matrix* sol_dst
);

b32 mat_rref(
    const matrix* eq, const matrix* sol_src, matrix* sol_dst
);

b32 mat_decomp(
    const matrix* equation, matrix* P, matrix* L, matrix* U
);

u32 mat_rank(const matrix* a);

b32 mat_identity(matrix* a);

void mat_transpose(matrix* a);

b32 mat_diagonal(matrix* a, matrix* vec);

u32 mat_sign(matrix* perm);