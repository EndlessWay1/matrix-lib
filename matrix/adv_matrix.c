// #include "adv_matrix.h"

/*

#include "../base.h"
#include "mat_base.h"
#include "mat_base.c"

const double EPS = 1E-9;

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
*/


/*
int main(void){
    
mem_arena* perm_arena = arena_create(GiB(1), MiB(1));
srand(42);

matrix* a = mat_create(perm_arena, 4, 4);
a->data[0] = 3.0f;
a->data[1] = 2.0f;
a->data[2] = -1.0f;
a->data[3] = 1.0f;
a->data[4] = 1.0f;
a->data[5] = 0.0f;
a->data[6] = 1.0f;
a->data[7] = 2.0f;
a->data[8] = 2.0f;
a->data[9] = 1.0f;
a->data[10] = 1.0f;
a->data[11] = -1.0f;
a->data[12] = 1.0f;
a->data[13] = 1.0f;
a->data[14] = 1.0f;
a->data[15] = 0.0f;

printf("Rank(A): %u\n", mat_rank(a));
printf("det(A): %Lf \n", mat_det(a));
printf("Matrix A:\n");

matrix* identity = mat_create(perm_arena, 4, 4);
matrix* lower_T = mat_create(perm_arena, 4, 4);
matrix* upper_T = mat_create(perm_arena, 4, 4);
mat_decomp(a, identity, lower_T, upper_T);

printf("Matrix A:\n");
print_mat(a);


matrix* Perm = mat_create(perm_arena, 6, 6);

mat_identity(Perm);
    _swap(Perm, 0, 1);
    _swap(Perm, 2, 4);
    _swap(Perm, 4, 3);
    
    printf("Matrix I 6x6:\n");
    print_mat(Perm);
    printf("det(Perm): %Lf\n", mat_det(Perm));
    
    printf("Matrix I:\n");
    mat_identity(identity);
    _swap(identity, 0, 3);
    _swap(identity, 0, 2);
    print_mat(identity);
    printf("det(Identity): %Lf\n", mat_det(identity));
    
    
    
    b32 r = mat_mul(a, lower_T, upper_T, 1, 0, 0);
    printf("%s\n", (r)? "Multiplied Succesfully": "Invalid Multiply of Matrix");
    printf("Matrix A:\n");
    print_mat(a);
    
    
    arena_destroy(perm_arena);
    
    return 0;
}
*/


f128 mat_mul_ij(matrix* a){
    
    f128 mul = 1.0f;
    u64 size = (u64)a->row*a->col;
    for(u64 i = 0; i < size; i++){
        mul *= a->data[i];
    }
    return mul;
}

b32 mat_diagonal(matrix* a, matrix* vec){
    if (a->row != a->col) return false;
    u64 size = (u64) vec->col*vec->row;
    if (a->row != size) return false;

    for (u64 i = 0; i < size; i++){
        vec->data[i] = a->data[i*a->col + i];
    }
    return true;
}


b32 mat_identity(matrix* a){
    if (a->col != a->row){
        return false;
    }

    mat_clear(a);
    for (u32 i = 0; i < a->row; i++){
        a->data[i*a->col + i] = 1.0f;
    }
    
    return true;
}

void _decomp_01(matrix* L, matrix* U){
    // ref part
    u32 h = 0; // height pivot
    u32 k = 0; // column pivot
    
    matrix* low = L;
    matrix* e = U;

    while (h < e->row && k < e->col){
        // choose a pivot with the highest value
        u32 h_th = _argmax(e, h, k, 1);

        // no pivot
        if (fabsl(e->data[h_th*e->col + k]) < EPS){
            k++;
        }
        else{
            // swap row
            _swap(e, h, h_th);
            _swap(low, h, h_th);

            for (u32 r = h + 1; r < e->row; r++){
                f128 f = (f128) e->data[e->col*r + k]/e->data[e->col*h + k];
                // modify the eq
                e->data[e->col*r + k] = 0.0f;

                for (u32 c = k + 1; c < e->col; c++){
                    e->data[e->col*r + c] -= e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                }

                for (u32 c = 0; c < low->col; c++){
                    low->data[r*low->col + c] -= low->data[h*low->col + c]*f;
                    low->data[low->col*r + c] = (fabsl(low->data[low->col*r + c]) < EPS) ? 0.0f:low->data[low->col*r + c];
                }
            }
            h++;
            k++;
        }
    }
}

void _decomp_11(matrix* P, matrix* L, matrix* U){
    // ref part
    u32 h = 0; // height pivot
    u32 k = 0; // column pivot
    
    matrix* low = L;
    matrix* e = U;

    while (h < e->row && k < e->col){
        // choose a pivot with the highest value
        u32 h_th = _argmax(e, h, k, 1);

        // no pivot
        if (fabsl(e->data[h_th*e->col + k]) < EPS){
            k++;
        }
        else{
            // swap row
            _swap(e, h, h_th);
            _swap(low, h, h_th);
            _swap(P, h, h_th);

            for (u32 r = h + 1; r < e->row; r++){
                f128 f = (f128) e->data[e->col*r + k]/e->data[e->col*h + k];
                // modify the eq
                e->data[e->col*r + k] = 0.0f;

                for (u32 c = k + 1; c < e->col; c++){
                    e->data[e->col*r + c] -= e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                }

                for (u32 c = 0; c < low->col; c++){
                    low->data[r*low->col + c] -= low->data[h*low->col + c]*f;
                    low->data[low->col*r + c] = (fabsl(low->data[low->col*r + c]) < EPS) ? 0.0f:low->data[low->col*r + c];
                }
            }
            h++;
            k++;
        }
    }
}

void _decomp_10(matrix* P, matrix* U){
    // ref part
    u32 h = 0; // height pivot
    u32 k = 0; // column pivot
    
    matrix* e = U;

    while (h < e->row && k < e->col){
        // choose a pivot with the highest value
        u32 h_th = _argmax(e, h, k, 1);

        // no pivot
        if (fabsl(e->data[h_th*e->col + k]) < EPS){
            k++;
        }
        else{
            // swap row
            _swap(e, h, h_th);
            _swap(P, h, h_th);

            for (u32 r = h + 1; r < e->row; r++){
                f128 f = (f128) e->data[e->col*r + k]/e->data[e->col*h + k];
                // modify the eq
                e->data[e->col*r + k] = 0.0f;

                for (u32 c = k + 1; c < e->col; c++){
                    e->data[e->col*r + c] -= e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                }
                
            }
            h++;
            k++;
        }
    }
}

void _decomp_00(matrix* U){
    // ref part
    u32 h = 0; // height pivot
    u32 k = 0; // column pivot
    
    matrix* e = U;

    while (h < e->row && k < e->col){
        // choose a pivot with the highest value
        u32 h_th = _argmax(e, h, k, 1);

        // no pivot
        if (fabsl(e->data[h_th*e->col + k]) < EPS){
            k++;
        }
        else{
            // swap row
            _swap(e, h, h_th);

            for (u32 r = h + 1; r < e->row; r++){
                f128 f = (f128) e->data[e->col*r + k]/e->data[e->col*h + k];
                // modify the eq
                e->data[e->col*r + k] = 0.0f;

                for (u32 c = k + 1; c < e->col; c++){
                    e->data[e->col*r + c] -= e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                }
                
            }
            h++;
            k++;
        }
    }
}

b32 mat_decomp(
    const matrix* equation, matrix* P, matrix* L, matrix* U
){
    b8 t_val = 0b00;
    if (P){
        t_val |= 0b10;
    }
    if (L){
        t_val |= 0b01;
    }
    
    
    if (equation->col != equation->row) return false;
    if (L){
        if (L->col != L->row) return false;
        if (L->col != equation->col) return false;
        
    }
    if (U){
        if (U->col != U->row) return false;
        if (U->col != equation->col) return false;
    }
    if (P){
        if (P->col != P->row) return false;
        if (P->col != equation->col) return false;
        
    }
    
    mem_arena_temp scratch = arena_scratch_get(NULL, 0);
        
    // copy equation
    matrix* e = mat_create(scratch.arena, equation->col, equation->row);
    matrix* low = (L) ? mat_create(scratch.arena, equation->col, equation->row): NULL;
    matrix* perm = (P) ? mat_create(scratch.arena, equation->col, equation->row): NULL;
    
    mat_copy(e, equation);
    if (P){
        (mat_identity(perm));
    }
    if (L){
        mat_identity(low);
    }

    // ref part
    switch (t_val)
    {
    case 0b11:
        _decomp_11(perm,low,e);
        break;
    case 0b10:
        _decomp_10(perm,e);
        break;
    case 0b01:
        _decomp_01(low,e);
        break;
    case 0b00:
        _decomp_00(e);
        break;
    }

    if (U){
        mat_copy(U, e);
    }
    if (L){

        mat_mul(low, low, perm, 1, 0, 1);
        mat_inverse(L, low);
    }
    if (P){
        mat_transpose(perm);
        mat_copy(P, perm);
    }

    arena_scratch_release(scratch);

    return true;
}

void mat_transpose(matrix* a){

    for (u32 i = 0; i < a->row; i++){
        for (u32 j = i; j < a->col; j++){
            f128 temp = a->data[i*a->col + j];
            a->data[i*a->col + j] = a->data[j*a->col + i];
            a->data[j*a->col + i] = temp;
        }
    }
    u32 temp = a->row;
    a->row = a->col;
    a->col = temp;
}


b32 mat_inverse(matrix* dst, const matrix* src){
    if (dst->col != src->col || dst->row != src->row) return false;
    if (dst->col != dst->row || src->col != src->row) return false;
    if (mat_rank(src) != src->row) return false;
    mem_arena_temp scratch = arena_scratch_get(NULL, 0);
    matrix* c = mat_create(scratch.arena, dst->row, dst->col);
    mat_identity(c);
    b32 r = mat_rref(src, c, dst);

    arena_scratch_release(scratch);
    return r;
}

u32 mat_sign(matrix* perm){
    if (perm->col != perm->row){
        return 0;
    }

    u32 n = 0;
    mem_arena_temp scratch = arena_scratch_get(NULL, 0);

    u32* visited = PUSH_ARRAY(scratch.arena, u32, perm->col);
    memset(visited, 0, sizeof(u32)*perm->col);

    u32 next_visited_ptr = 0;
    u32 curr_ptr = next_visited_ptr;

    visited[next_visited_ptr] = 1;
    u32 length = 0;
    u32 zeros = 0;
    while (next_visited_ptr < perm->col && length < perm->col && zeros < perm->col){
        for (u32 i = 0; i < perm->col; i++){
            if (perm->data[curr_ptr*perm->col + i] == 1.0f){
                // found curr ptr in loop
                if (next_visited_ptr == i){
                    // update next_visited_ptr
                    while(next_visited_ptr < perm->col && visited[next_visited_ptr] != 0){
                        next_visited_ptr++;
                    }
                    if (next_visited_ptr < perm->col){
                        curr_ptr = next_visited_ptr;
                        visited[curr_ptr] = 1;
                    }

                    n += length;
                    length = 0;
                    break;
                }
                
                visited[i] = 1;
                length++;
                curr_ptr = i;
                zeros = 0;
            };

        }

        zeros++;

        
    }

    arena_scratch_release(scratch);
    return n;
}

b32 mat_adj(matrix* dst, matrix* src){
    if(dst->row != dst->col) return false;
    if(src->row != src->col) return false;
    if(src->row != dst->col || src->col != dst->row) return false;

    f128 det = mat_det(src);

    mat_copy(dst, src);
    mat_inverse(dst, dst);

    mat_scale(dst, det);

    return true;
}


f128 mat_det(matrix* a){
    if (a->col != a->row) return 0.0f;

    // copy everything
    mem_arena_temp scrath = arena_scratch_get(NULL, 0);

    matrix* e = mat_create(scrath.arena, a->row, a->col);
    matrix* P = mat_create(scrath.arena, a->row, a->col);

    mat_identity(P);
    mat_copy(e, a);
    
    mat_decomp(e, P, NULL, e);
    // ref 
    // check if the rank is correct
    if (mat_rank(e) < e->row){
        return 0;
    }

    // do the opposite
    u32 h = e->row - 1;
    u32 k = e->col - 1;
    while (h > 0 && k > 0){
        if (fabsl(e->data[h*e->col + k]) < EPS){
            k--;
        }
        else{
            for (i32 r = h - 1; r > -1; r--){
                f128 f = (float) e->data[e->col*r + k]/e->data[h*e->col + k];
                e->data[e->col*r + k] = 0.0f;
                for (i32 c = k - 1; c > -1; c--){
                    // printf("height: %i, column: %i\n", r, c);
                    e->data[e->col*r + c] = e->data[e->col*r + c] - e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0:e->data[e->col*r + c];
                }
            }
        }
        h--;
        k--;
    }

    matrix* unknown = mat_create(scrath.arena, e->col, 1);
    mat_diagonal(e, unknown);

    return mat_mul_ij(unknown)* ((mat_sign(P) % 2) ? -1 : 1);
}


b32 mat_rref(
    const matrix* eq, const matrix* sol_src, matrix* sol_dst
){

    if (!sol_src || !sol_dst) return false;
    if (sol_src->col != sol_dst->col || sol_src->row != sol_dst->row) return false;
    if (eq->col != sol_src->row) return false;


    // copy everything
    mem_arena_temp scrath = arena_scratch_get(NULL, 0);

    matrix* e = mat_create(scrath.arena, eq->row, eq->col);
    matrix* s = mat_create(scrath.arena, sol_src->row, sol_src->col);
    mat_copy(e, eq);
    mat_copy(s, sol_src);
    
    mat_ref(e, s, e, s);
    // ref 
    // check if the rank is correct
    if (mat_rank(e) < e->row){
        return false;
    }

    // do the opposite
    u32 h = e->row - 1;
    u32 k = e->col - 1;
    while (h > 0 && k > 0){
        if (fabsl(e->data[h*e->col + k]) < EPS){
            k--;
        }
        else{
            for (i32 r = h - 1; r > -1; r--){
                f128 f = (float) e->data[e->col*r + k]/e->data[h*e->col + k];
                e->data[e->col*r + k] = 0.0f;
                for (i32 c = k - 1; c > -1; c--){
                    // printf("height: %i, column: %i\n", r, c);
                    e->data[e->col*r + c] = e->data[e->col*r + c] - e->data[e->col*h + c]*f;
                    e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0:e->data[e->col*r + c];
                }
                
                for (u32 c = 0; c < s->col; c++){
                    s->data[r*s->col + c] = s->data[r*s->col + c] - s->data[h*s->col + c]*f;
                    s->data[s->col*r + c] = (fabsl(s->data[s->col*r + c]) < EPS) ? 0:s->data[s->col*r + c];
                }
            }
        }
        h--;
        k--;
    }

    // matrix* unknown = mat_create(scrath.arena, e->col, 1);
    // mat_diagonal(e, unknown);
    
    // for (u32 r = 0; r < unknown->row; r++){
    //     f128 dat = unknown->data[r*unknown->col];
    //     for (u32 c = 0; c < s->col; c++){
    //         s->data[r*s->col + c] *= 1.0f/dat;
    //     }
    // }

    mat_copy(sol_dst, s);

    arena_scratch_release(scrath);
    return true;
}

u32 _argmax(matrix* mat, u32 r, u32 c, b8 absolute){
    u32 r_max = r;
    u32 col = mat->col;
    if (absolute){
        for (u32 rw = r; rw < mat->row; rw++){
            f128 data = fabsl(mat->data[col*rw + c]);
            f128 comp_data = fabsl(mat->data[col*r_max + c]);
            if (comp_data < data){
                r_max = rw;
            }
        }
    }
    else{
        for (u32 rw = r; rw < mat->row; rw++){
            if (mat->data[col*rw + c] > mat->data[col*r_max + c]){
                r_max = rw;
            }
        }
    }
    return r_max;
}

void _swap(matrix* mat, u32 r1, u32 r2){
    u32 col = mat->col;
    for (u32 i = 0; i < col; i++){
        f128 temp = mat->data[col*r1 + i];
        mat->data[col*r1 + i] = mat->data[col*r2 + i];
        mat->data[col*r2 + i] = temp;
    }
}



b32 mat_ref(
    const matrix* eq_src, const matrix* sol_src,
    matrix* eq_dst, matrix* sol_dst
){
    b8 there_is_sol = 1;
    if (!sol_src || !sol_dst) there_is_sol = 0;
    if (eq_src->col != eq_dst->col || eq_src->row != eq_dst->row) return false;
    if (there_is_sol){
        if (sol_src->col != sol_dst->col || sol_src->row != sol_dst->row) return false;
        if (eq_src->col != sol_src->row) return false;
    }
    

    mem_arena_temp scrath = arena_scratch_get(NULL, 0);

    matrix* e = mat_create(scrath.arena, eq_src->row, eq_src->col);
    mat_copy(e, eq_src);
    
    // print_mat(e);
    u32 h = 0; // height pivot
    u32 k = 0; // column pivot
    if(there_is_sol){
        matrix* s = mat_create(scrath.arena, sol_src->row, sol_src->col);
        mat_copy(s, sol_src);
        while (h < e->row && k < e->col){
            // choose a pivot with the highest value
            u32 h_th = _argmax(e, h, k, 1);
    
            // no pivot
            if (fabsl(e->data[h_th*e->col + k]) < EPS){
                k++;
            }
            else{
                // swap row
                _swap(e, h, h_th);
                _swap(s, h, h_th);
    
                for (u32 r = h + 1; r < e->row; r++){
                    f128 f = (float) e->data[e->col*r + k]/e->data[e->col*h + k];
                    // modify the eq
                    e->data[e->col*r + k] = 0.0f;
                    for (u32 c = k + 1; c < e->col; c++){
                        e->data[e->col*r + c] -= e->data[e->col*h + c]*f;
                        e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                    }
                    // modify the sol
                    for (u32 c = 0; c < s->col; c++){
                        s->data[r*s->col + c] -= s->data[h*s->col + c]*f;
                        s->data[s->col*r + c] = (fabsl(s->data[s->col*r + c]) < EPS) ? 0.0f:s->data[s->col*r + c];
                    }
                }
                // reduce it further
                f128 pivot = e->data[e->col*h + k];
                if (!(fabsl(pivot) < EPS)){
                    e->data[e->col*h + k] = 1.0f;
                    for (u32 c = k + 1; c < e->col; c++){
                        e->data[e->col*h + c] /= pivot;
                    }
                    
                    for (u32 c = 0; c < s->col; c++){
                        s->data[s->col*h + c] /= pivot;
                    }
                }
                
                h++;
                k++;
                
            }
        }
        mat_copy(sol_dst, s);
    }
    else{
        while (h < e->row && k < e->col){
            // choose a pivot with the highest value
            u32 h_th = _argmax(e, h, k, 1);
    
            // no pivot
            if (fabsl(e->data[h_th*e->col + k] )< EPS){
                k++;
            }
            else{
                // swap row
                _swap(e, h, h_th);
    
                for (u32 r = h + 1; r < e->row; r++){
                    f128 f = (float) e->data[e->col*r + k]/e->data[e->col*h + k];
                    // modify the eq
                    e->data[e->col*r + k] = 0.0f;
                    for (u32 c = k + 1; c < e->col; c++){
                        e->data[e->col*r + c] = e->data[e->col*r + c] - e->data[e->col*h + c]*f;
                        e->data[e->col*r + c] = (fabsl(e->data[e->col*r + c]) < EPS) ? 0.0f:e->data[e->col*r + c];
                    }
                    // modify the sol
                }
                f128 pivot = e->data[e->col*h + k];
                if (!(fabsl(pivot) < EPS)){
                    e->data[e->col*h + k] = 1.0f;
                    for (u32 c = k + 1; c < e->col; c++){
                        e->data[e->col*h + c] /= pivot;
                    }
                }
                h++;
                k++;
                
            }
        }
        
    }
    mat_copy(eq_dst, e);
    
    arena_scratch_release(scrath);
    
    return true;

}

u32 mat_rank(const matrix* a){
    mem_arena_temp scratch  = arena_scratch_get(NULL, 0);
    matrix* b = mat_create(scratch.arena, a->row, a->col);
    mat_ref(a, NULL, b, NULL);

    u64 N_rank = 0;
    for (i32 r = b->row - 1; r > -1 ; r--){
        u32 zero = 0;
        for (u32 c = 0; c < b->col; c++){
            if (fabsl(b->data[r*b->col + c]) < EPS) zero++;
        }
        if (zero == b->col){
            N_rank ++;
        }
    }

    arena_scratch_release(scratch);
    
    return a->row - N_rank;
}