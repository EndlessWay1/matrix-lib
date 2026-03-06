
matrix_cplx *mat_c_create(mem_arena* arena, u32 row, u32 col){
    matrix_cplx* c = PUSH_STRUCT(arena, matrix_cplx);
    c->col = col;
    c->row = row;
    return c;
}

b32 mat_c_to_r(matrix* dst, const matrix_cplx* src){
    if (dst->col != src->col || dst->row != src->row) return false;
    u64 size = (u64)dst->col*dst->row;
    for (u64 i = 0; i < size; i++){
        dst->data[i] = src->data[i].v[0];
    }
    return true;
}

b32 mat_r_to_c(matrix_cplx* dst, const matrix* src){
    if (dst->col != src->col || dst->row != src->row) return false;
    u64 size = (u64)dst->col*dst->row;
    for (u64 i = 0; i < size; i++){
        dst->data[i].v[0] = src->data[i];
    }
    return true;
}


void print_matc(const matrix_cplx* mat){
    for(u64 i= 0; i < mat->row; i++){
        printf("[ %Lf %Lfj", mat->data[i*mat->col].r, mat->data[i*mat->col].i);
        for(u64 j = 1; j < mat->col; j++){
            printf(", %Lf %Lfj", mat->data[i*mat->col].r, mat->data[i*mat->col].i);
        }
        printf("]\n");
    }
}

void mat_clearc(matrix_cplx* mat){
    memset(mat->data, 0, sizeof(complex_data)*(u64) mat->col*mat->row);
}

b32 mat_copyc(matrix_cplx* dst, const matrix_cplx* src){
    if (dst->col != src->col || dst->row != src->row) return false;
    memcpy(dst->data, src->data, sizeof(complex_data)*(u64) src->col*dst->row);
    return true;
}

void mat_fillc(matrix_cplx* mat, f128 r, f128 i){
    u64 size = (u64)mat->col*mat->row;
    for (u64 c = 0; c < size; c++){
        mat->data[c].v[0] = r;
        mat->data[c].v[1] = i;
    }
}

void mat_fillc_rand(matrix_cplx* mat, f128 lower_r, f128 upper_r, f128 lower_i, f128 upper_i){
    u64 size = (u64)mat->col*mat->row;
    f128 l_r = MIN(upper_r, lower_r);
    f128 l_i = MIN(upper_i, lower_i);
    f128 u_r = MAX(upper_r, lower_r);
    f128 u_i = MAX(upper_i, lower_i);
    for (u64 i = 0; i < size; i++){
        mat->data[i].v[0] = l_r + (f128)(rand()/(f128)(RAND_MAX/(u_r - l_r)));
        mat->data[i].v[1] = l_i + (f128)(rand()/(f128)(RAND_MAX/(u_i - l_i)));
    }
}

u64 mat_argmaxc(matrix_cplx* mat, b8 use_i){
    use_i = (use_i)? 1 : 0;
    u64 size = (u64) mat->col* mat->row;
    u64 max_i = 0;
    for (u64 i = 0; i < size; i++){
        if (mat->data[i].v[use_i] > mat->data[max_i].v[use_i]){
            max_i = i;
        }
    }
    return max_i;
}

complex_data mat_sumc(matrix_cplx* mat){
    
    complex_data sums = (complex_data){
        .i = 0.0f,
        .r = 0.0f,
    };

    u64 size = (u64) mat->col* mat->row;
    for (u64 i = 0; i < size; i++){
        sums.i += mat->data[i].v[1];
        sums.r += mat->data[i].v[0];
    }
    return sums;
}

b32 mat_addc(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    if (a->row != b->row || a->col != b->col) return false;
    if (a->row != out->row || a->col != out->col) return false;

    u64 size = (u64) a->col*a->row;
    for(u64 i = 0; i < size; i++){
        out->data[i].v[0] = a->data[i].v[0] + b->data[i].v[0];
        out->data[i].v[1] = a->data[i].v[1] + b->data[i].v[1];
    }

    return true;
}

b32 mat_subc(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    if (a->row != b->row || a->col != b->col) return false;
    if (a->row != out->row || a->col != out->col) return false;

    u64 size = (u64) a->col*a->row;
    for(u64 i = 0; i < size; i++){
        out->data[i].v[0] = a->data[i].v[0] - b->data[i].v[0];
        out->data[i].v[1] = a->data[i].v[1] - b->data[i].v[1];
    }

    return true;
}

void mat_scalec(matrix_cplx* mat, f128 scale){
    
    u64 size = (u64) mat->col*mat->row;
    for(u64 i = 0; i < size; i++){
        mat->data[i].v[0] *= scale;
        mat->data[i].v[1] *= scale;
    }

}


#define CMPLX_SUB(o, a, b) \
    (o).i += ((a).i - (b).i);\
    (o).r += ((a).r - (b).r)\

#define CMPLX_ADD(o, a, b) \
    (o).i += ((a).i + (b).i);\
    (o).r += ((a).r + (b).r)\

#define CMPLX_MUL(o, a, b) \
    (o).r += ((((a).r)*((b).r)) - (((a).i)*((b).i)));\
    (o).i += ((((a).r)*((b).i)) + (((b).r)*((a).i)))\


void _mat_mulc_nn(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 k = 0; k < a->col; k++){
            for (u64 j = 0; j < out->col; j++){
                CMPLX_MUL(out->data[i*out->col + j], a->data[i*a->col + k], b->data[k*out->col + j]);
            }
        }
    }
}
void _mat_mulc_tn(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    for (u64 k = 0; k < a->col; k++){
        for (u64 i = 0; i < out->row; i++){
            for (u64 j = 0; j < out->col; j++){
                CMPLX_MUL(out->data[i*out->col + j], a->data[k*a->col + i], b->data[k*out->col + j]);
            }
        }
    }
}

void _mat_mulc_nt(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 j = 0; j < out->col; j++){
            for (u64 k = 0; k < a->col; k++){
                CMPLX_MUL(out->data[i*out->col + j], a->data[i*a->col + k], b->data[j*out->col + k]);
            }
        }
    }
}
void _mat_mulc_tt(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 j = 0; j < out->col; j++){
            for (u64 k = 0; k < a->col; k++){
                CMPLX_MUL(out->data[i*out->col + j], a->data[k*a->col + i], b->data[j*out->col + k]);
            }
        }
    }
}


b32 mat_mulc(
    matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b, 
    b8 zero_out, b8 transpose_a, b8 transpose_b
){
    mem_arena_temp scratch = arena_scratch_get(NULL, 0);
    matrix_cplx* copy_a = mat_c_create(scratch.arena, a->row, a->col);
    matrix_cplx* copy_b = mat_c_create(scratch.arena, b->row, b->col);

    mat_copyc(copy_a, a);
    mat_copyc(copy_b, b);
    a = copy_a;
    b = copy_b;

    u32 a_row = transpose_a ? a->col : a->row;
    u32 a_col = transpose_a ? a->row : a->col;
    u32 b_row = transpose_b ? b->col : b->row;
    u32 b_col = transpose_b ? b->row : b->col;
    if (a_col != b_row || out->row != a_row || out->col != b_col) return false;
    

    if (zero_out){
        mat_clearc(out);
    }

    u32 transpose = (transpose_a << 1) | transpose_b;

    switch (transpose){
        case 0b00: { _mat_mulc_nn(out, a, b); } break; 
        case 0b10: { _mat_mulc_tn(out, a, b); } break; 
        case 0b01: { _mat_mulc_nt(out, a, b); } break; 
        case 0b11: { _mat_mulc_tt(out, a, b); } break; 
    }

    arena_scratch_release(scratch);
    
    return true;    
}

