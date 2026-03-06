typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef i8 b8;
typedef i32 b32;
typedef float f32;
typedef double f64;
typedef long double f128;

typedef enum {
    TY_STRING,
    TY_CHAR,
    TY_I8,
    TY_I16,
    TY_I32,
    TY_I64,
    TY_U8,
    TY_U16,
    TY_U32,
    TY_U64,
    TY_F32,
    TY_B8,
    TY_B32,
    TY_F64,
    TY_F128
} dtypes;

u32 size_data_types(dtypes d){
    switch (d) {
        case TY_STRING:{return sizeof(char);}
        case TY_CHAR:{return sizeof(char);}
        case TY_I8:{return sizeof(i8);}
        case TY_I16:{return sizeof(i16);}
        case TY_I32: {return sizeof(i32);}
        case TY_I64: {return sizeof(i64);}
        case TY_U8:{return sizeof(u8);}
        case TY_U16:{return sizeof(u16);}
        case TY_U32: {return sizeof(u32);}
        case TY_U64: {return sizeof(u64);}
        case TY_F32: {return sizeof(f32);}
        case TY_B8: {return sizeof(b8);}
        case TY_B32: {return sizeof(b32);}
        case TY_F64: {return sizeof(f64);}
        default:
            return false;
    }
}

u64 _get_size(dtypes d, void* src){
    if (!src) return 0;
    
    if (d == TY_STRING){ 
        return (strlen((char*)src) + 1) * sizeof(char);
    }
    return size_data_types(d);
}


#define TO_TYPE(ptr, type) \
(\
    (type == TY_STRING) ? (char*) ptr:\
    (type == TY_CHAR) ? *(char *) ptr:\
    (type == TY_I8) ? *(i8*) ptr:\
    (type == TY_I16) ? *(i16*) ptr:\
    (type == TY_I32) ? *(i32*) ptr:\
    (type == TY_I64) ? *(i64*) ptr:\
    (type == TY_U8) ? *(u8*) ptr:\
    (type == TY_U16) ? *(u16*) ptr:\
    (type == TY_U32) ? *(u32*) ptr:\
    (type == TY_U64) ? *(u64*) ptr:\
    (type == TY_F32) ? *(f32*) ptr:\
    (type == TY_B8) ? *(b8*) ptr:\
    (type == TY_B32) ? *(b32*) ptr:\
    (type == TY_F64) ? *(f64*) ptr:\
    (ptr)\
)\

