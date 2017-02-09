// ===================================================================
// @last-modified 2016-06-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <tmmintrin.h> // Shuffle bytes
#include <wmmintrin.h> // AES NI

#include "kiasu.h"

// ---------------------------------------------------------------------
// Load, Store, Helpers
// ---------------------------------------------------------------------

#define loadu(p)       _mm_loadu_si128((__m128i*)p)
#define load(p)        _mm_load_si128((__m128i*)p)
#define storeu(p,x)    _mm_storeu_si128((__m128i*)p, x)
#define store(p,x)     _mm_store_si128((__m128i*)p, x)
#define vxor(x,y)      _mm_xor_si128(x,y)
#define vxor3(x,y,z)   _mm_xor_si128(_mm_xor_si128(x,y),z)

#define zero           _mm_setzero_si128()
#define one            _mm_set_epi32(0, 0, 0, 1)

// ---------------------------------------------------------------------

#define vaesdec(x,k)       _mm_aesdec_si128(x, k)
#define vaesdeclast(x,k)   _mm_aesdeclast_si128(x, k)
#define vaesenc(x,k)       _mm_aesenc_si128(x, k)
#define vaesenclast(x,k)   _mm_aesenclast_si128(x, k)
#define vaesimc(x)         _mm_aesimc_si128(x)

// ---------------------------------------------------------------------

#ifdef DEBUG
static void print_hex(const char label[], 
                      const uint8_t* x, 
                      const size_t n) {
    printf("%s ", label);

    for (size_t i = 0; i < n; ++i) {
        printf("%02x", x[i]);
    }

    puts("");
}

// ---------------------------------------------------------------------

static void print_block(const char label[], const uint8_t* x) {
    print_hex(label, x, KIASU_BLOCKLEN);
}
#endif

// ---------------------------------------------------------------------

static void mixcolumns(block out, const block in) {
    __m128i state = load(in);
    state = _mm_aesdeclast_si128(state, zero);
    state = _mm_aesenc_si128(state, zero);
    storeu(out, state);
}

// ---------------------------------------------------------------------

void kiasu_encrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key) {
    __m128i state = load(in);
    state = vaesenc(state, vxor(load(round_key), load(tweak)));
    storeu(out, state);
}

// ---------------------------------------------------------------------

void kiasu_encrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {
    const __m128i t = load(tweak); 
    const __m128i* k = ctx->enc_key;
    __m128i state = load(in);

    // If #rounds != KIASU_ROUNDS, we start directly from the colliding round
    // where we do not want to break the collision from XORing the tweak to the
    // state before the first round.

    if (num_rounds == KIASU_ROUNDS) {
        state = vxor3(state, t, k[0]);
        state = vaesenc(state, vxor(k[1], t));
    } else {
        state = vxor(state, k[1]);
    }

    for (uint8_t i = 2; i <= KIASU_ROUNDS; ++i) {
        state = vaesenc(state, vxor(k[i], t));
    }

    storeu(out, state);
}

// ---------------------------------------------------------------------

void kiasu_decrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key) {
    __m128i state = load(in);
    state = vaesdec(state, vxor(load(round_key), load(tweak)));
    storeu(out, state);
}

// ---------------------------------------------------------------------

void kiasu_decrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {

    const uint8_t first_round = KIASU_ROUNDS - num_rounds;
    const __m128i t = load(tweak); 
    const __m128i t_imc = vaesimc(t);
    const __m128i* k = ctx->dec_key;

    __m128i state = vxor3(load(in), k[KIASU_ROUNDS], t);
    state = vaesimc(state);

    for (uint8_t i = KIASU_ROUNDS-1; i > first_round; --i) {
        state = vaesdec(state, vxor(k[i], t_imc));
    }

    state = vaesdeclast(state, vxor(ctx->enc_key[first_round], t));
    storeu(out, state);
}

// ---------------------------------------------------------------------
// AES
// ---------------------------------------------------------------------

static __m128i aes_keygen_assist(__m128i temp1, __m128i temp2)
{
    __m128i temp3;
    temp2 = _mm_shuffle_epi32(temp2, 0xff);
    temp3 = _mm_slli_si128(temp1, 0x4);
    temp1 = vxor(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = vxor(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = vxor(temp1, temp3);
    temp1 = vxor(temp1, temp2);
    return temp1;
}

// ---------------------------------------------------------------------

#define aes_expand_round_key(round, rcon) \
    tmp = _mm_aeskeygenassist_si128(enc_key[round-1], rcon); \
    enc_key[round] = aes_keygen_assist(enc_key[round-1], tmp)

// ---------------------------------------------------------------------

static void aes_revert_key(const __m128i* enc_key, __m128i* dec_key)
{
    dec_key[0] = enc_key[0];
    dec_key[1] = _mm_aesimc_si128(enc_key[1]);
    dec_key[2] = _mm_aesimc_si128(enc_key[2]);
    dec_key[3] = _mm_aesimc_si128(enc_key[3]);
    dec_key[4] = _mm_aesimc_si128(enc_key[4]);
    dec_key[5] = _mm_aesimc_si128(enc_key[5]);
    dec_key[6] = enc_key[6];
}

// ---------------------------------------------------------------------

static void aes_expand_key(__m128i* enc_key, const block key) {
    __m128i tmp;

    enc_key[0] = load(key);
    aes_expand_round_key(1, 0x01);
    aes_expand_round_key(2, 0x02);
    aes_expand_round_key(3, 0x04);
    aes_expand_round_key(4, 0x08);
    aes_expand_round_key(5, 0x10);
    aes_expand_round_key(6, 0x20);

    // We only store at most 6 round keys
    // aes_expand_round_key(7, 0x40);
    // aes_expand_round_key(8, 0x80);
    // aes_expand_round_key(9, 0x1b);
    // aes_expand_round_key(10, 0x36);
}

// ---------------------------------------------------------------------

void kiasu_key_setup(kiasu_ctx* ctx, const block key) {
    aes_expand_key(ctx->enc_key, key);
    aes_revert_key(ctx->enc_key, ctx->dec_key);
}

// ---------------------------------------------------------------------

static const uint8_t TWEAK_SHIFTS_AFTER_SBOX[KIASU_SBOXLEN] = {
    0x00, 0x1f, 0x14, 0x18, 0x91, 0x08, 0x0c, 0xa6, 0x53, 0x62, 0x04, 0x48, 0x9d, 0xb4, 0xc8, 0x15, 
    0xa9, 0xe1, 0xaa, 0x1e, 0x99, 0x3a, 0x24, 0x93, 0xce, 0xb7, 0xc1, 0xcc, 0xff, 0xc7, 0x11, 0xa3, 
    0xd4, 0x9e, 0xf0, 0x45, 0x55, 0x5c, 0x94, 0xaf, 0x57, 0xc6, 0x86, 0x92, 0x12, 0xbb, 0x52, 0x76, 
    0x67, 0xa4, 0x40, 0xa0, 0x7b, 0xf5, 0x66, 0xf9, 0x64, 0x71, 0xe3, 0x81, 0x88, 0x44, 0xd1, 0x16, 
    0x6a, 0xe0, 0x4f, 0x79, 0x78, 0x0d, 0x39, 0xc3, 0x31, 0x58, 0xb5, 0xd0, 0x4a, 0x80, 0x4c, 0xe7, 
    0x30, 0xb2, 0x63, 0x8e, 0x43, 0x9f, 0xd2, 0x38, 0x09, 0xa8, 0xdd, 0x5a, 0x29, 0x2f, 0x3b, 0xac, 
    0xb3, 0x8c, 0xc9, 0x98, 0x20, 0x2e, 0x50, 0xe6, 0x26, 0x9a, 0x61, 0x1c, 0x33, 0x5f, 0xfc, 0xcb, 
    0x32, 0xc0, 0x23, 0xec, 0xf1, 0xfe, 0x5b, 0x96, 0xdf, 0xd5, 0xb9, 0x42, 0x73, 0x9c, 0x90, 0xb1, 
    0xae, 0x6f, 0x70, 0x8f, 0x3c, 0xf4, 0x27, 0x74, 0xa7, 0xc4, 0x1d, 0x5e, 0x07, 0x3e, 0x7a, 0x10, 
    0x03, 0xe2, 0x2c, 0xbf, 0x41, 0x49, 0xf3, 0xeb, 0x25, 0x8d, 0xdb, 0x77, 0xbd, 0x3d, 0x68, 0xb8, 
    0x83, 0x51, 0x59, 0x69, 0x2a, 0x65, 0x47, 0x3f, 0xa1, 0xb0, 0xcf, 0x01, 0xf2, 0xf6, 0x87, 0x1a, 
    0x84, 0xab, 0x54, 0x0e, 0xee, 0xb6, 0x2d, 0xca, 0x0f, 0x35, 0x97, 0x89, 0x06, 0x19, 0xcd, 0x6b, 
    0xd9, 0x1b, 0x46, 0x4d, 0x7f, 0xc5, 0xd7, 0xa5, 0x8b, 0xbe, 0x17, 0x7c, 0x28, 0xde, 0xe8, 0xe9, 
    0x13, 0x5d, 0xd6, 0x05, 0x2b, 0x60, 0x95, 0x6d, 0x02, 0x56, 0x34, 0xda, 0xe5, 0xa2, 0x7e, 0xfd, 
    0x82, 0x9b, 0xfb, 0x72, 0x0a, 0xba, 0xed, 0xf7, 0xf8, 0x7d, 0xe4, 0x8a, 0xad, 0x36, 0x4b, 0xbc, 
    0xef, 0xc2, 0xea, 0x6e, 0xdc, 0x85, 0x21, 0x0b, 0x22, 0xfa, 0x4e, 0x6c, 0xd3, 0x37, 0xd8, 0x75
};

// ---------------------------------------------------------------------

void prepare_delta(block out, const uint8_t tweak_shift) {
    const uint8_t tweak_shift_after_sbox = TWEAK_SHIFTS_AFTER_SBOX[tweak_shift];

    block DELTA_IN = { 
        tweak_shift_after_sbox, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 
    };
    mixcolumns(out, DELTA_IN);
    out[0] ^= tweak_shift;
}

// ---------------------------------------------------------------------
// Reserved for later key recovery
// ---------------------------------------------------------------------

// static void invert_mixcolumns(block out, const block in) {
//     __m128i state = load(in);
//     // state = vaesenclast(state, zero);
//     // state = vaesdec(state, zero);
//     state = vaesimc(state);
//     storeu(out, state);
// }

// ---------------------------------------------------------------------

// static void invert_shiftrows(block out, const block in) {
//     const __m128i shuffle = _mm_setr_epi8(
//          0, 13, 10,  7,  4,  1, 14, 11,  8,  5,  2, 15, 12,  9,  6,  3
//     );
//     // 0x03, 0x06, 0x09, 0x12, 
//     // 0x15, 0x02, 0x05, 0x08, 
//     // 0x11, 0x14, 0x01, 0x04, 
//     // 0x07, 0x10, 0x13, 0x00
//     __m128i state = load(in);
//     state = _mm_shuffle_epi8(state, shuffle);
//     storeu(out, state);
// }
