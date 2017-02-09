#pragma once

#include <stdint.h>
#include <emmintrin.h>

// ---------------------------------------------------------------------

#define KIASU_SBOXLEN  256
#define KIASU_BLOCKLEN  16
#define KIASU_ROUNDS     6
#define KIASU_KEYS       7
#define KIASU_KEYLEN    16

typedef uint8_t block[KIASU_BLOCKLEN];

typedef struct {
    __m128i enc_key[KIASU_KEYS];
    __m128i dec_key[KIASU_KEYS];
} kiasu_ctx;

// ---------------------------------------------------------------------

void kiasu_encrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key);

// ---------------------------------------------------------------------

void kiasu_encrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds);

// ---------------------------------------------------------------------

void kiasu_decrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key);

// ---------------------------------------------------------------------

void kiasu_decrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds);

// ---------------------------------------------------------------------

void kiasu_key_setup(kiasu_ctx* ctx, const block key);

// ---------------------------------------------------------------------

void prepare_delta(block out, const uint8_t tweak_shift);
