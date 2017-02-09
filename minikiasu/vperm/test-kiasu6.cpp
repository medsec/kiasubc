// ===================================================================
// @last-modified 2016-06-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "kiasu.h"

// ---------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

static void run_tests(const size_t num_texts, 
                      const block k, 
                      const block* expected_p, 
                      const block* tweaks,
                      const block* expected_c) {
    kiasu_ctx ctx;
    kiasu_key_setup(&ctx, k);
    int num_errors = 0;

    for (size_t i = 0; i < num_texts; ++i) {
        block c;
        kiasu_encrypt(&ctx, c, expected_p[i], tweaks[i], KIASU_ROUNDS);

        if (memcmp(expected_c[i], c, KIASU_BLOCKLEN)) {
            printf("Kiasu encryption test %zu failed\n", i+1);
            print_block("Expected", expected_c[i]);
            print_block("But was ", c);
            num_errors++;
        }
    }

    printf("%zu/%zu encryption tests correct\n", num_texts-num_errors, num_texts);
    num_errors = 0;

    for (size_t i = 0; i < num_texts; ++i) {
        block p;
        kiasu_decrypt(&ctx, p, expected_c[i], tweaks[i], KIASU_ROUNDS);
        
        if (memcmp(expected_p[i], p, KIASU_BLOCKLEN)) {
            printf("Kiasu decryption test %zu failed\n", i+1);
            print_block("Expected", expected_p[i]);
            print_block("But was ", p);
            num_errors++;
        }
    }

    printf("%zu/%zu decryption tests correct\n", num_texts-num_errors, num_texts);
}

// ---------------------------------------------------------------------

static void test_kiasu() {
    const block k = { 0x7d, 0x89, 0x55, 0x51, 0xa4, 0x5f, 0xa7, 0xce };

#define NUM_KIASU_TEST_TEXTS 4
    const block expected_p[NUM_KIASU_TEST_TEXTS] = { 
        { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }, 
        { 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF }, 
        { 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE }, 
        { 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0 }
    };

    const block expected_c[NUM_KIASU_TEST_TEXTS] = { 
        { 0x8b, 0x81, 0x19, 0x8f, 0x5b, 0x21, 0x1b, 0x36 }, 
        { 0x50, 0x81, 0x81, 0xe8, 0xbd, 0x8f, 0xb0, 0x59 },
        { 0x3c, 0x48, 0x55, 0x9d, 0x31, 0x6b, 0xfb, 0x55 },
        { 0xb3, 0x47, 0x16, 0xdd, 0xb5, 0x6a, 0x11, 0xc6 }
    };

    const block tweaks[NUM_KIASU_TEST_TEXTS] = { 
        { 0x3a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x3a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x3a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x3a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }
    };

    run_tests(NUM_KIASU_TEST_TEXTS, k, expected_p, tweaks, expected_c); 
#undef NUM_KIASU_TEST_TEXTS
}

// ---------------------------------------------------------------------

static void test_kiasu_with_found_quartets() {
    const block k = { 0x7d, 0x89, 0x55, 0x51, 0xa4, 0x5f, 0xa7, 0xce };

#define NUM_KIASU_TEST_TEXTS 8
    const block expected_p[NUM_KIASU_TEST_TEXTS] = { 
        { 0x0b, 0xfc, 0x80, 0x5a, 0x35, 0x57, 0x8f, 0x20 }, 
        { 0x0b, 0xfc, 0x80, 0x5a, 0x35, 0x67, 0x8f, 0x2a }, 
        { 0x0b, 0xfc, 0x80, 0x5a, 0x35, 0x57, 0x8f, 0x27 }, 
        { 0x0b, 0xfc, 0x80, 0x5a, 0x35, 0x67, 0x8f, 0x23 }, 
        { 0x97, 0xcc, 0xdd, 0xbf, 0xd5, 0x2c, 0x2d, 0xe9 }, 
        { 0x77, 0xcc, 0xdd, 0xbf, 0xd5, 0xbc, 0x2d, 0xe5 }, 
        { 0xaa, 0x3f, 0x0c, 0xf4, 0xea, 0x97, 0x05, 0x91 }, 
        { 0x4a, 0x3f, 0x0c, 0xf4, 0xea, 0xc7, 0x05, 0x9c }
    };

    const block expected_c[NUM_KIASU_TEST_TEXTS] = { 
        { 0x5a, 0x84, 0x6e, 0x0a, 0x6c, 0xc2, 0xbc, 0x26 }, 
        { 0x3f, 0xf4, 0x5b, 0x9b, 0x48, 0x72, 0xe1, 0x1f }, 
        { 0x38, 0xac, 0x7f, 0xbc, 0xd7, 0x06, 0x31, 0x8d }, 
        { 0x5e, 0xd1, 0xe1, 0xab, 0x69, 0x9c, 0xff, 0xf7 }, 
        { 0x31, 0xe0, 0x21, 0x22, 0x3a, 0xe2, 0xd5, 0x05 }, 
        { 0xc5, 0x35, 0x91, 0x2f, 0x18, 0xed, 0xfb, 0x83 }, 
        { 0xde, 0x20, 0xeb, 0x10, 0xbb, 0xa2, 0xbe, 0xaa }, 
        { 0x6b, 0x7c, 0x32, 0x86, 0x1d, 0xba, 0x81, 0x86 }, 
    };

    const block tweaks[NUM_KIASU_TEST_TEXTS] = { 
        { 0xaa, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x9a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x9a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0xaa, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0xba, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x8a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0x8a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }, 
        { 0xba, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 }
    };

    run_tests(NUM_KIASU_TEST_TEXTS, k, expected_p, tweaks, expected_c);
#undef NUM_KIASU_TEST_TEXTS
}

// ---------------------------------------------------------------------

int main() {
    test_kiasu();
    test_kiasu_with_found_quartets();
    return 0;
}