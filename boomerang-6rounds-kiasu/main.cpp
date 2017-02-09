// ===================================================================
// @last-modified 2016-09-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "kiasu.h"

// ---------------------------------------------------------------------
// Static Vars and Typedefs
// ---------------------------------------------------------------------

#ifdef DEBUG
static const block BASE_KEY = { 
    0x7d, 0x89, 0x55, 0x51, 0xa4, 0x5f, 0xa7, 0xce, 
    0x80, 0x25, 0x56, 0xe5, 0xff, 0x76, 0xf1, 0xcf
};
#endif
static const block BASE_PLAINTEXT = { 
    0x20, 0xa7, 0xd3, 0x35, 0x5d, 0x4a, 0x4e, 0xa4, 
    0x0b, 0xfc, 0x8f, 0x5a, 0x35, 0xe7, 0x8f, 0x2d
};
static const block BASE_TWEAK = { 
    0xe3, 0x23, 0x00, 0x00, 0xe0, 0x11, 0x00, 0x00, 
    0x84, 0x6d, 0x00, 0x00, 0x6f, 0x87, 0x00, 0x00
};
// The delta difference at the end of our boomerang trails.
static const uint64_t MASKS[7][2] = {
    { 0L, 0L }, 
    { 0L, 0L }, 
    { 0L, 0L }, 
    { 0x00000000FFFFFFFFL, 0xFFFFFFFFFFFFFFFFL }, 
    { 0x00FFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL }, 
    { 0x00FFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL }, 
    { 0x00FFFFFFFF00FFFFL, 0xFFFF00FFFFFFFF00L }
};
static const uint8_t TWEAK_SHIFT = 0x01;
static const uint8_t NUM_QUARTETS = 16;
static block DELTA;

typedef struct {
    block p;
    block q;
    uint64_t q1_64;
    uint64_t q2_64;
    block t;
    block c;
    block d;
} pair;

// ---------------------------------------------------------------------

typedef struct {
    block p;
    block p_;
    block q;
    block q_;
    block t;
    block t_;
    block c;
    block c_;
    block d;
    block d_;
} quartet;

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

static void xor_block(block out, 
                      const block a, 
                      const block b) {
    for (uint32_t i = 0; i < KIASU_BLOCKLEN; ++i) {
        out[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

static void get_random(void* buffer, size_t num_bytes) {
    FILE* fp;
    fp = fopen("/dev/urandom", "r");
    fread(buffer, 1, num_bytes, fp);
    fclose(fp);
}

// ---------------------------------------------------------------------

static inline uint64_t to_uint64(const uint8_t* x) {
    return ((uint64_t)x[0] << 56) 
         | ((uint64_t)x[1] << 48)
         | ((uint64_t)x[2] << 40)
         | ((uint64_t)x[3] << 32)
         | ((uint64_t)x[4] << 24)
         | ((uint64_t)x[5] << 16)
         | ((uint64_t)x[6] <<  8)
         |  (uint64_t)x[7];
}

// ---------------------------------------------------------------------

static inline int compare_plaintexts(const void* lhs, 
                                     const void* rhs) {
    // Return -1 if B < A
    // Return  0 if B = A
    // Return +1 if B > A
    // Compares texts for the 48 bits that are NOT on the main diagonal.
    // It actually does not matter how we sort since our interest limits to 
    // the identification of pairs which are equal in all 48 bits NOT on the
    // main diagonal for finding quartets.
    const pair* a_ = (const pair*)lhs;
    const pair* b_ = (const pair*)rhs;
    if(a_->q1_64 > b_->q1_64)  return(+1);
    if(a_->q1_64 < b_->q1_64)  return(-1);
    if(a_->q2_64 > b_->q2_64)  return(+1);
    if(a_->q2_64 < b_->q2_64)  return(-1);
    return (0);
}

// ---------------------------------------------------------------------

static inline int are_equal(const pair* a, const pair* b) {
    return (a->q1_64 == b->q1_64) && (a->q2_64 == b->q2_64);
}

// ---------------------------------------------------------------------

static void create_plaintext(block p, 
                             const uint32_t i, 
                             const uint8_t num_rounds) {
    if (num_rounds == KIASU_ROUNDS) {
        p[ 0] = ((i >> 24) & 0xFF);
        p[ 5] = ((i >> 16) & 0xFF); 
        p[10] = ((i >>  8) & 0xFF);
        p[15] = i & 0xFF;
    } else {
        p[ 0] = ((i >> 24) & 0xFF);
        p[ 1] = ((i >> 16) & 0xFF); 
        p[ 2] = ((i >>  8) & 0xFF);
        p[ 3] = i & 0xFF;
    }
}

// ---------------------------------------------------------------------

static void create_pair(pair* target, 
                        const block p, 
                        const block c, 
                        const block q, 
                        const block d, 
                        const block t, 
                        const uint64_t mask1, 
                        const uint64_t mask2) {
    memcpy(target->p, p, KIASU_BLOCKLEN);
    memcpy(target->c, c, KIASU_BLOCKLEN);
    memcpy(target->q, q, KIASU_BLOCKLEN);
    memcpy(target->d, d, KIASU_BLOCKLEN);
    memcpy(target->t, t, KIASU_BLOCKLEN);
    target->q1_64 = to_uint64(q) & mask1;
    target->q2_64 = to_uint64(q+8) & mask2;
}

// ---------------------------------------------------------------------

static void print_quartet(const pair* lhs, const pair* rhs) {
    print_block("P ", lhs->p);
    print_block("Q ", lhs->q);
    print_block("T ", lhs->t);
    print_block("C ", lhs->c);
    print_block("D ", lhs->d);
    
    print_block("P'", rhs->p);
    print_block("Q'", rhs->q);
    print_block("T'", rhs->t);
    print_block("C ", rhs->c);
    print_block("D ", rhs->d);
    puts("");
}

// ---------------------------------------------------------------------

static void find_num_quartets(pair* pairs, 
                              const uint32_t num_pairs, 
                              uint32_t* num_quartets) {
    uint32_t num_equal = 0;
    pair* current;
    pair* previous;

    for (uint32_t i = 1; i < num_pairs; ++i) {
        previous = &(pairs[i-1]);
        current = &(pairs[i]);

        if (are_equal(current, previous)) {
            num_equal++;

            if (num_equal + *num_quartets < NUM_QUARTETS) {
                // print_quartet(current, previous);
            }
        } else if (num_equal > 0) {
            *num_quartets += ((num_equal + 1) * num_equal) / 2;
            num_equal = 0;
        }
    }

    if (num_equal > 0) {
        *num_quartets += ((num_equal + 1) * num_equal) / 2;
        num_equal = 0;
    }
}

// ---------------------------------------------------------------------

static void find_quartets(pair* pairs, 
                          const uint32_t num_pairs, 
                          uint32_t* num_quartets) {
    qsort(pairs, num_pairs, sizeof(pair), compare_plaintexts);
    find_num_quartets(pairs, num_pairs, num_quartets);
}

// ---------------------------------------------------------------------

static void get_key(block key) {
#ifdef DEBUG
    memcpy(key, BASE_KEY, KIASU_KEYLEN);
#else
    // Initialize Kiasu-BC with a random key.
    get_random(key, KIASU_KEYLEN);
#endif
}

// ---------------------------------------------------------------------

static void perform_experiment(block key, 
                               const uint32_t num_sets, 
                               const uint32_t num_texts_per_set, 
                               const uint32_t num_structures,
                               const uint8_t num_rounds,
                               uint32_t* num_quartets) {
    kiasu_ctx ctx;
    kiasu_key_setup(&ctx, key);

    block t;
    block shifted_tweak;
    memcpy(t, BASE_TWEAK, KIASU_BLOCKLEN);
    memcpy(shifted_tweak, BASE_TWEAK, KIASU_BLOCKLEN);

    const uint32_t num_pairs = num_sets * num_texts_per_set;
    pair* pairs = (pair*)malloc(num_pairs * sizeof(pair));
    const uint64_t mask1 = MASKS[num_rounds][0];
    const uint64_t mask2 = MASKS[num_rounds][1];

    puts("------------------");
    print_block("Base K ", key);
    print_block("Base P ", BASE_PLAINTEXT);
    print_block("Base T ", t);
    print_block("DELTA  ", DELTA);

    // Build texts
    for (uint32_t k = 0; k < num_structures; ++k) {
        // ---------------------------------------------------------------
        // Setup base plaintext: random for each structure
        // ---------------------------------------------------------------

        block p;
        block q;
        block c;
        block d;
        get_random(p, KIASU_BLOCKLEN);
        
        for (uint32_t i = 0; i < num_texts_per_set; ++i) {
            // Encode the counter into the first column.
            create_plaintext(p, i, num_rounds);

            // Encrypt P under all considered tweaks
            for(uint32_t j = 0; j < num_sets; ++j) {
                t[0] = j & 0xFF;
                shifted_tweak[0] = t[0] ^ TWEAK_SHIFT;

                kiasu_encrypt(&ctx, c, p, t, num_rounds);
                xor_block(d, c, DELTA);
                kiasu_decrypt(&ctx, q, d, shifted_tweak, num_rounds);

                pair target;
                create_pair(&target, p, c, q, d, t, mask1, mask2);
                pairs[i * num_sets + j] = target;
#ifdef DEBUG
                print_block("p ", p);
                print_block("q ", q);
                printf("q_ %016lx%016lx\n", target.q1_64, target.q2_64);
                printf("m_ %016lx%016lx\n", mask1, mask2);
                print_block("c ", c);
                print_block("d ", d);
                print_block("t ", t);
                print_block("t'", shifted_tweak);
                puts("------------------");
#endif  
            }
        }

        // ---------------------------------------------------------------
        // Correct quartets must be in the same structure. 
        // For num_rounds < KIASU_ROUNDS, we could also perform quartet 
        // search for each plaintext separately (since for those, we start 
        // from the end of Round 1 where P = P'). 
        // For num_rounds = KIASU_ROUNDS, we start from begin of Round 1.
        // ---------------------------------------------------------------

        find_quartets(pairs, num_pairs, num_quartets);

        printf("Structure %2u/%2u: %u quartets \n", 
            k+1, num_structures, *num_quartets);
    }

    free(pairs);
}

// ---------------------------------------------------------------------

static void perform_experiments(const uint32_t num_keys, 
                                const uint32_t num_texts_per_set, 
                                const uint32_t num_sets, 
                                const uint32_t num_structures,
                                const uint8_t num_rounds) {
    block key;
    
    for (uint32_t i = 0; i < num_keys; ++i) {
        get_key(key);
        uint32_t num_quartets = 0;
        perform_experiment(
            key, num_sets, num_texts_per_set, 
            num_structures, num_rounds, &num_quartets
        );
        puts("Experiment finished");
        print_block("Base K ", key);
    }
}

// ---------------------------------------------------------------------
// User interface, main, and options.
// ---------------------------------------------------------------------

static void print_usage(const char* name) {
    puts("Tests the 6-round boomerang distinguisher we have for Kiasu-BC "\
         "with structures of 2^{p} * 2^{t} plaintexts each.");
    printf("Usage: %s <#keys> <#texts per set> <#sets> [<#rounds>] [<#structures>]\n", name);
    puts("- <#keys>: #Random keys to test. E.g. 1 or 20");
    puts("- <#texts per set>: #Texts per set. E.g. 2^{16}");
    puts("- <#sets>: #Sets. Should be <= 2^8");
    puts("- <#rounds>: #Rounds of the distinguisher. Should be 1-6. Default: 5");
    puts("- <#structures>: #Structures. Default: 1");
}

// ---------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        exit(1);
    }

    // for (uint32_t i = 0x01; i < 256; ++i) {
        // TWEAK_SHIFT = i;
    prepare_delta(DELTA, TWEAK_SHIFT);
    
    printf(     "Delta T %02x\n", TWEAK_SHIFT);

    // ---------------------------------------------------------------------
    // Interface
    // ---------------------------------------------------------------------

    const uint32_t num_keys = atoi(argv[1]);
    const uint32_t num_texts_per_set = atoi(argv[2]);
    const uint32_t num_sets = atoi(argv[3]);
    uint32_t num_structures = 1;
    uint8_t num_rounds = 5;

    if (argc >= 5) {
        num_rounds = (uint8_t)(atoi(argv[4]) & 0xFF);
    }

    if (argc >= 6) {
        num_structures = atoi(argv[5]);
    }

    perform_experiments(
        num_keys, num_texts_per_set, num_sets, num_structures, num_rounds
    );
    // }
    
    return 0;
}
