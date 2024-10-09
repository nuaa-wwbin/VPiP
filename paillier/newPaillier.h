#pragma once
// #include "ophelib/paillier_fast.h"
#include "ophelib/paillier_base.h"
#include "ophelib/packing.h"
using namespace ophelib;

namespace opaillierlib
{
    /**
     * Paillier private key. p^2 and p^q are not stored here,
     * they are precomputed in the Paillier class if needed.
     */
    class NewPrivateKey{// : public PrivateKey
    public:
        size_t key_size_bits;
        size_t a_bits;

        Integer p;
        Integer q;
        /**
         * Not all implementations use this
         */
        Integer a;
        //新变量
        Integer inv_a;//用于私钥加密 Enc(m) = g^m * r^n mod n^2 = (n * k_g * inv_a * m + 1) * r^n mod n^2
        Integer nkg; //用于加密. n * k_g = 2^lambda - 1 mod n^2
        Integer inv_kg; // 用于解密 g^a - 1 = (c^a - 1) / (g^a - 1) = ((c^a - 1) / n) * inv_kg mod n

        NewPrivateKey(const size_t key_size_bits_, const Integer &p_, const Integer &q_);
        NewPrivateKey(const size_t key_size_bits_, const size_t a_bits_, const Integer &p_, const Integer &q_, const Integer &a_);
        NewPrivateKey();

        // NewPrivateKey &operator=(const NewPrivateKey &priv_);
        bool operator==(const NewPrivateKey &input) const;

        const std::string to_string(const bool brief = true) const;

        NewPrivateKey &operator=(const NewPrivateKey &priv);
    };

    /**
     * Wrapper helper class
     */
    class NewKeyPair : public KeyPair{
    public:
        PublicKey pub;
        NewPrivateKey priv;

        NewKeyPair(const PublicKey &pub, const NewPrivateKey &priv);
        NewKeyPair();

        bool operator==(const NewKeyPair &input) const;

        const std::string to_string(const bool brief = true) const;
    };

    class NewPaillier : public ophelib::PaillierBase
    {
        /**
         * Basic, slow Ciphertext randomizer
         */
        class Randomizer {
        protected:
            const NewPaillier *paillier;
            Integer g_pow_n;
            const Integer r() const;
            bool precomputed = false;

        public:
            Randomizer(const NewPaillier *paillier);

            virtual void precompute();
            /**
             * Get a random value to randomize the ciphertext with
             */
            virtual const Integer get_noise() const;
            virtual const std::string to_string(const bool brief = true) const;
        };

        /**
         * Fast randomizer using a random cache/lookup table
         */
        class FastRandomizer: Randomizer {
            /**
             * Size of lookup table, i.e.
             * number of random `(g^n)^r` to generate
             */
            const size_t r_lut_size;

            /**
             * How many values to select randomly from the
             * lookup table each time get_noise is called
             */
            const size_t r_use_count;

            /**
             * Lookup table
             */
            std::vector<Integer> gn_pow_r;
        public:
            FastRandomizer(const NewPaillier *paillier, const size_t r_lut_size, const size_t r_use_count);

            /**
             * Fill random cache
             */
            void precompute();
            const Integer get_noise() const;
            const std::string to_string(const bool brief = true) const;
        };
        
    private:
        const size_t a_bits;
        const size_t r_bits;
        FastRandomizer randomizer;

        /**
         * pubkey precomputation
         */
        Integer n2;

        /**
         * privkey precomputation
         */
        Integer mu;

        Ciphertext precomputed_zero;
    public:
        NewPrivateKey priv;
        PublicKey pub;
    protected:
        NewPaillier() = delete;
        NewPaillier(const size_t key_size_bits_, const size_t a_bits_, const size_t r_bits_);

        void check_valid_key_size(size_t key_size_bits) const;
        void check_valid_r_bits(size_t r_bits) const;
        size_t param_a_bits(size_t key_size_bits) const;
        size_t param_r_bits(size_t key_size_bits) const;
        size_t param_r_lut_size(size_t r_bits) const;
        size_t param_r_use_count(size_t r_bits) const;
        void precompute();
        
        Integer check_plaintext(const Integer &plaintext) const;
        
    public:
        NewPaillier(const size_t key_size_bits);
        NewPaillier(const PublicKey &pub);
        NewPaillier(const PublicKey &pub, const NewPrivateKey &priv);
        NewPaillier(const NewKeyPair &pair);
        ~NewPaillier();
        /**
         * @brief Get the pack Count of each ciphertext 
         * 
         * @param plaintext_bits 
         * @return size_t 
         */
        size_t get_packCount(const size_t plaintext_bits);

        /**
         * @brief generate_keys_fast 相比generate_keys函数，本函数的g=n+1；
         * 这可以提升加密过程和同态计算。(n+1)^m mod n^2 = nm+1 mod n^2
         */
        void generate_keys();

        Integer decrypt(const Ciphertext &ciphertext) const;
        Ciphertext encrypt(const Integer &plaintext) const;
        Ciphertext zero_ciphertext();

        const std::string to_string(const bool brief = true) const final;

        ophelib::PackedCiphertext encrypt_pack(const Integer *plaintexts_begin, const Integer *plaintexts_end, const size_t plaintext_bits);
        ophelib::PackedCiphertext encrypt_pack(const Vec<Integer> &plaintexts, size_t plaintext_bits);
        ophelib::PackedCiphertext encrypt_halfPack(const Vec<Integer> &plaintexts, size_t plaintext_bits);
        
        void decrypt_pack(const PackedCiphertext &ciphertext, Integer *plaintexts_begin, Integer *plaintexts_end);
        Vec<Integer> decrypt_pack(const PackedCiphertext &ciphertext);


        Integer decrypt_test(const Ciphertext &ciphertext) const;
    };
}