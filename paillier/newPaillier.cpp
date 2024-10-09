#include "newPaillier.h"
#include "ophelib/random.h"
#include "ophelib/omp_wrap.h"
#include <stdio.h>
using namespace std;
using namespace ophelib;
using namespace opaillierlib;

/*benchmark需要的头文件和自定义*/
#include <chrono>
#define TIMER_START(name) auto start_##name = std::chrono::steady_clock::now();
#define TIMER_STOP(name) \
    auto stop_##name = std::chrono::steady_clock::now(); \
    std::chrono::duration<double, std::milli> time_##name = stop_##name - start_##name; \
    printf(#name " time: %f ms\n", time_##name.count());

NewPaillier::NewPaillier(const size_t key_size_bits_, const size_t a_bits_, const size_t r_bits_) : 
            PaillierBase(key_size_bits_),
            a_bits(a_bits_),
            r_bits(r_bits_),
            randomizer(this, param_r_lut_size(r_bits), param_r_use_count(r_bits)){
}

NewPaillier::NewPaillier(const size_t key_size_bits_)
        : NewPaillier(key_size_bits_,
            param_a_bits(key_size_bits_),
            param_r_bits(key_size_bits_)) { }

NewPaillier::NewPaillier(const PublicKey &pub_)
        : NewPaillier(pub_.key_size_bits,
            param_a_bits(pub_.key_size_bits),
            param_r_bits(pub_.key_size_bits)) {
    pub = pub_;
    have_pub = true;
    precompute();
}

NewPaillier::NewPaillier(const PublicKey &pub_, const NewPrivateKey &priv_)
        : NewPaillier(pub_.key_size_bits,
            param_a_bits(pub_.key_size_bits),
            param_r_bits(pub_.key_size_bits)) {
    if(priv_.a == 0 || priv_.a_bits == 0)
        error_exit("invalid private key, not from a NewPaillier instance!");
    pub = pub_;
    priv = priv_;
    have_priv = true;
    have_pub = true;
    precompute();
}

NewPaillier::NewPaillier(const NewKeyPair &pair)
        : NewPaillier(pair.pub.key_size_bits,
                        param_a_bits(pair.pub.key_size_bits),
                        param_r_bits(pair.pub.key_size_bits)) {
    if(pair.priv.a == 0 || pair.priv.a_bits == 0)
        error_exit("invalid private key, not from a NewPaillier instance!");
    pub = pair.pub;
    priv = pair.priv;
    have_priv = true;
    have_pub = true;
    precompute();
}
NewPaillier::~NewPaillier()
{
}

void NewPaillier::check_valid_key_size(size_t key_size_bits_) const {
    #ifdef DEBUG
    if(key_size_bits_ == 1024)
        std::cerr << "WARNING: Key size of 1024 used, insecure for production!\n";
    #endif

    if(key_size_bits_ != 1024 &&
        key_size_bits_ != 2048 &&
        key_size_bits_ != 3072 &&
        key_size_bits_ != 4096 &&
        key_size_bits_ != 7680)
        error_exit("supported key_size_bits are: 2048, 3072, 4096, 7680!");
}

void NewPaillier::check_valid_r_bits(size_t r_bits_) const {
    if(r_bits_ != 80 &&
        r_bits_ != 112 &&
        r_bits_ != 128 &&
        r_bits_ != 140 &&
        r_bits_ != 192)
        error_exit("supported r_bits are: 112, 128, 140, 192!");
}

Integer NewPaillier::check_plaintext(const Integer &plaintext) const {
    if(!have_pub)
        error_exit("don't have a public key!");

    if(plaintext < 0) {
        return pub.n + plaintext;
    } else {
        return plaintext;
    }
}

size_t NewPaillier::param_a_bits(size_t key_size_bits_) const {
    check_valid_key_size(key_size_bits_);
    switch(key_size_bits_) {
        case 1024: return 320;
        case 2048: return 512;//原为512，改为320bit可以提高解密性能，且不影响安全性。
        case 3072: return 512;
        case 4096: return 512;
        case 7680: return 1024;
        default: return 0;
    }
}

size_t NewPaillier::param_r_bits(size_t key_size_bits_) const {
    check_valid_key_size(key_size_bits_);
    switch(key_size_bits_) {
        case 1024: return 80;
        case 2048: return 112;
        case 3072: return 128;
        case 4096: return 140;
        case 7680: return 192;
        default: return 0;
    }
}

size_t NewPaillier::param_r_lut_size(size_t r_bits_) const {
    check_valid_r_bits(r_bits_);
    switch(r_bits_) {
        case 80: return 256;
        case 112: return 4096;
        case 128: return 4096;
        case 140: return 8192;
        case 192: return 16384;
        default: return 0;
    }
}

size_t NewPaillier::param_r_use_count(size_t r_bits_) const {
    check_valid_r_bits(r_bits_);
    switch(r_bits_) {
        case 80: return 15;
        case 112: return 12;
        case 128: return 14;
        case 140: return 14;
        case 192: return 18;
        default: return 0;
    }
}

void NewPaillier::precompute() {
    if(!have_pub)
        error_exit("don't have a public key!");

    n2 = pub.n * pub.n;
    n2_shared = std::make_shared<Integer>(n2);

    if(have_priv) {
        fast_mod = std::make_shared<FastMod>(priv.p, priv.q, priv.p * priv.p, priv.q * priv.q, pub.n, n2);
        mu = Integer::L(fast_mod.get()->pow_mod_n2(pub.g, priv.a), pub.n).inv_mod_n(pub.n);
    }

    pos_neg_boundary = pub.n / 2;
    plaintxt_upper_boundary = pos_neg_boundary;
    plaintxt_lower_boundary = -pos_neg_boundary;

    randomizer.precompute();
    precomputed_zero = encrypt(0);
}

void NewPaillier::generate_keys(){
    Integer p, q, n, g, a;
        const size_t prime_size_bits = key_size_bits / 2 - a_bits;
        Random& rand = Random::instance();

        do {
            Integer cp, cq;

            cp = rand.rand_int_bits(prime_size_bits);
            if(cp.size_bits() != prime_size_bits)
                continue;
            cq = rand.rand_int_bits(prime_size_bits);
            if(cq.size_bits() != prime_size_bits)
                continue;

            a = rand.rand_prime(a_bits);

            p = a * cp + 1;
            while(!p.is_prime())
                p = p + a;

            q = a * cq + 1;
            while(!q.is_prime())
                q = q + a;

            n = p * q;
        }
        while(n.size_bits() != key_size_bits || p == q);

        if(p > q)
            swap(p, q);

        Integer lambda = Integer::lcm(p - 1, q - 1);
        g = Integer(2).pow_mod_n(lambda / a, n*n);
        // g = n + 1;

        priv = NewPrivateKey(key_size_bits, a_bits, p, q, a);
        pub = PublicKey(key_size_bits, n, g);

        have_priv = true;
        have_pub = true;
        this->precompute();
}

Ciphertext NewPaillier::encrypt(const Integer &plaintext) const {
    if(!have_pub)
        error_exit("don't have a public key!");

    Integer m = check_plaintext(plaintext),
            tmp;

    if(have_priv) {
        // tmp = fast_mod.get()->pow_mod_n2(pub.g, m) * randomizer.get_noise();
        
        // // enc(m) = g^m * r^n = (nkg * inv_a * m + 1) * r^n mod n^2
        tmp = (priv.nkg * priv.inv_a * m + 1) * randomizer.get_noise();
    } else {
        tmp = pub.g.pow_mod_n(m, n2) * randomizer.get_noise();
    }
    return Ciphertext(tmp % n2, n2_shared, fast_mod);
}

Integer NewPaillier::decrypt(const Ciphertext &ciphertext) const {
    if(!have_priv)
        error_exit("don't have a private key!");

    #ifdef DEBUG
    /* If they have the same pointer, they are the same. If not, it might
        * still be the same number, but initialized seperately. */
    if(ciphertext.n2_shared && this->n2_shared.get() != ciphertext.n2_shared.get() &&
            *(this->n2_shared.get()) != *(ciphertext.n2_shared.get()))
        error_exit("cannot decrypt a ciphertext from another n!");
    #endif

    // Integer ret = (
    //     Integer::L(
    //         fast_mod.get()->pow_mod_n2(ciphertext.data, priv.a),
    //         pub.n
    //     ) * mu
    // ) % pub.n;

    // (L(c^a - 1) / n) * inv_kg mod n
    Integer ret = ((Integer::L(fast_mod.get()->pow_mod_n2(ciphertext.data, priv.a),pub.n)) * priv.inv_kg) % pub.n;

    if(ret > pos_neg_boundary) {
        ret -= pub.n;
    }

    return ret;
}
Integer NewPaillier::decrypt_test(const Ciphertext &ciphertext) const{
    Integer ret = ((Integer::L(fast_mod.get()->pow_mod_n2(ciphertext.data, pub.n-1),pub.n)) / (Integer::L(fast_mod.get()->pow_mod_n2(pub.g, pub.n-1),pub.n))) % pub.n;

    // if(ret > pos_neg_boundary) {
    //     ret -= pub.n;
    // }

    return ret;
}
Ciphertext NewPaillier::zero_ciphertext() {
    if(!have_pub)
        error_exit("don't have a public key!");
    // return precomputed_zero;
    return Ciphertext(randomizer.get_noise(), n2_shared, fast_mod);
}

const std::string NewPaillier::to_string(bool brief) const {
    std::ostringstream o("");

    o << "<Paillier[" << key_size_bits << ", " << a_bits <<"]";
    o << " n2=" << n2.to_string(brief);
    if(have_pub) {
        o << " pub=" << pub.to_string(brief);
    } else {
        o << " have_pub=" << have_pub;
    }
    if(have_priv) {
        o << " priv=" << priv.to_string(brief);
        o << " mu=" << mu.to_string(brief);
    } else {
        o << " have_priv=" << have_priv;
    }
    o << " a_bits=" << a_bits;
    o << " r_bits=" << r_bits;
    o << " randomizer=" << randomizer.to_string(brief);
    o << ">";

    return o.str();
}


/*------------------- Pack enc and dec ---------------------*/
size_t NewPaillier::get_packCount(const size_t plaintext_bits) {
    return plaintext_size_bits() / plaintext_bits;
    // return (pai.plaintext_size_bits() / plaintext_bits) / 2;
}

PackedCiphertext NewPaillier::encrypt_pack(const Integer *plaintexts_begin, const Integer *plaintexts_end, const size_t slot_bits) {
    const size_t shift = slot_bits;
    const size_t n_plaintexts = (size_t) (plaintexts_end - plaintexts_begin);
    if(n_plaintexts > get_packCount(slot_bits))
        error_exit("trying to pack too many elements!\n");

    // TIMER_START(encode)
    //编码
    Integer sum = n_plaintexts > 0 ? *plaintexts_begin : 0;
    for (auto iter = plaintexts_begin + 1; iter < plaintexts_end; iter++) {
        if(iter->size_bits() > slot_bits)
            error_exit("plaintext size too large!\n");
        sum <<= shift;
        sum += *iter;
    }
    // TIMER_STOP(encode)
    //加密
    return PackedCiphertext(
            encrypt(sum),
            (const size_t) n_plaintexts,
            slot_bits
    );
}

PackedCiphertext NewPaillier::encrypt_pack(const Vec<Integer> &plaintexts, size_t plaintext_bits) {
    return encrypt_pack(plaintexts.begin(), plaintexts.end(), plaintext_bits);
}

PackedCiphertext NewPaillier::encrypt_halfPack(const Vec<Integer> &plaintexts, size_t slot_bits) {
    size_t pack_counts = get_packCount(slot_bits);
    if (pack_counts/2 < plaintexts.length()){
        error_exit("plaintexts is too much!");
    }
    Vec<Integer> expand;
    expand.SetLength(plaintexts.length()*2);
    size_t end = plaintexts.length();
    for (size_t i = 0; i < end; i++)
    {
        expand[i] = 0;
    }
    for (size_t i = 0; i < plaintexts.length(); i++)
    {
        expand[i+end] = plaintexts[i];
    }
    return encrypt_pack(expand.begin(), expand.end(), slot_bits);
}
void NewPaillier::decrypt_pack(const PackedCiphertext &ciphertext, Integer *plaintexts_begin, Integer *plaintexts_end) {
    const long n_plaintexts = ciphertext.n_plaintexts;
    const size_t plaintext_bits = ciphertext.plaintext_bits;
    if(n_plaintexts > get_packCount(plaintext_bits))
        error_exit("trying to unpack too many elements!\n");
    const size_t shift = plaintext_bits;

    if(n_plaintexts != (plaintexts_end - plaintexts_begin))
        error_exit("trying to unpack too many elements!\n");

    const Integer mask_plus_1 = Integer(1) << shift;
    const Integer mask = mask_plus_1 - 1;

    Integer sum = this->decrypt(ciphertext.data);

    for(auto i = n_plaintexts; i-- != 0;) {
        mpz_and(plaintexts_begin[i].get_mpz_t(), sum.get_mpz_t(), mask.get_mpz_t());
        const bool is_neg = (bool)mpz_tstbit(sum.get_mpz_t(), shift - 1);
        if(is_neg) {
            plaintexts_begin[i] -= mask_plus_1;
            sum -= plaintexts_begin[i];
        }
        sum >>= shift;
    }
}
Vec<Integer> NewPaillier::decrypt_pack(const PackedCiphertext &ciphertext) {
    Vec<Integer> result;
    const size_t n_plaintexts = ciphertext.n_plaintexts;

    if(n_plaintexts > get_packCount(ciphertext.plaintext_bits))
        error_exit("trying to unpack too many elements!");

    result.SetLength(n_plaintexts);
    this->decrypt_pack(ciphertext, result.begin(), result.end());
    return result;
}


/*------------------ New Private key -----------------*/

NewPrivateKey::NewPrivateKey(const size_t key_size_bits_, const Integer &p_, const Integer &q_)
        : key_size_bits(key_size_bits_),
            p(p_),
            q(q_) { }

NewPrivateKey::NewPrivateKey(const size_t key_size_bits_, const size_t a_bits_, const Integer &p_, const Integer &q_, const Integer &a_)
        : key_size_bits(key_size_bits_),
            a_bits(a_bits_),
            p(p_),
            q(q_),
            a(a_){
    // mpz_invert(inv_a.get_mpz_t(), a.get_mpz_t(),(p*q).get_mpz_t());
    inv_a = a.inv_mod_n(p*q);// a<n -> mod n

    Integer lambda = Integer::lcm(p - 1, q - 1);
    auto n = q*p;   
    auto n2 = n*n;
    mpz_powm(nkg.get_mpz_t(),Integer(2).get_mpz_t(),lambda.get_mpz_t(),n2.get_mpz_t());
    nkg = (nkg - 1) % n2;
    inv_kg = (nkg/n).inv_mod_n(n);
    // inv_nkg = nkg.inv_mod_n(n);// 代确认: 是否直接 mod n 即可
}

NewPrivateKey::NewPrivateKey() : key_size_bits(0) { }

bool NewPrivateKey::operator==(const NewPrivateKey &input) const {
    return p == input.p && q == input.q && a == input.a;
}

NewPrivateKey &NewPrivateKey::operator=(const NewPrivateKey &priv_) {
    this->key_size_bits = priv_.key_size_bits;
    this->a_bits = priv_.a_bits;
    this->p = priv_.p;
    this->q = priv_.q;
    this->a = priv_.a;
    this->inv_a = priv_.inv_a;
    this->nkg = priv_.nkg;
    this->inv_kg = priv_.inv_kg;
    return *this;
}

const std::string NewPrivateKey::to_string(const bool brief) const {
    std::ostringstream o("");
    o << "<NewPrivateKey[";
    o << key_size_bits << "]";
    o << " p=" << p.to_string(brief);
    o << " q=" << q.to_string(brief);
    o << " a=" << a.to_string(brief);
    o << " inv_a=" << inv_a.to_string(brief);
    o << " nkg=" << nkg.to_string(brief);
    o << ">";
    return o.str();
}


/*------------ random --------------*/

const Integer NewPaillier::Randomizer::get_noise() const {
    if(!precomputed)
        error_exit("lookup table not precomputed!");

    if(paillier->fast_mod) {
        return paillier->fast_mod.get()->pow_mod_n2(g_pow_n, r());
    } else {
        return g_pow_n.pow_mod_n(r(), paillier->n2);
    }
}

NewPaillier::Randomizer::Randomizer(const NewPaillier *paillier_)
        : paillier(paillier_) { }

void NewPaillier::Randomizer::precompute() {
    if(paillier->fast_mod) {
        g_pow_n = paillier->fast_mod.get()->pow_mod_n2(paillier->pub.g, paillier->pub.n);
    } else {
        g_pow_n = paillier->pub.g.pow_mod_n(paillier->pub.n, paillier->n2);
    }
    precomputed = true;
}

const Integer NewPaillier::Randomizer::r() const {
    return Random::instance().rand_int_bits(paillier->r_bits);
}

const std::string NewPaillier::Randomizer::to_string(const bool brief) const {
    std::ostringstream o("");
    o << "<Randomizer";
    o << " g_pow_n=" << g_pow_n.to_string(brief);
    o << " precomputed=" << precomputed;
    o << ">";

    return o.str();
}

NewPaillier::FastRandomizer::FastRandomizer(const NewPaillier *paillier_, const size_t r_lut_size_, const size_t r_use_count_)
        : Randomizer(paillier_),
            r_lut_size(r_lut_size_),
            r_use_count(r_use_count_) { }

void NewPaillier::FastRandomizer::precompute() {
    Randomizer::precompute();

    #ifdef DEBUG
    std::cerr << "PaillierFast::FastRandomizer::precompute: precomputing " << r_lut_size<< std::endl;
    #endif

    gn_pow_r.reserve(r_lut_size);

    omp_declare_lock(writelock);
    omp_init_lock(&writelock);
    if(paillier->fast_mod) {
        paillier->fast_mod.get()->pow_mod_n2(g_pow_n, r());

        #pragma omp parallel for
        for(auto i = 0u; i < r_lut_size; i++) {
            const auto rand = paillier->fast_mod.get()->pow_mod_n2(g_pow_n, r());

            omp_set_lock(&writelock);
            gn_pow_r.push_back(rand);
            omp_unset_lock(&writelock);
        }
    } else {
        #pragma omp parallel for
        for(auto i = 0u; i < r_lut_size; i++) {
            const auto rand = g_pow_n.pow_mod_n(r(), paillier->n2);

            omp_set_lock(&writelock);
            gn_pow_r.push_back(rand);
            omp_unset_lock(&writelock);
        }
    }
    omp_destroy_lock(&writelock);

    #ifdef DEBUG
    size_t size = 0;
    for(auto i = 0u; i < r_lut_size; i++) {
        size += gn_pow_r[i].size_bits();
    }
    std::cerr << "PaillierFast::FastRandomizer::precompute: computed r, size=" << size / 8 << " bytes" << std::endl;
    #endif
}

const Integer NewPaillier::FastRandomizer::get_noise() const {
    if(!precomputed)
        error_exit("lookup table not precomputed!");
    Integer ret = 1;
    Random &rand = Random::instance();
    for(auto i = 0u; i < r_use_count; i++) {
        unsigned int ix = rand.rand_int(r_lut_size).to_uint();
        ret = (ret * gn_pow_r[ix]) % paillier->n2;
    }

    return ret;
}

const std::string NewPaillier::FastRandomizer::to_string(const bool brief) const {
    std::ostringstream o("");
    o << "<FastRandomizer";
    o << " g_pow_n=" << g_pow_n.to_string(brief);
    o << " r_lut_size=" << r_lut_size;
    o << " r_use_count=" << r_use_count;
    o << " precomputed=" << precomputed;
    o << ">";

    return o.str();
}