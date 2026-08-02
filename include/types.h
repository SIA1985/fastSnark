#ifndef TYPES_H
#define TYPES_H

#include "nlohmann/json.hpp"
#include <openssl/sha.h>
#include <openssl/evp.h>

#define MCL_USE_GMP 0
#define MCL_USE_VINT 1
#include <mcl/bn.hpp>

namespace snrk {
    typedef nlohmann::json json_t;
    class CircutParams;

    typedef int GateType;
}

namespace in {

typedef std::size_t address_t;

typedef unsigned long witness_t;

class witnesses_t: public std::vector<witness_t>
{
private:
    witnesses_t() = default;
    using std::vector<witness_t>::operator=;

    friend witnesses_t genWitnesses(witness_t start, std::size_t count, witness_t wStep);
    friend class snrk::CircutParams;

    friend void to_json(snrk::json_t& j, const witnesses_t& ws);
    friend void from_json(const snrk::json_t& j, witnesses_t& ws);
};

typedef mcl::Fr value_t;
unsigned long getUL(value_t v);
void negative(value_t &v);

typedef std::vector<value_t> values_t;

typedef value_t X_t;
typedef value_t Y_t;

template<typename X = X_t, typename Y = Y_t>
struct DotType{X x; Y y;};

typedef DotType<X_t, Y_t> dot_t;
bool operator<(const dot_t &a, const dot_t &b);
typedef std::vector<dot_t> dots_t;

typedef mcl::bn::G1 G1;
typedef mcl::bn::G2 G2;
typedef mcl::bn::GT GT;
typedef G1 commit_t;

typedef std::vector<commit_t> commits_t;
typedef std::vector<G1> keys_t;

typedef std::vector<X_t> xs_t;

class InterpolationPolynom;
typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;
typedef InterpolationPolynom W_t;

typedef std::array<unsigned char, SHA256_DIGEST_LENGTH> hash_t;
std::string operator+(const hash_t &a, const hash_t &b);
std::string operator+(const std::string &s, const hash_t &h);
std::string operator+(const hash_t &h, const std::string &s);
typedef std::vector<hash_t> hashes_t;

class hasher {
public:
    hasher();
    ~hasher();

    hash_t operator()(std::string toHash);
    hash_t operator()(hash_t a, hash_t b);

private:
    EVP_MD_CTX *m_context;
};

hash_t hash(std::string toHash);
value_t cut(value_t r, value_t max);

}

namespace snrk {
    class Jsonable
    {
    public:
        virtual snrk::json_t toJson() const;
        virtual bool fromJson(const snrk::json_t &json);
    };

    class MerkleTree {
    public:
        class Proof_t : public Jsonable
        {
        public:
            virtual snrk::json_t toJson() const override;
            virtual bool fromJson(const snrk::json_t &json) override;

            in::hashes_t path;
            std::size_t index;
        } ;

        class MultiProof_t : public Jsonable
        {
        public:
            virtual snrk::json_t toJson() const override;
            virtual bool fromJson(const snrk::json_t &json) override;

            //todo: оптимизация
            std::vector<Proof_t> proofs;
            std::vector<std::size_t> indicies;
        };


    public:
        MerkleTree(const std::vector<std::string> &data);

        in::hash_t root() const;

        std::optional<Proof_t> proof(std::size_t index) const;

        static bool verify(const Proof_t &proof, in::hash_t leaf, in::hash_t root);

        std::optional<MultiProof_t> multiProof(const std::vector<std::size_t> &indices) const;

        static bool multiVerify(const MultiProof_t &proof, in::hashes_t leafs, in::hash_t root);

    private:
        std::vector<in::hashes_t> m_tree;
    };

    struct ProverParams {in::keys_t keys;};
    struct VerifierParams {in::G1 g1; in::G2 g2, tG2;};

    template<typename T>
    std::optional<T> fromString(const std::string &string);

    std::string toString(const in::value_t &value);
    template<>
    std::optional<in::value_t> fromString<in::value_t>(const std::string &string);

    std::string toString(const in::commit_t &commit);
    template<>
    std::optional<in::commit_t> fromString<in::commit_t>(const std::string &string);

    std::string toString(const in::hash_t &hash);
    template<>
    std::optional<in::hash_t> fromString<in::hash_t>(const std::string &string);
}

#endif // TYPES_H
