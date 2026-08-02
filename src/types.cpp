#include "types.h"

#include <set>

namespace in
{

bool operator<(const dot_t &a, const dot_t &b)
{
    return a.x <= b.x; // != 1;
}

unsigned long getUL(value_t v)
{
    uint8_t buf[64] = {0};

    size_t written = v.getLittleEndian(buf, sizeof(buf));

    if (written == 0) return 0;

    unsigned long result = 0;
    std::memcpy(&result, buf, sizeof(result));

    return result;
}

void negative(value_t &v)
{
    mcl::Fr::neg(v, v);
}

hasher::hasher()
{
    m_context = EVP_MD_CTX_new();
}

hasher::~hasher()
{
    EVP_MD_CTX_free(m_context);
}

hash_t hasher::operator()(std::string toHash)
{
    hash_t result;
    unsigned int length = 0;

    EVP_DigestInit_ex(m_context, EVP_sha256(), nullptr);
    EVP_DigestUpdate(m_context, toHash.c_str(), toHash.size());
    EVP_DigestFinal_ex(m_context, result.data(), &length);

    return result;
}

hash_t hasher::operator()(hash_t a, hash_t b)
{
    hash_t result;
    unsigned int length = 0;

    EVP_DigestInit_ex(m_context, EVP_sha256(), nullptr);

    EVP_DigestUpdate(m_context, a.data(), a.size());
    EVP_DigestUpdate(m_context, b.data(), b.size());

    EVP_DigestFinal_ex(m_context, result.data(), &length);

    return result;
}

hash_t hash(std::string toHash)
{
    hasher h;

    return h(toHash);
}

std::string operator+(const in::hash_t &a, const in::hash_t &b)
{
    return std::string{a.begin(), a.end()}.append({b.begin(), b.end()});
}

std::string operator+(const std::string &s, const in::hash_t &h)
{
    return std::string{s}.append({h.begin(), h.end()});
}

std::string operator+(const in::hash_t &h, const std::string &s)
{
    return std::string{h.begin(), h.end()}.append(s);
}


value_t cut(value_t r, value_t max)
{
    return getUL(r) % (getUL(max) + 1);
}

}


namespace snrk {

snrk::json_t Jsonable::toJson() const
{
    return "{}";
}

bool Jsonable::fromJson(const snrk::json_t &json)
{
    return false;
}

const std::string Path = "Path";
const std::string Index = "Index";

snrk::json_t MerkleTree::Proof_t::toJson() const
{
    snrk::json_t json;

    snrk::json_t::array_t arr;
    for(const auto &h : path) {
        arr.push_back(toString(h));
    }
    json[Path] = arr;

    json[Index] = index;

    return json;
}

bool MerkleTree::Proof_t::fromJson(const snrk::json_t &json)
{
    if ((!json.contains(Path) || !json[Path].is_array()) &&
        (!json.contains(Index) || !json[Index].is_number_integer())) {
        return false;
    }

    json_t::array_t items = json[Path];
    path.clear();

    path.reserve(items.size());
    for(const auto &item : items) {
        if (!item.is_string()) {
            return false;
        }

        auto opt = fromString<in::hash_t>(item.get<std::string>());
        if (!opt) {
            return false;
        }

        path.push_back(opt.value());
    }

    index = json[Index];

    return true;
}

const std::string Proofs = "Proofs";
const std::string Indicies = "Indicies";

snrk::json_t MerkleTree::MultiProof_t::toJson() const
{
    snrk::json_t json;

    snrk::json_t::array_t arr;
    for(const auto &proof : proofs) {
        arr.push_back(proof.toJson());
    }
    json[Proofs] = arr;

    json[Indicies] = indicies;

    return json;
}

bool MerkleTree::MultiProof_t::fromJson(const snrk::json_t &json)
{
    if ((!json.contains(Proofs) || !json[Proofs].is_array()) &&
        (!json.contains(Indicies) || !json[Indicies].is_array())) {
        return false;
    }

    json_t::array_t items = json[Proofs];

    proofs.clear();
    proofs.reserve(items.size());
    for(const auto &item : items) {
        Proof_t proof;
        if (!proof.fromJson(item)) {
            return false;
        }

        proofs.push_back(std::move(proof));
    }

    for(const auto& idx : json[Indicies]) {
        if (!idx.is_number_integer()) {
            return false;
        }
    }

    indicies = json[Indicies].get<std::vector<std::size_t>>();

    return true;
}

MerkleTree::MerkleTree(const std::vector<std::string> &data)
{
    if (data.empty()) {
        return;
    }


    in::hasher h;
    m_tree.reserve(std::ceil(std::log2(data.size())) + 1);

    in::hashes_t lay;
    for(const auto &d : data) {
        lay.push_back(in::hash(d));
    }

    m_tree.push_back(lay);

    while(lay.size() > 1) {
        if (lay.size() % 2 == 1) {
            lay.push_back(lay.back());
        }

        in::hashes_t temp;
        for(std::size_t i = 0; i < lay.size(); i += 2) {
            temp.push_back(h(lay[i], lay[i + 1]));
        }

        lay = temp;
        m_tree.push_back(lay);
    }

}

in::hash_t MerkleTree::root() const
{
    return m_tree.back().back();
}

std::optional<MerkleTree::Proof_t> MerkleTree::proof(std::size_t index) const
{
    if (m_tree.empty()) {
        return std::nullopt;
    }

    Proof_t proof;

    proof.index = index;

    for (std::size_t i = 0; i < m_tree.size() - 1; ++i) {
        const auto &layer = m_tree[i];

        std::size_t pairIndex = (index % 2 == 0) ? index + 1 : index - 1;

        if (pairIndex >= layer.size()) {
            pairIndex = index;
        }

        proof.path.push_back(layer[pairIndex]);
        index /= 2;
    }

    return proof;
}

bool MerkleTree::verify(const Proof_t &proof, in::hash_t leaf, in::hash_t root)
{
    in::hash_t result = leaf;
    auto index = proof.index;

    for(const auto &h : proof.path) {
        if (index % 2 == 0) {
            result = in::hash(in::operator+(result, h));
        } else {
            result = in::hash(in::operator+(h, result));
        }

        index /= 2;
    }

    return root == result;
}

std::optional<MerkleTree::MultiProof_t> MerkleTree::multiProof(const std::vector<std::size_t> &indices) const
{
    if (m_tree.empty() || indices.empty()) {
        return std::nullopt;
    }

    std::vector<std::size_t> unique;
    std::set<std::size_t> seen;

    for(auto i : indices) {
        if (seen.count(i) != 0) {
            continue;
        }

        unique.push_back(i);
        seen.insert(i);
    }

    MultiProof_t mproof;
    mproof.indicies = indices;

    for(auto indx : unique) {
        auto opt = proof(indx);
        if (!opt) {
            return std::nullopt;
        }

        mproof.proofs.push_back(opt.value());
    }

    return mproof;

}

bool MerkleTree::multiVerify(const MultiProof_t &mproof, in::hashes_t leafs, in::hash_t root)
{
    if (mproof.indicies.size() != leafs.size() || mproof.indicies.empty()) {
        return false;
    }

    std::set<std::size_t> visited;
    in::hashes_t unique;
    for(std::size_t i = 0; i < mproof.indicies.size(); ++i) {
        auto indx = mproof.indicies[i];

        if (visited.count(indx) != 0) {
            continue;
        }

        unique.push_back(leafs[i]);
        visited.insert(indx);
    }

    bool result = true;
    for(std::size_t i = 0; i < unique.size(); i++) {
        result = result && MerkleTree::verify(mproof.proofs[i], unique[i], root);
    }

    return result;
}


std::string toString(const in::value_t &value)
{
    std::string result;

    value.getStr(result);

    return result;
}

template<>
std::optional<in::value_t> fromString<in::value_t>(const std::string &string)
{
    if (string.empty()) {
        return std::nullopt;
    }

    in::value_t value;

    value.setStr(string);

    return {value};
}

constexpr size_t G1_BINARY_SIZE = 64;
std::string toString(const in::commit_t &commit)
{
    return commit.getStr();
}

template<>
std::optional<in::commit_t> fromString<in::commit_t>(const std::string &string)
{
    if (string.empty()) {
        return std::nullopt;
    }

    in::commit_t result;

    result.setStr(string);

    return {result};
}

std::string toString(const in::hash_t &hash)
{
    std::stringstream ss;
    ss << std::hex << std::setfill('0');
    for (unsigned char byte : hash) {
        ss << std::setw(2) << static_cast<int>(byte);
    }
    return ss.str();
}

template<>
std::optional<in::hash_t> fromString<in::hash_t>(const std::string &string)
{
    if (string.length() != SHA256_DIGEST_LENGTH * 2) {
        return std::nullopt;
    }

    in::hash_t result;

    for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        std::string byteString = string.substr(i * 2, 2);

        unsigned int byteValue;
        std::stringstream ss;
        ss << std::hex << byteString;
        ss >> byteValue;

        result[i] = static_cast<unsigned char>(byteValue);
    }

    return result;
}

}
