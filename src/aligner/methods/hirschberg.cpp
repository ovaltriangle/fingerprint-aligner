/**
 * hirschberg
 * @author astaroth
 * @date 10/10/20
 *
 * @brief
 */

#include "hirschberg.hpp"

namespace hirschberg {
    std::vector<long> nw_score(const std::string &a, const std::string &b, const WeightTable &weights) {
        auto alen {a.size() + 1}, blen {b.size() + 1};

        /// do something similar: two lines only. at the end, return the last one. this makes it possible to have a
        /// 2*max(n, m) = max(n, m) linear space
        std::vector<long> score_h (alen), score_l (alen);

        for (std::size_t i {1}; i < alen; ++i) {  // decreasing first
            score_h.at(i) = score_h.at(i - 1) + weights.insert_cost;
        }

        for (std::size_t i {1}; i < blen; ++i) {
            // decreasing column
            score_l.at(0) = score_h.at(0) + weights.delete_cost;
            for (std::size_t j {1}; j < alen; ++j) {
                // take the confronting nts
                auto asel {a.at(j - 1)}, bsel {b.at(i - 1)};

                // diag
                auto sub {score_h.at(j - 1) + ((asel == bsel) ? weights.match_cost : weights.substitution_cost)};
                // left
                auto ins {score_l.at(j - 1) + weights.insert_cost};
                // up
                auto del {score_h.at(j) + weights.delete_cost};

                score_l.at(j) = std::max(sub, std::max(ins, del));
            }
            score_h.swap(score_l);
        }

        return score_h;
    }

    std::pair<std::string, std::string> align(const std::string &a, const std::string &b, const WeightTable &weights) {
        std::string z {}, w {};

        if (a.empty()) {
            for (auto const &ch: b) {
                z += '-';
                w += ch;
            }
        } else if (b.empty()) {
            for (auto const &ch: a) {
                z += ch;
                w += '-';
            }
        } else if (a.size() == 1 or b.size() == 1) {
            // NW special case emulation
            if (a.size() == 1 and b.size() == 1) {
                return {a, b};
            } else if (a.size() == 1) {
                for (auto const &ch: b) {
                    z += '-';
                    w += ch;
                }
                std::size_t search {b.find(a)};
                std::size_t index {(search == std::string::npos) ? 0 : search};
                z.at(index) = a.at(0);
            } else if (b.size() == 1) {
                for (auto const &ch: a) {
                    z += ch;
                    w += '-';
                }
                std::size_t search {a.find(b)};
                std::size_t index {(search == std::string::npos) ? 0 : search};
                w.at(index) = b.at(0);
            }
        } else {
            std::size_t amid {a.size() / 2};

            std::string rev_b {b};
            std::reverse(rev_b.begin(), rev_b.end());

            std::string aL {a.substr(0, amid)}, aR {a.substr(amid)};

            std::reverse(aR.begin(), aR.end());

            /// CHECK: are these launched together?
            auto scoreL {nw_score(b, aL, weights)}, scoreR {nw_score(rev_b, aR, weights)};

            std::reverse(scoreR.begin(), scoreR.end());

            std::size_t max {(scoreL.size() > scoreR.size()) ? scoreL.size() : scoreR.size()};
            std::vector<long> vectorial_sum (max);

            for (size_t i {0}; i < scoreL.size(); ++i) {
                vectorial_sum.at(i) += scoreL.at(i);
            }
            for (size_t i {0}; i < scoreR.size(); ++i) {
                vectorial_sum.at(i) += scoreR.at(i);
            }

            /// Save up some space, man!
            scoreL.clear(); scoreR.clear();

            auto bmid {utils::argmax(vectorial_sum)};
            //auto bmid {std::distance(vectorial_sum.begin(), std::max_element(vectorial_sum.begin(), vectorial_sum.end()))};

            /// Even more
            vectorial_sum.clear();

            std::string bL {b.substr(0, bmid)}, bR {b.substr(bmid)};

            std::reverse(aR.begin(), aR.end());

            auto hirsL {align(aL, bL, weights)}, hirsR {align(aR, bR, weights)};

            z = hirsL.first + hirsR.first;
            w = hirsL.second + hirsR.second;
        }

        return {z, w};
    }

    double score(const std::string &a, const std::string &b, const WeightTable &weights) {
        /// Perform Hirschberg alignment
        auto hirsch{align(a, b, weights)};

        /// Calculate the score of the alignment
        std::size_t score {0};
        for (std::size_t i {0}; i < hirsch.first.size(); ++i) {
            if (hirsch.first.at(i) == '-' and hirsch.second.at(i) == '-')
                score += 1;
            else if (hirsch.first.at(i) == hirsch.second.at(i))
                score += 2;
        }

        return static_cast<double>(score) / (static_cast<double>(a.size()) * 2);
    }

    AlignerResults score_full(const std::string &a, const std::string &b, const WeightTable &weights) {
        /// Prepare the returned value
        AlignerResults ar;

        /// Perform Hirschberg alignment
        auto hirsch{align(a, b, weights)};
        ar.aligned_a = hirsch.first;
        ar.aligned_b = hirsch.second;

        /// Calculate the result of the alignment
        for (std::size_t i{0}; i < hirsch.first.size(); ++i) {
            if (hirsch.first.at(i) == '-' and hirsch.second.at(i) == '-')
                ar.score += 1;
            else if (hirsch.first.at(i) == hirsch.second.at(i)) {
                ar.score += 2;
            }
        }

        /// Calculate the ratio
        ar.score /= static_cast<double>(a.size() * 2);

        return ar;
    }
}  // namespace aligners