 //   struct TermWithCount {
 //       spinor_product term;
//        int count;

//        bool operator<(const TermWithCount &other) const {
//            return count < other.count; // or > for descending
//        }
//    };

//struct polynomial {
//        std::vector<TermWithCount> terms;

//        polynomial(const spinor_sum &sum, const spinor &var) {
//           for (const auto& p : sum.terms) {
//                int count = 0;
//                for (const auto& s : p.numerator) {
//                    if (s == var) ++count;
//                }
//                terms.push_back({p, count});
//            }
//        std::sort(terms.begin(), terms.end());
//        }
//    };

//    bool is_long_divisable(spinor_sum &sum, spinor_sum &factor, const spinor &var){
//        polynomial ordered_polynomial(sum,var);
//        int highest_order = ordered_polynomial.terms[0].count;
//        for(int i = 0; i < highest_order; ++i){

//        }

//    }

    factorized_spinor factor(factorized_spinor &sum, std::vector<factor_struct> possible_factors){
        if (possible_factors.empty()) {
            return sum;
        }
        factor_struct current_factor = possible_factors[0];

        if (current_factor.factor_count == sum.original_sum_length){
            spinor_sum true_factor(current_factor.factor);
            sum.factor.push_back(true_factor);
            for(spinor_product &p : sum.remainder.terms){
                p = p / current_factor.factor;
                reduce_product(p);
                if (std::find(p.numerator.begin(),p.numerator.end(),current_factor.factor) == p.numerator.end()){
                    current_factor.factor_count -= 1;
                }
                if (current_factor.factor_count == 0){
                    possible_factors.erase(possible_factors.begin());
                }
                else{
                    std::sort(possible_factors.begin(),possible_factors.end());
                }
            }
            return factor(sum,possible_factors);
        }


        int N = 1 << (possible_factors.size() - 1);
        for (int mask = 1; mask < N - 1; ++mask) {
            spinor_sum new_factor(possible_factors[0].factor);
            spinor_sum remainder;
            for (int i = 1; i < possible_factors.size(); ++i) {
                if(mask & (1 << (i-1))){
                   new_factor.terms.push_back(possible_factors[i].factor);
                }
                else{
                    remainder.terms.push_back(possible_factors[i].factor);
                }
            }
            std::cout << "Factor:" << new_factor << " , " << "Remainder" << remainder << std::endl;
            if(new_factor * remainder == sum.remainder){
                    sum.factor.push_back(new_factor);
                    sum.remainder = remainder;
                    return factor(sum,possible_factors);
            }
        }
        if(sum.order == 0){
            possible_factors.erase(possible_factors.begin());
            sum.order = sum.original_sum_length;
        }
        else{
            std::rotate(possible_factors.begin(), possible_factors.end() - 1, possible_factors.end());
            sum.order -= 1;
        }
        
        return factor(sum,possible_factors);
    }


    spinor_sum factor_out(const spinor_sum & sum, const spinor s){
        spinor_sum result;
        for (spinor_product p : sum){
            int idx;
            spinor_product temp;
            temp.denominator = p.denominator;
            temp.factor = p.factor;
            auto iter std::find(p.numerator.begin(),p.numerator.end(),s);
            if (iter != p.numerator.end()){
                idx = std::distance(p.begin(),iter);
                for (int i = 0; p.numerator.size(); ++i){
                    if(i != idx){
                        temp *= p.numerator[i];
                    }
                }
                result += temp;
            }
        }

        return result;
    }