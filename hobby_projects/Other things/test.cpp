    bool is_factor(const spinor_sum & sum, const spinor_sum &factor, spinor_sum leftover_sum){
        spinor var = factor.terms[0].numerator;
        Polynomial poly(sum,var);
        spinor_sum remainder = sum;
        
        int i = 0;
        while(poly.highest_order > i){
            remainder = remainder - (poly[poly.highest_order - i]/var) * factor;
            leftover_sum = leftover_sum + (poly[poly.highest_order - i]/var);

            if(remainder.is_null == true){
                return true;
            }
            i++;
        }
        return false;
    }

    factorized_spinor factor(factorized_spinor &sum, std::vector<factor_struct> possible_factors){
        if (possible_factors.empty() or sum.remainder.terms.size() <= 1) {
            return sum;
        }
        factor_struct current_factor = possible_factors[0];

        int N = 1 << possible_factors.size();
        for (int mask = 1; mask < N; ++mask) {
            spinor_sum candidate_factor(current_factor.factor);
            for (int i = 1; i < possible_factors.size(); ++i) {
                if(mask & (1 << i)){
                   candidate_factor = candidate_factor + possible_factors[i].factor;
                }
            }
            spinor_sum remainder = sum.remainder;
            spinor_sum leftover_factors;

            if (is_factor(remainder,candidate_factor,leftover_factors)){
                sum.factor.push_back(candidate_factor);
                sum.remainder = leftover_factors;
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




    /// ----------------------
    spinor_sum factor_out(const spinor_sum & sum, const spinor &s){
        spinor_sum result;
        for (spinor_product p : sum.terms){
            int idx;
            spinor_product temp;
            temp.denominator = p.denominator;
            temp.factor = p.factor;
            auto it = std::find(p.numerator.begin(),p.numerator.end(),s);
            if (it != p.numerator.end()){
                idx = std::distance(p.numerator.begin(),it);
                for (int i = 0; i < p.numerator.size(); ++i){
                    if(i != idx){
                        temp = temp * p.numerator[i];
                    }
                }
                result = result + temp;
            }
        }

        return result;
    }

    factorized_spinor factor(factorized_spinor &sum, std::vector<factor_struct> possible_factors){
        if (possible_factors.empty() or sum.remainder.terms.size() <= 1) {
            return sum;
        }
        factor_struct current_factor = possible_factors[0];

        int N = 1 << possible_factors.size();
        for (int mask = 1; mask < N; ++mask) {
            spinor_sum candidate_factor(current_factor.factor);
            for (int i = 1; i < possible_factors.size(); ++i) {
                if(mask & (1 << i)){
                   candidate_factor = candidate_factor + possible_factors[i].factor;
                }
            }
            spinor_sum remainder = sum.remainder;
            spinor_sum leftover_factors;
            for (spinor_product p : candidate_factor.terms){
                spinor temp = p.numerator[0];
                spinor_sum factored_out = factor_out(sum.remainder,temp);
                remainder = remainder - temp * factored_out;
                leftover_factors = factored_out;
            }
            if(remainder.is_null){
                    sum.factor.push_back(candidate_factor);
                    sum.remainder = leftover_factors;
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