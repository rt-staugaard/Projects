#include "spinor_header.h"
#include <numeric>

namespace Spinor{

    struct Polynomial {
        std::map<int, spinor_sum> terms;
        int highest_order = 0;
        spinor var;
        
        Polynomial(const spinor_sum &sum, const spinor &var) {
            this->var = var;
            for (spinor_product prod : sum.terms) {
                int order = 0;
                for(spinor s : prod.numerator){
                    if (is_similar_spinor(s,var)){
                        order += 1;
                    }
                }

                terms[order].terms.push_back(prod);

                if (order > highest_order){
                    highest_order = order;
                }
            }
        }

        Polynomial operator-(const spinor_sum &sum){
            Polynomial result = *this;
            for (spinor_product prod : sum.terms) {
                int order = 0;
                for(spinor s : prod.numerator){
                    if (is_similar_spinor(s,var)){
                        order += 1;
                    }
                }

                result.terms[order] = result.terms[order] - prod;
            }
            return result;
        }

        bool is_null(){
            for(int i = 0; i < highest_order; i++){
                if(!terms[i].is_null){
                    return false;
                }
            }
            return true;
        }
    };

    // ---------- Addition ----------
    spinor_sum operator+(const spinor &s1, const spinor &s2){
        spinor_sum result(s1,s2);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor_product &p1, const spinor &s2){
        spinor_product p2(s2);
        spinor_sum result(p1,p2);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor &s1,const spinor_product &p2){
        spinor_product p1(s1);
        spinor_sum result(p1,p2);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor_product &p1, const spinor_product &p2){
        spinor_sum result(p1,p2);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor_sum &sum1, const spinor_product &p2){
        spinor_sum result = sum1;
        result.terms.push_back(p2);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor_product &p1, const spinor_sum &sum2){
        spinor_sum result = sum2;
        result.terms.push_back(p1);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor_sum &sum1, const spinor &s2){
        spinor_sum result = sum1;
        spinor_product temp(s2);
        result.terms.push_back(temp);
        cancellation(result);
        return result;
    }
    spinor_sum operator+(const spinor &s1,const spinor_sum &sum2){
        spinor_sum result = sum2;
        spinor_product temp(s1);
        result.terms.push_back(temp);
        cancellation(result);
        return result;
    }
    
    spinor_sum operator+(const spinor_sum &sum1, const spinor_sum &sum2){
        spinor_sum result = sum1;
        for (int i = 0; i < sum2.terms.size(); ++i){
            result.terms.push_back(sum2.terms[i]);
        }
        cancellation(result);
        return result;
    }

    spinor_fraction operator+(const spinor_fraction &f1, const spinor &s2){
        spinor_fraction f2(s2);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor &s1, const spinor_fraction &f2){
        spinor_fraction f1(s1);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_product &p2){
        spinor_fraction f2(p2);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor_product &p1, const spinor_fraction &f2){
        spinor_fraction f1(p1);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_sum &sum2){
        spinor_fraction f2(sum2);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor_sum &sum1, const spinor_fraction &f2){
        spinor_fraction f1(sum1);
        spinor_fraction result = f1 + f2;
        return result;
    }
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_fraction &f2){
        spinor_fraction result;
        result.denominator = f1.denominator * f2.denominator;
        result.numerator = (f1.numerator*f2.denominator) + (f2.numerator*f1.denominator);
        return result;
    }

// ---------- Subtraction ----------
    spinor_sum operator-(const spinor &s1, const spinor &s2){
        spinor_sum result(s1,s2);
        result.terms[1].factor *= -1;
        cancellation(result);
        return result;
    }
    spinor_sum operator-(const spinor_product &p1, const spinor &s2){
        spinor_product p2(s2);
        spinor_sum result(p1,p2);
        result.terms[1].factor *= -1;
        cancellation(result);
        return result;
    }
    spinor_sum operator-(const spinor &s1,const spinor_product &p2){
        spinor_product p1(s1);
        spinor_sum result(p1,p2);
        result.terms[1].factor *= -1;
        cancellation(result);
        return result;
    }
    spinor_sum operator-(const spinor_product &p1, const spinor_product &p2){
        spinor_sum result(p1,p2);
        result.terms[1].factor *= -1;
        cancellation(result);
        return result;
    }
    spinor_sum operator-(const spinor_sum &sum1, const spinor_product &p2){
        spinor_sum sum2(p2);
        spinor_sum result = sum1 - sum2;
        return result;
    }
    spinor_sum operator-(const spinor_product &p1, const spinor_sum &sum2){
        spinor_sum sum1(p1);
        spinor_sum result = sum1 - sum2;
        return result;
    }
    spinor_sum operator-(const spinor_sum &sum1, const spinor &s2){
        spinor_sum sum2(s2);
        spinor_sum result = sum1 - sum2;
        return result;
    }
    spinor_sum operator-(const spinor &s1,const spinor_sum &sum2){
        spinor_sum sum1(s1);
        spinor_sum result = sum1 - sum2;
        return result;
    }

    spinor_sum operator-(const spinor_sum &sum1, const spinor_sum &sum2){
        spinor_sum result = sum1;
        spinor_product temp;
        for (int i = 0; i < sum2.terms.size(); ++i){
            temp = sum2.terms[i];
            temp.factor *= -1;
            result.terms.push_back(temp);
        }
        cancellation(result);
        return result;
    }
    
    spinor_fraction operator-(const spinor_fraction &f1, const spinor &s2){
        spinor_fraction f2(s2);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor &s1, const spinor_fraction &f2){
        spinor_fraction f1(s1);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_product &p2){
        spinor_fraction f2(p2);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor_product &p1, const spinor_fraction &f2){
        spinor_fraction f1(p1);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_sum &sum2){
        spinor_fraction f2(sum2);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor_sum &sum1, const spinor_fraction &f2){
        spinor_fraction f1(sum1);
        spinor_fraction result = f1 - f2;
        return result;
    }
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_fraction &f2){
        spinor_fraction result;
        result.denominator = f1.denominator * f2.denominator;
        result.numerator = (f1.numerator*f2.denominator) - (f2.numerator*f1.denominator);
        return result;
    }

// ---------- Multiplication ----------

    spinor_product operator*(const spinor &s1, const spinor &s2){
        spinor_product result{s1,s2};
        result.factor = s1.factor * s2.factor;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator*(const spinor_product &p1, const spinor &s2){
        spinor_product result = p1;
        result.numerator.push_back(s2);
        result.factor *= s2.factor;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator*(const spinor &s1, const spinor_product &p2){
        spinor_product result = p2;
        result.numerator.push_back(s1);
        result.factor *= s1.factor;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator*(const spinor_product &p1, const spinor_product &p2){
        spinor_product result = p1;
        result.numerator.insert(result.numerator.end(), p2.numerator.begin(), p2.numerator.end());
        result.factor *= p2.factor;
        order(result);
        reduce_product(result);
        return result;
    }

    spinor operator*(const double &c, const spinor &s){
        spinor result = s;
        result.factor *= c;
        return result;
    }

    spinor operator*(const spinor &s, const double &c){
        spinor result = s;
        result.factor *= c;
        return result;
    }

    spinor_product operator*(const double &c, const spinor_product &p){
        spinor_product result = p;
        result.factor *= c;
        return result;
    }

    spinor_product operator*(const spinor_product &p, const double &c){
        spinor_product result = p;
        result.factor *= c;
        return result;
    }

    spinor_sum operator*(const spinor_product &p1,const spinor_sum &sum2){
        spinor_sum sum1(p1);
        spinor_sum result = sum1 * sum2; 
        return result;
    }
    spinor_sum operator*(const spinor_sum &sum1,const spinor_product &p2){
        spinor_sum sum2(p2);
        spinor_sum result = sum1 * sum2; 
        return result;
    }

    spinor_sum operator*(const spinor_sum &sum1,const spinor_sum &sum2){
        spinor_sum result;
        spinor_product temp;
        if (sum1.terms.empty() or sum2.terms.empty()){
            if (sum1.is_null or sum2.is_null){
                result.is_null = true;
            }
            else if(sum1.terms.empty()){
                result = sum2;
            }
            else{
                result = sum1;
            }
            return result;
        }

        for (int i = 0; i < sum1.terms.size(); ++i){
            for (int j = 0; j < sum2.terms.size(); ++j){
                temp = sum1.terms[i]*sum2.terms[j];
                reduce_product(temp);
                result.terms.push_back(temp);
            }
        }
        cancellation(result);
        return result;
    }

    spinor_sum operator*(const spinor &s1,const spinor_sum &sum2){
        spinor_sum sum1(s1);
        spinor_sum result = sum1 * sum2; 
        return result;
    }

    spinor_sum operator*(const spinor_sum &sum1,const spinor &s2){
        spinor_sum sum2(s2);
        spinor_sum result = sum1 * sum2; 
        return result;
    }

    spinor_sum operator*(const double &c, const spinor_sum &sum){
        spinor_sum result = sum;
        for (int i = 0; i < sum.terms.size(); ++i){
            result.terms[i].factor *= c;
        }
        return result;
    }
    spinor_sum operator*(const spinor_sum &sum,const double &c){
        spinor_sum result = sum;
        for (int i = 0; i < sum.terms.size(); ++i){
            result.terms[i].factor *= c;
        }
        return result;
    }

    spinor_fraction operator*(const spinor_fraction &f1, const spinor &s2){
        spinor_fraction f2(s2);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor &s1, const spinor_fraction &f2){
        spinor_fraction f1(s1);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_product &p2){
        spinor_fraction f2(p2);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor_product &p1, const spinor_fraction &f2){
        spinor_fraction f1(p1);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_sum &sum2){
        spinor_fraction f2(sum2);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor_sum &sum1, const spinor_fraction &f2){
        spinor_fraction f1(sum1);
        spinor_fraction result = f1 * f2;
        return result;
    }
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_fraction &f2){
        spinor_fraction result;
        if (f2.denominator.terms.empty()){
            result.numerator = f1.numerator * f2.numerator;
            result.denominator = f1.denominator;
        }
        else if (f2.numerator.terms.empty()){
            result.numerator = f1.numerator;
            result.denominator = f1.denominator * f2.denominator;
        }
        else{
            result.numerator = f1.numerator * f2.numerator;
            result.denominator = f1.denominator * f2.denominator;
        }

        return result;
    }
    spinor_fraction operator*(const double &c, const spinor_fraction &f){
        spinor_fraction result = f;
        result.numerator = result.numerator * c;
        return result;
    }
    spinor_fraction operator*(const spinor_fraction &f, const double &c){
        spinor_fraction result = f;
        result.numerator = result.numerator * c;
        return result;
    }

// ---------- Division ----------

    spinor_product operator/(const spinor &s1, const spinor &s2){
        spinor_product p1(s1);
        spinor_product p2(s2);
        spinor_product result = p1 / p2;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator/(const spinor_product &p1, const spinor &s2){
        spinor_product p2(s2);
        spinor_product result = p1 / p2;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator/(const spinor &s1, const spinor_product &p2){
        spinor_product p1(s1);
        spinor_product result = p1 / p2;
        order(result);
        reduce_product(result);
        return result;
    }
    spinor_product operator/(const spinor_product &p1, const spinor_product &p2){
        spinor_product result = p1;
        result.numerator.insert(result.numerator.end(),p2.denominator.begin(),p2.denominator.end());
        result.denominator.insert(result.denominator.end(),p2.numerator.begin(),p2.numerator.end());
        result.factor /= p2.factor;
        order(result);
        reduce_product(result);
        return result;
    }

    spinor_product operator/(const double &c, const spinor &s){
        spinor_product result;
        result.denominator.push_back(s);
        result.factor *= c;
        return result;
    }

    spinor_product operator/(const double &c, const spinor_product &p){
        spinor_product result;
        result.denominator = p.numerator;
        result.numerator = p.denominator;
        result.factor *= c;
        return result;
    }

    spinor operator/(const spinor &s,const double &c){
        if (std::abs(c) < 1e-12) {
            throw std::runtime_error("Division by zero scalar");
        }

        spinor result = s;
        result.factor /= c;
        return result;
    }

    spinor_product operator/(const spinor_product &p,const double &c){
        if (std::abs(c) < 1e-12) {
            throw std::runtime_error("Division by zero scalar");
        }

        spinor_product result = p;
        result.factor /= c;
        return result;
    }
    spinor_sum operator/(const spinor_sum &sum1,const spinor &s2){
        spinor_product p2(s2);
        spinor_sum result;
        for (int i = 0; i < sum1.terms.size(); ++i){
            result.terms.push_back(sum1.terms[i]/p2);
        }
        return result;
    }

    spinor_fraction operator/(const spinor_product &p1,const spinor_sum &sum2){
        spinor_sum sum1(p1);
        spinor_fraction result(sum1,sum2);
        reduce_fraction(result);
        return result;
    }
    spinor_sum operator/(const spinor_sum &sum1,const spinor_product &p2){
        spinor_sum result;
        for (int i = 0; i < sum1.terms.size(); ++i){
            result.terms.push_back(sum1.terms[i]/p2);
        }
        return result;
    }
    spinor_fraction operator/(const spinor_sum &sum1,const spinor_sum &sum2){
        spinor_fraction result(sum1,sum2);
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor &s1,const spinor_sum &sum2){
        spinor_sum sum1(s1);
        spinor_fraction result(sum1,sum2);
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const double &c, const spinor_sum &sum){
        spinor_fraction result;
        result.denominator = sum;
        for (int i = 0; i < sum.terms.size(); ++i){
            result.denominator.terms[i].factor /= c;
        }
        reduce_fraction(result);
        return result;
    }
    spinor_sum operator/(const spinor_sum &sum,const double &c){    
        if (std::abs(c) < 1e-12) {
            throw std::runtime_error("Division by zero scalar");
        }

        spinor_sum result = sum;
        for (int i = 0; i < sum.terms.size(); ++i){
            result.terms[i].factor /= c;
        }
        return result;
    }
    spinor_fraction operator/(const spinor_fraction &f1, const spinor &s2){
        spinor_fraction f2(s2);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor &s1, const spinor_fraction &f2){
        spinor_fraction f1(s1);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_product &p2){
        spinor_fraction f2(p2);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_product &p1, const spinor_fraction &f2){
        spinor_fraction f1(p1);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_sum &sum2){
        spinor_fraction f2(sum2);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_sum &sum1, const spinor_fraction &f2){
        spinor_fraction f1(sum1);
        spinor_fraction result = f1 / f2;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_fraction &f2){
        spinor_fraction result;
        if (f2.denominator.terms.empty()){
            result.numerator = f1.numerator;
            result.denominator = f1.denominator * f2.numerator;
        }
        else if(f2.numerator.terms.empty()){
            result.numerator = f1.numerator * f2.denominator;
            result.denominator = f1.denominator;
        }
        else{
            result.numerator = f1.numerator * f2.denominator;
            result.denominator = f1.denominator * f2.numerator;
        }
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const double &c, const spinor_fraction &f){
        if (f.numerator.terms.empty()) {
            throw std::runtime_error("Division by zero fraction");
        }

        spinor_fraction result;
        result.numerator = f.denominator;
        result.denominator = f.numerator;
        result.numerator = c * result.numerator;
        reduce_fraction(result);
        return result;
    }
    spinor_fraction operator/(const spinor_fraction &f, const double &c){
        if (std::abs(c) < 1e-12) {
            throw std::runtime_error("Division by zero scalar");
        }

        spinor_fraction result = f;
        result.numerator = result.numerator / c;
        reduce_fraction(result);
        return result;
    }

// ---------- Order ----------

    bool operator<(const spinor& a, const spinor& b){
        if (a.i < b.i){return true;}
        else if (a.i > b.i){return false;}
        else{
            if (a.j < b.j){return true;}
            else {return false;}
        }
    }

    bool operator<(const spinor_product& p1, const spinor_product& p2) {
        std::vector<spinor> list1 = p1.numerator;
        std::vector<spinor> list2 = p2.numerator;
        list1.insert(list1.end(),p1.denominator.begin(),p1.denominator.end());
        list2.insert(list2.end(),p2.denominator.begin(),p2.denominator.end());

        return std::lexicographical_compare(list1.begin(),list1.end(),list2.begin(),list2.end());
    }

    // Sort by descending factor_count (higher factor_count comes first in std::sort)
    bool operator<(const factor_struct &fs1, const factor_struct &fs2){
        return fs1.factor_count > fs2.factor_count;
    }

    bool operator==(const spinor& s1, const spinor& s2) {
        return s1.factor == s2.factor && s1.bracket == s2.bracket && s1.i == s2.i && s1.j == s2.j;
    }

    bool operator!=(const spinor& s1, const spinor& s2) {
        return !(s1 == s2);
    }

    bool operator==(const std::vector<spinor> &vec1, const std::vector<spinor> &vec2){
        auto temp1 = vec1;
        auto temp2 = vec2;
        if (temp1.size() == temp1.size()){
            std::sort(temp1.begin(), temp1.end());
            std::sort(temp2.begin(), temp2.end());
            
            for (int i = 0; i < temp1.size(); ++i){
                if(!is_similar_spinor(temp1[i],temp2[i])){
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    bool operator==(const spinor_product& p1, const spinor_product& p2) {
        if(p1.factor == p2.factor and p1.numerator == p2.numerator and p1.denominator == p2.denominator){
            return true;
        }
        return false;
    }

    bool operator==(const spinor_sum& sum1, const spinor_sum& sum2) {
        if (sum1.terms.size() != sum2.terms.size()){
            return false;
        }

        auto sort1 = sum1.terms;
        auto sort2 = sum2.terms;

        std::sort(sort1.begin(),sort1.end());
        std::sort(sort2.begin(),sort2.end());

        return std::equal(sort1.begin(), sort1.end(), sort2.begin());
    }

    bool operator!=(const spinor_product& p1, const spinor_product& p2) {
        return !(p1 == p2);
    }

    bool operator!=(const spinor_sum& sum1, const spinor_sum& sum2) {
        return !(sum1 == sum2);
    }

    void order(spinor_product &p){
        std::sort(p.numerator.begin(), p.numerator.end());
        std::sort(p.denominator.begin(), p.denominator.end());
    }

    bool is_similar_spinor(const spinor s1, const spinor s2){
         return s1.bracket == s2.bracket and s1.i == s2.i and s1.j == s2.j;
    }

    bool is_similar_spinor_vector(const std::vector<spinor> &vec1, const std::vector<spinor> &vec2){
        auto temp1 = vec1;
        auto temp2 = vec2;
        if (temp1.size() == temp1.size()){
            if (temp1.empty() or temp2.empty()){
                return true;
            }

            std::sort(temp1.begin(), temp1.end());
            std::sort(temp2.begin(), temp2.end());
            
            for (int i = 0; i < temp1.size(); ++i){
                if(!is_similar_spinor(temp1[i],temp2[i])){
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    bool is_similar(const spinor_product p1, const spinor_product p2){
        if(is_similar_spinor_vector(p1.numerator,p2.numerator) and is_similar_spinor_vector(p1.denominator,p2.denominator)){
            return true;
        }
        else {return false;}
    }

// ---------- Simplify ----------

    void cancellation(spinor_sum &sum){
        std::vector<spinor_product> result;
        std::vector<bool> used(sum.terms.size(), false);
        if (sum.terms.size() <= 1){
            return;
        }

        for (int i = 0; i < sum.terms.size(); ++i){
            if (used[i]) {continue;}

            for (int j = i + 1; j < sum.terms.size(); ++j){
                if (used[j]) {continue;}

                if (is_similar(sum.terms[i],sum.terms[j])){
                    if (std::abs(sum.terms[i].factor + sum.terms[j].factor) < 1e-9){
                        used[i] = true;
                        used[j] = true;
                    }
                    else{
                        used[j] = true;
                        sum.terms[i].factor = sum.terms[i].factor + sum.terms[j].factor;
                    }
                }
            }
            if (!used[i]){
                result.push_back(sum.terms[i]);
            }
        }
        sum.terms = result;
        if (sum.terms.empty()){
            sum.is_null = true;
        }
    }

    void reduce_product(spinor_product &p){
        spinor_product temp;
        std::vector<int> indices{};
        
        if (!p.denominator.empty()){
            for (spinor &s : p.denominator){
                auto it = std::find(p.numerator.begin(),p.numerator.end(),s);

                if (it != p.numerator.end()){
                    indices.push_back(distance(p.numerator.begin(),it));
                }
                else {temp.denominator.push_back(s);}

                s.factor = 1;
            }

            if (!indices.empty()){
                for (int j = indices.size() - 1; j >= 0; --j){
                    p.numerator.erase(p.numerator.begin() + indices[j]);
                }
            }

            p.denominator = temp.denominator;
        }
        if(!p.numerator.empty()){
            for(spinor &s : p.numerator){
                s.factor = 1;
            }
        }

    }
// ---------- Factorization (only works for linear expressions with integer coefficients) ----------

    bool is_spinor_in_numerator(const spinor_product &p, const spinor &s){
        if(std::find(p.numerator.begin(),p.numerator.end(),s) != p.numerator.end()){
            return true;
        }
        return false;
    }

    std::vector<factor_struct> find_common_spinors(const spinor_sum &sum){
        std::vector<factor_struct> result;
        std::vector<spinor> all_terms;
        spinor_product current_product;

        // Appending all terms with no duplicates in specific product 
        for (const spinor_product &p : sum.terms){
            current_product = p;
            std::sort(current_product.numerator.begin(),current_product.numerator.end());
            auto it = std::unique(current_product.numerator.begin(),current_product.numerator.end());
            current_product.numerator.erase(it,current_product.numerator.end());

            all_terms.insert(all_terms.end(),current_product.numerator.begin(),current_product.numerator.end());
        }

        // Sorting vector to count duplicates accross products (most likely factors)
        std::sort(all_terms.begin(),all_terms.end());

        spinor current_spinor = all_terms[0];
        factor_struct current_factor;
        int count = 1;

        for (int i = 1; i < all_terms.size(); ++i) {
            if (is_similar(current_spinor,all_terms[i])) {
                count += 1;
            } 
            else {
                if (count > 1) {
                    current_spinor.factor = 1;
                    current_factor = factor_struct(current_spinor, count);
                    result.push_back(current_factor);
                }
            current_spinor = all_terms[i];
            count = 1;
            }
        }

        if (count > 1) {
            current_spinor.factor = 1;
            current_factor = factor_struct(current_spinor, count);
            result.push_back(current_factor);
        }
        // Sorting the product again to have highest count first
        std::sort(result.begin(),result.end());
        return result;
    }
    
    int gcd(int a, int b) {
        while (b != 0) {
            int tmp = b;
            b = a % b;
            a = tmp;
            }
        return a;
    }

    std::vector<int> find_unique_gcds(const spinor_sum &sum){
        std::vector<int> result;

        for(int i = 0; i < sum.terms.size(); ++i){
            double gcd_value;
            for (int j = i + 1; j < sum.terms.size(); ++j){
                gcd_value = gcd(int(sum.terms[i].factor),int(sum.terms[j].factor));

                if (std::find(result.begin(),result.end(),gcd_value) == result.end()){
                    result.push_back(gcd_value);
                }
            }
        }
        return result;
    }

    std::set<std::vector<int>> make_gcd_combinations(const int &gcd_number){
        std::vector<int> arrangement;
        for(int i = 0; i < gcd_number; ++i){
            arrangement.push_back(i);
        }
        
        std::set<std::vector<int>> permutations;
        do {
            permutations.insert(arrangement);
        } 
        while (std::next_permutation(arrangement.begin(), arrangement.end()));
        return permutations;
    }

    bool is_perfect_product(const spinor_sum &sum, const int &highest_order, spinor_sum &factor,spinor_sum &leftover_sum){
        spinor_sum perfect_product; 
        std::vector<int> coefficients;

        for(spinor_product f : factor.terms){
            spinor_sum temp = sum;
            for (int i = 0; i < highest_order; ++i){
                temp = temp/f.numerator[0];
            }

            for(spinor_product term : temp.terms){
                if (term.denominator.empty()){
                    coefficients.push_back(int(std::exp(std::log(term.factor)/highest_order)));
                }
            }
        }
       
        spinor_sum new_factor;
        int i = 0;
        for(spinor_product f : factor.terms){
            new_factor = new_factor + f.numerator[0] * coefficients[i];
            ++i;
        }

        for (int i = 0; i < highest_order; ++i){
            perfect_product = perfect_product * new_factor;
        }
        
        spinor_sum result = sum - perfect_product;
        if (result.is_null){
            for (int i = 0; i < highest_order - 1; ++i){
                leftover_sum = leftover_sum * new_factor;
            }
            factor = new_factor;
            return true;
        }
        return false;
    }

    bool is_factor(Polynomial poly, const spinor_sum & sum, const spinor_sum &factor, spinor_sum &leftover){        
        spinor_sum leftover_sum;

        int i = 0;
        while(poly.highest_order > i){
            leftover_sum = leftover_sum + (poly.terms[poly.highest_order - i]/factor.terms[0]);
            poly = poly - ((poly.terms[poly.highest_order - i]/factor.terms[0]) * factor);

            if(poly.is_null()){
                leftover = leftover_sum;
                return true;
            }
            i++;
        }
        return false;
    }

    factorized_spinor factor(factorized_spinor &sum, std::vector<factor_struct> possible_factors){
        if (possible_factors.empty() or sum.remainder.terms.size() <= 2) {
            return sum;
        }

        factor_struct current_factor = possible_factors[0];
        std::vector<int> gcd_values = find_unique_gcds(sum.remainder);
        std::set<std::vector<int>> permutations = make_gcd_combinations(gcd_values.size());

        // Estimate amount of possible factors (I assume about half)
        int N = 1 << int(possible_factors.size()/2 + 1);

        for (int mask = 1; mask < N; ++mask) {
            spinor_sum candidate_factor(current_factor.factor);

            for (int i = 1; i < possible_factors.size(); ++i) {
                if(mask & (1 << i)){
                   candidate_factor = candidate_factor + possible_factors[i].factor;
                }
            }

            // Make the current entry into polynomial form to do long division
            spinor var = candidate_factor.terms[0].numerator[0];
            var.factor = candidate_factor.terms[0].factor;
            Polynomial poly(sum.remainder,var);
            
            // Try factoring into perfect square
            spinor_sum leftover;
            spinor_sum candidate = candidate_factor;
            if(is_perfect_product(sum.remainder,poly.highest_order,candidate,leftover)){      
                sum.factor.push_back(candidate);
                sum.remainder = leftover;
                return factor(sum,possible_factors);
            }

            // Try all permutations of GCD assignments
            for (std::vector<int> perms : permutations){
                spinor_sum permuted_candidate = candidate_factor;

                for(size_t j = 0; j < permuted_candidate.terms.size(); ++j){
                    permuted_candidate.terms[j].factor = gcd_values[perms[j]];
                }

                if (is_factor(poly, sum.remainder, permuted_candidate, leftover)){
                    sum.factor.push_back(permuted_candidate);
                    sum.remainder = leftover;
                    return factor(sum,possible_factors);
                }
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

    void flatten_fraction(spinor_fraction &f){
        if (!f.numerator.terms.empty()){
            for (spinor_product &p : f.numerator.terms){
                if (p.denominator.size() != 0){
                    for (spinor elem : p.denominator){
                        f.numerator = elem * f.numerator ;
                        f.denominator = elem * f.denominator;
                    }
                    p.denominator.clear();
                }
            }
        cancellation(f.numerator);
        }

        if (!f.denominator.terms.empty()){
            for (spinor_product &p : f.denominator.terms){
                if (p.denominator.size() != 0){
                    for (spinor elem : p.denominator){
                        f.numerator = elem * f.numerator ;
                        f.denominator = elem * f.denominator;
                    }
                p.denominator.clear();
                }
            }
        cancellation(f.denominator);
        }
    }


    void reduce_fraction(spinor_fraction &f){
        if (f.denominator.terms.size() == 0 or f.numerator.terms.size() == 0){
            return;
        }

        flatten_fraction(f);

        if (f.denominator.terms.size() == 1){
            spinor_sum temp = f.numerator; 
            temp = temp / f.denominator.terms[0];
            for (auto &t : temp.terms){
                reduce_product(t);
            }
            spinor_fraction result(temp);
            f = result;
        }
        else{
            std::vector<factor_struct> cand_num = find_common_spinors(f.numerator);
            std::vector<factor_struct> cand_denom = find_common_spinors(f.denominator);
            factorized_spinor num(f.numerator);
            factorized_spinor denom(f.denominator);

            factorized_spinor factored_num = factor(num,cand_num);   
            factorized_spinor factored_denom = factor(denom,cand_denom);

            factored_num.factor.push_back(factored_num.remainder);
            factored_num.remainder = {};

            factored_denom.factor.push_back(factored_denom.remainder);
            factored_denom.remainder = {};

            spinor_fraction result;

            std::vector<spinor_sum> remaining_denom_factors = factored_denom.factor;
            
            for (const auto& fact : factored_num.factor) {
                auto it = std::find(remaining_denom_factors.begin(), remaining_denom_factors.end(), fact);
                if (it != remaining_denom_factors.end()) {
                    remaining_denom_factors.erase(it);
                } 
                else {
                    result.numerator = result.numerator * fact;
                }
            }
            
            for (const auto& fact : remaining_denom_factors) {
                result.denominator = result.denominator * fact;
            }
            f = result;
        }
    }


// ---------- Display ----------

    std::string display_spinor(const spinor& s) {
        if (s.bracket == 'a') {
            return "(" + std::to_string(s.i) + "," + std::to_string(s.j) + ")";
        } 
        else {
            return "[" + std::to_string(s.i) + "," + std::to_string(s.j) + "]";
        }
    }

    std::ostream& operator<<(std::ostream& os, const spinor& s) {
        if (s.factor < 0){
            os << " - ";
        }

        if (std::abs(s.factor) != 1){
            os << std::to_string(std::abs(s.factor));
        }

        os << display_spinor(s);

        return os;
    }

    std::ostream& operator<<(std::ostream& os, const spinor_product& p) {
        if (p.factor < 0){
            os << " - ";
        }

        if (std::abs(p.factor) != 1){
            os << std::to_string(std::abs(p.factor)) + " ";
        }

        if (!p.numerator.empty()){
            for (int i = 0; i < p.numerator.size(); ++i){
                os << display_spinor(p.numerator[i]);
            }

            if (!p.denominator.empty()){
                os << " / ( ";

                for (int i = 0; i < p.denominator.size(); ++i){
                    os << display_spinor(p.denominator[i]);
                }
                os << " )";
            }
            return os;
        }

        else{
            if (!p.denominator.empty()){
                if (p.factor == 1){
                    os << "1 ";
                }

                os << "/ ( ";

                for (int i = 0; i < p.denominator.size(); ++i){
                    os << display_spinor(p.denominator[i]);
                }
                os << " )";
            }
            return os;
        }
    }

    std::ostream& operator<<(std::ostream& os,const spinor_sum& sum) {
        if (sum.terms.empty()){
            os << "0";
        }
        else{
            for (int i = 0; i < sum.terms.size(); ++i){
                if (i > 0 and sum.terms[i].factor > 0){
                    os << " + ";
                }
                os << sum.terms[i];
            }
        }
        return os;
    }

    std::ostream& operator<<(std::ostream& os,const spinor_fraction& f) {
        if (f.numerator.terms.empty() and f.denominator.terms.empty()){
            if (f.numerator.is_null){
                os << "0";
            }
            else{
                os << "1";
            }
                os << "1";
        }
        else if(f.numerator.terms.empty()){
            os << "1 / ";
            os << f.denominator;
        }
        else if(f.denominator.terms.empty()){
            os << f.numerator;
        }
        else{
            os << f.numerator;
            os << " / (";
            os << f.denominator;
            os << " )";
        }

        return os;
    }

    std::ostream& operator<<(std::ostream& os,const factorized_spinor& f) {
        if (f.factor.empty()){
            os << "factor: 1 ,";
        }
        else {
            for (spinor_sum facs : f.factor){
                os << facs << " , ";
            }
        }
        os << "Remainder: " << f.remainder;
        return os;
    }

    std::vector<std::string> tokenize(const std::string &input){
        std::vector<std::string> tokens;
        std::string current_string;

        for (char c : input){  
            if (std::isspace(c)) continue;

            if (std::isdigit(c) or c == ',' or c == '[' or c == ']') {
                current_string += c;
            }
            else if (c == '(' ){
                if (!current_string.empty()){
                    tokens.push_back(current_string);
                    current_string.clear();
                }
                current_string += c;
            } 
            else if (c == ')' ){
                if (!current_string.empty() and current_string.back() == ')'){
                    tokens.push_back(current_string);
                    current_string.clear();
                }
                current_string += c;
            } 

            else{
                if (!current_string.empty()) {
                    tokens.push_back(current_string);
                    current_string.clear();
                }

                if (c == '*' or c == '/' or c == '+' or c == '-' or '.') {
                    tokens.emplace_back(1, c);
                }
            }
        }
        tokens.push_back(current_string);
        return tokens;
    }

    spinor parse_spinor(const std::string &token){
        spinor s;

        if (token.front() == '('){
            s.bracket = 'a';
        }
        else {s.bracket = 'b';}

        auto it = std::find(token.begin(),token.end(),',');

        int pos = std::distance(token.begin(), it);
    
        s.i = std::stoi(token.substr(1,pos - 1));
        s.j = std::stoi(token.substr(pos + 1, token.length() - pos - 2));

        return s;
    }

}

 int main(){
    using namespace Spinor;
    spinor spinor1(2,5,'a');
    spinor spinor2(2,6,'a');
    spinor spinor3(2,7,'b');
    spinor spinor4(1,2,'a');
    spinor spinor5(3,4,'a');
    spinor spinor6(1,4,'b');
    auto new_spinor = ((2*spinor1+3*spinor2)*(2*spinor1+3*spinor2)*(2*spinor1+3*spinor2))/((2*spinor1 + 3*spinor2)*(spinor3+4*spinor4));
    //((2*spinor2 + spinor3) * (3*spinor1 + spinor4) )/(2 * spinor2 + spinor3);
    //((2*spinor2*spinor4 + 3*spinor3 * spinor1 + spinor3 * spinor4+ 6*spinor2*spinor1))/(2 * spinor2 + spinor3);
    std::string line = "(1,34)";

    auto outputs = tokenize(line);

    std::cout << new_spinor << std::endl;

    //for (int i = 0; i < outputs.size(); ++i){
    //    std::cout << outputs[i] << " ";
    //}

    //spinor spinor5 = parse_spinor(line);

    //std::cout << spinor5 << std::endl;
    //factorized_spinor spinf(new_spinor);

    //std::vector<factor_struct> factors = find_common_spinors(new_spinor);
    // auto factored_spinor = factor(spinf,factors);
    //std::cout << factored_spinor << std::endl;



    //for(factor_struct f :factors){
    //    std::cout << f.factor << ", count: " << f.factor_count << std::endl;
    //}




    //std::cout << new_spinor.terms[0] << std::endl;

    return 0;
    }

